function aj_automate_spm_results(GLM_dir)
    % Démarrer SPM12
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    % Définir les sous-dossiers et les contrastes associés
    subfolders = {'MTsat_GM', 'MTsat_WM', 'PDmap_GM', 'PDmap_WM', 'R1map_GM', 'R1map_WM', 'R2starmap_GM', 'R2starmap_WM'};
    contrast_type = containers.Map;
    contrast_type('MTsat_') = 'Decrease';
    contrast_type('PDmap_') = 'Increase';
    contrast_type('R1map_') = 'Decrease';
    contrast_type('R2starmap_') = 'Increase';

    % Parcourir chaque sous-dossier
    for i = 1:length(subfolders)
        subfolder = subfolders{i};
        spm_mat_path = fullfile(GLM_dir, subfolder, 'SPM.mat');

        if ~exist(spm_mat_path, 'file')
            warning('SPM.mat non trouvé dans %s', fullfile(input_dir, subfolder));
            continue;
        end

        % Charger SPM.mat
        loaded_data = load(spm_mat_path);
        SPM = loaded_data.SPM;

        % Mettre à jour SPM.swd avec le bon chemin
        SPM.swd = fileparts(spm_mat_path);

        % Déterminer le type de contraste à utiliser
        contrast_key = [subfolder(1:find(subfolder=='_', 1)-1), '_'];
        contrast_type_name = contrast_type(contrast_key);

        % Trouver l'index du contraste
        contrast_idx = [];
        for j = 1:length(SPM.xCon)
            if contains(SPM.xCon(j).name, contrast_type_name)
                contrast_idx = j;
                break;
            end
        end

        if isempty(contrast_idx)
            % Afficher les noms des contrastes disponibles
            fprintf('Contrast names in %s:\n', spm_mat_path);
            for j = 1:length(SPM.xCon)
                fprintf('%d: %s\n', j, SPM.xCon(j).name);
            end
            warning('Contraste contenant "%s" non trouvé dans %s', contrast_type_name, spm_mat_path);
            continue;
        end

        % Déterminer si c'est AR+ ou AR-
        if contains(contrast_type_name, 'Increase')
            ar_sign = 'AR+';
        else
            ar_sign = 'AR-';
        end

        % Configurer les paramètres de visualisation
        xSPM = struct();
        xSPM.swd = SPM.swd;
        xSPM.Ic = contrast_idx;
        xSPM.Im = [];
        xSPM.Ex = [];
        xSPM.thresDesc = 'FWE';
        xSPM.u = 0.05;
        xSPM.k = 0;

        % Obtenir les résultats
        [hReg, xSPM] = spm_getSPM(xSPM);
        if isempty(hReg)
            warning('Impossible de récupérer les résultats pour %s', spm_mat_path);
            continue;
        end

        % Lire les données de la carte statistique
        Z = spm_read_vols(xSPM.Vspm);

        % Calculer le seuil de manière alternative
        try
            u = spm_u(xSPM.u, xSPM.df, xSPM.STAT);
        catch
            % Si spm_u échoue, utiliser une autre méthode pour obtenir le seuil
            if strcmp(xSPM.STAT, 'T')
                u = spm_invCdft(xSPM.df, 1 - xSPM.u/2); % seuil bilatéral pour T
            else
                u = spm_invCdff(xSPM.df(1), xSPM.df(2), 1 - xSPM.u); % seuil pour F
            end
        end
        
        % Créer un masque binaire des clusters significatifs
        binary_mask = Z > u;

        % Sauvegarder le masque binaire
        mask_name = sprintf('mask_FWE005_%s_%s.nii', subfolder, ar_sign);
        V = struct( ...
            'fname', fullfile(GLM_dir, subfolder, mask_name), ...
            'dim', xSPM.DIM, ...
            'dt', [spm_type('uint8') spm_platform('bigend')], ...
            'mat', xSPM.M, ...
            'pinfo', [1; 0; 0], ...
            'desc', 'binary mask' ...
        );
        V = spm_create_vol(V);
        V = spm_write_vol(V, uint8(binary_mask));

        % Créer une figure et enregistrer l'image
        fig = figure;
        spm_orthviews('addcolouredimage', 1, Z, xSPM.XYZmm, xSPM.M);
        spm_orthviews('redraw');
        fig_name = sprintf('%s_%s_FWE005.png', subfolder, ar_sign);
        saveas(fig, fullfile(GLM_dir, subfolder, fig_name));

        % Fermer la fenêtre des résultats
        close all;
    end
end