% Main script to compute the brain 1D profile analysis. 
%
% First part of the script. The plots contain
% the unsmoothed and smoothed using TWS and TSPOON signals extracted along
% one selected axis in the quantitative 3D maps. The tissue density and TWS
% and TSPOON denominator signals are also extracted as the WTA masks.
%
% Second part of the script: The plots contain the WTA and merged
% unsmoothed and smoothed using TWS and TSPOON signals. Merging is done
% through GM and WM by applying the WTA masks to GM and WM signals.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Save plot after managing it manually (full screen)
dir_out = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\qMRI_results\Brain1D';
pth_out = fullfile(dir_out, sprintf('%s_%s.%s.%s_S%d_HR.eps',qMRI_param,coord{1},coord{2},coord{3},subject));
exportgraphics(figureHandle, pth_out, 'Resolution', 450); % file format defined in pth_out

%% Plot unsmoothed + TWS/TSPOON GM/WM signals & TWS/TSPOON GM/WM denominator signals & GM/WM population masks
close all; clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
TCs = {'GM', 'WM', 'CSF'};

% Choose the one to plot
qMRI_params = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
qMRI_param_num = 4;
qMRI_param = qMRI_params{qMRI_param_num};
subject = 1; % (1:138)
coord = {'100', '', '100'}; % [181 217 181]
text_size = 12;

% Get the subject & qMRI parameter -specific MNI path
MNI_dir = fullfile(ds_dir,'derivatives','SPM12_dartel');
MNI_sub_dir = fullfile(MNI_dir, sprintf('sub-%03d', subject), 'anat');
MNI_pattern = sprintf('^.*space-MNI_%s\\.nii$', qMRI_param);
MNI_path = spm_select('ExtFPList', MNI_sub_dir, MNI_pattern);

% Get the subject-specific GM/WM density path
GM_density_pattern = sprintf('^.*label-%s_probseg\\.nii$', TCs{1});
WM_density_pattern = sprintf('^.*label-%s_probseg\\.nii$', TCs{2});
GM_density_path = spm_select('ExtFPList', MNI_sub_dir, GM_density_pattern);
WM_density_path = spm_select('ExtFPList', MNI_sub_dir, WM_density_pattern);

% Get the subject & qMRI parameter -specific GM/WM TWS path
TWS_dir = fullfile(ds_dir,'derivatives','AJ-TWS');
TWS_sub_dir = fullfile(TWS_dir, sprintf('sub-%03d', subject), 'anat');
TWS_GM_pattern = sprintf('^TWS_%s.*%s\\.nii$', TCs{1}, qMRI_param);
TWS_WM_pattern = sprintf('^TWS_%s.*%s\\.nii$', TCs{2}, qMRI_param);
TWS_GM_path = spm_select('ExtFPList', TWS_sub_dir, TWS_GM_pattern);
TWS_WM_path = spm_select('ExtFPList', TWS_sub_dir, TWS_WM_pattern);

% Get the subject-specific smoothed GM/WM density path (non thresholded TWS denominator)
sGM_pattern = sprintf('^gs.*label-%s_probseg\\.nii$', TCs{1});
sWM_pattern = sprintf('^gs.*label-%s_probseg\\.nii$', TCs{2});
sGM_path = spm_select('ExtFPList', TWS_sub_dir, sGM_pattern);
sWM_path = spm_select('ExtFPList', TWS_sub_dir, sWM_pattern);

% Get the subject & qMRI parameter -specific GM/WM TSPOON path
TSPOON_dir = fullfile(ds_dir,'derivatives','AJ-TSPOON');
TSPOON_sub_dir = fullfile(TSPOON_dir, sprintf('sub-%03d', subject), 'anat');
TSPOON_GM_pattern = sprintf('^TSPOON_%s.*%s\\.nii$', TCs{1}, qMRI_param);
TSPOON_WM_pattern = sprintf('^TSPOON_%s.*%s\\.nii$', TCs{2}, qMRI_param);
TSPOON_GM_path = spm_select('ExtFPList', TSPOON_sub_dir, TSPOON_GM_pattern);
TSPOON_WM_path = spm_select('ExtFPList', TSPOON_sub_dir, TSPOON_WM_pattern);

% Get the subject-specific smoothed GM/WM masks path (non thresholded TSPOON denominator)
GM_sMask_pattern = sprintf('^gs_Mask_gs.*label-%s_probseg\\.nii$', TCs{1});
WM_sMask_pattern = sprintf('^gs_Mask_gs.*label-%s_probseg\\.nii$', TCs{2});
GM_sMask_path = spm_select('ExtFPList', TSPOON_sub_dir, GM_sMask_pattern);
WM_sMask_path = spm_select('ExtFPList', TSPOON_sub_dir, WM_sMask_pattern);

if isempty(GM_sMask_path) || isempty(WM_sMask_path)
    gs_TCmap_paths = spm_select('ExtFPList', TSPOON_sub_dir, '^gs.*label-.*_probseg\.nii$');
    sMasks_paths = aj_get_sMask(cellstr(gs_TCmap_paths)', TSPOON_sub_dir);
    GM_sMask_path = sMasks_paths{2}{2};
    WM_sMask_path = sMasks_paths{2}{3};
    file2delete = char(sMasks_paths{1}{:},sMasks_paths{2}{1},sMasks_paths{3}{:});
    for i = 1:text_size(file2delete,1)
        delete(file2delete(i,:));
    end
end

% Get the population GM/WM masks (winner takes all)
deriv_dir = fullfile(ds_dir,'derivatives');
WTAM_paths = cell(1,2);
WTAM_paths{1} = spm_select('ExtFPList', deriv_dir, '^atlas-GM_space-MNI_mask\.nii$');
WTAM_paths{2} = spm_select('ExtFPList', deriv_dir, '^atlas-WM_space-MNI_mask\.nii$');

% Get 1D profiles for the first subplot
subplot1_paths = cell(1,5);
subplot1_paths{1} = MNI_path;
subplot1_paths{2} = TWS_GM_path;
subplot1_paths{3} = TWS_WM_path;
subplot1_paths{4} = TSPOON_GM_path;
subplot1_paths{5} = TSPOON_WM_path;

% Get 1D profiles for the second subplot
subplot2_paths = cell(1,6);
subplot2_paths{1} = GM_density_path;
subplot2_paths{2} = WM_density_path;
subplot2_paths{3} = sGM_path;
subplot2_paths{4} = sWM_path;
subplot2_paths{5} = GM_sMask_path;
subplot2_paths{6} = WM_sMask_path;

% Obtenir les profils 1D pour le premier sous-graphe
[subplot1_profile, subplot1_coords] = aj_get_profile(subplot1_paths, coord);

% Obtenir les profils 1D pour le deuxième sous-graphe
[subplot2_profile, subplot2_coords] = aj_get_profile(subplot2_paths, coord);

% Obtenir les profils 1D pour le troisième sous-graphe
[subplot3_profile, subplot3_coords] = aj_get_profile(WTAM_paths, coord);

% Déterminer les coordonnées MNI utilisées
mni_coords = subplot1_coords{1}; % Supposé identique pour toutes les données extraites

% Group profiles (inchangé)
x = mni_coords; % Utiliser les coordonnées MNI au lieu des indices
y1 = [subplot1_profile{1}(:), subplot1_profile{2}(:), subplot1_profile{4}(:)];
y2 = [subplot2_profile{1}(:), subplot2_profile{3}(:), subplot2_profile{5}(:)];
y3 = [subplot3_profile{1}(:), subplot3_profile{2}(:)];
y4 = [subplot2_profile{2}(:), subplot2_profile{4}(:), subplot2_profile{6}(:)];
y5 = [subplot1_profile{1}(:), subplot1_profile{3}(:), subplot1_profile{5}(:)];

% titles = {'GM Signal Intensity', 'GM Class', 'WTA Masks', 'WM Class', 'WM Signal Intensity'};
if qMRI_param_num == 1
    qMRI = 'MTsat';
    ylabels = {'Intensity [AU]', 'Intensity [/]', 'Binary Mask Value', 'Intensity [/]', 'Intensity [AU]'};
elseif qMRI_param_num == 2
    qMRI = 'PD';
    ylabels = {'Intensity [%]', 'Intensity [/]', 'Binary Mask Value', 'Intensity [/]', 'Intensity[ %]'};
elseif qMRI_param_num == 3
    qMRI = 'R1';
    ylabels = {'Intensity [s^{-1}]', 'Intensity [/]', 'Binary Mask Value', 'Intensity [/]', 'Intensity [s^{-1}]'};
elseif qMRI_param_num == 4
    qMRI = 'R2*';
    ylabels = {'Intensity [s^{-1}]', 'Intensity [/]', 'Binary Mask Value', 'Intensity [/]', 'Intensity [s^{-1}]'};
end

legends = cell(1,5);
legends{1} = {'Not Smoothed', 'TWS GM', 'TSPOON GM'};
legends{2} = {'GM Density', 'TWS sGM', 'TSPOON sMask GM'};
legends{3} = {'GM', 'WM'};
legends{4} = {'WM Density', 'TWS sWM', 'TSPOON sMask WM'};
legends{5} = {'Not Smoothed', 'TWS WM', 'TSPOON WM'};

% Transition detection
transition_indices_top = []; % Transitions from the first curve of y3
transition_indices_bottom = []; % Transitions from the second curve of y3

% Detect transitions in each curve of y3
for col = 1:size(y3, 2)
    diff_y3 = diff(y3(:, col));
    trans = find(diff_y3 ~= 0);
    
    for t = 1:length(trans)
        idx = trans(t);
        if y3(idx + 1, col) == 1
            % Transition 0 -> 1: mark at idx + 1
            marked_idx = idx + 1;
        else
            % Transition 1 -> 0: mark at idx
            marked_idx = idx;
        end
        
        % Assign transition index based on curve
        if col == 1
            transition_indices_top = [transition_indices_top; marked_idx];
        elseif col == 2
            transition_indices_bottom = [transition_indices_bottom; marked_idx];
        end
    end
end

% Remove duplicate indices and sort
transition_indices_top = unique(transition_indices_top);
transition_indices_bottom = unique(transition_indices_bottom);

% Plot
figureHandle = figure;
for i = 1:5
    subplot(5, 1, i);
    
    switch i
        case 1
            plot(x, y1, 'LineWidth', 1.5);
            hold on;
        case 2
            plot(x, y2, 'LineWidth', 1.5);
            hold on;
        case 3
            plot(x, y3, 'LineWidth', 1.5);
            hold on;
        case 4
            plot(x, y4, 'LineWidth', 1.5);
            hold on;
        case 5
            plot(x, y5, 'LineWidth', 1.5);
            hold on;
    end
    
    % Ajouter des lignes verticales basées sur les indices de transition
    if i <= 3  % Subplots du haut
        for t = transition_indices_top'
            if t >= 1 && t <= length(x)
                x_val = x(t);
                plot([x_val, x_val], ylim, '--', 'LineWidth', 1.2, 'Color', [0.5, 0.5, 0.5]);
            end
        end
    end
    
    if i >= 3  % Subplots du bas
        for t = transition_indices_bottom'
            if t >= 1 && t <= length(x)
                x_val = x(t);
                plot([x_val, x_val], ylim, '--', 'LineWidth', 1.2, 'Color', [0.5, 0.5, 0.5]);
            end
        end
    end
    
    legend(legends{i}, 'Location', 'best', 'FontSize', text_size-3);
    xlabel('Length in MNI Space [mm]', 'FontSize', text_size);
    ylabel(ylabels{i}, 'FontSize', text_size);
    grid on;
    hold off;
    
    % Ajouter les annotations A/P si l'axe correspond à y
    if i == 3 && isnan(str2double(coord{2}))
%         xline(0, '--k', 'Ligne Médiane'); % Indique le plan médian
        xline(0, '--k');
        text(min(x), max(ylim), 'Anterior', 'HorizontalAlignment', 'left', 'FontSize', text_size);
        text(max(x), max(ylim), 'Posterior', 'HorizontalAlignment', 'right', 'FontSize', text_size);
    end
end

mainTitle = sprintf('Brain 1D Profile: Subject %d, qMRI Parameter %s, Profil (%s,%s,%s)', ...
    subject, qMRI, coord{1}, coord{2}, coord{3});
sgtitle(mainTitle, 'FontSize', text_size);

%% Save plot after managing it manually (full screen)
dir_out = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\qMRI_results\Brain1D';
pth_out = fullfile(dir_out, sprintf('%s_%s.%s.%s_S%d_merged_HR.eps',qMRI_param,coord{1},coord{2},coord{3},subject));
exportgraphics(figureHandle, pth_out, 'Resolution', 450); % file format defined in pth_out

%% Plot unsmoothed + TWS/TSPOON signals & GM/WM population masks
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
TCs = {'GM', 'WM', 'CSF'};

% Choose the one to plot
qMRI_params = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
qMRI_param_num = 4;
qMRI_param = qMRI_params{qMRI_param_num};
subject = 1; % (1:138)
coord = {'100', '', '100'}; % [181 217 181]
text_size = 12;

% Get the subject & qMRI parameter -specific MNI path
MNI_sub_dir = fullfile(ds_dir, 'derivatives', 'SPM12_dartel', sprintf('sub-%03d', subject), 'anat');
MNI_path = spm_select('ExtFPList', MNI_sub_dir, sprintf('^.*space-MNI_%s\\.nii$', qMRI_param));

% Get the subject & qMRI parameter -specific TWS path
TWS_sub_dir = fullfile(ds_dir, 'derivatives', 'AJ-TWS', sprintf('sub-%03d', subject), 'anat');
TWS_path = spm_select('ExtFPList', TWS_sub_dir, sprintf('^merged_.*%s\\.nii$', qMRI_param));

if isempty(TWS_path)
    dir = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives';
    x_paths = cell(1, 2);
    x_paths{1} = spm_select('ExtFPList', TWS_sub_dir, sprintf('^.*_GM.*%s\\.nii$', qMRI_param));
    x_paths{2} = spm_select('ExtFPList', TWS_sub_dir, sprintf('^.*_WM.*%s\\.nii$', qMRI_param));
    exMask_paths = cell(1, 2);
    exMask_paths{1} = spm_select('ExtFPList', dir, '^mask_WTA_GM_mask\\.nii$');
    exMask_paths{2} = spm_select('ExtFPList', dir, '^mask_WTA_WM_mask\\.nii$');
    TWS_path = aj_TS2One(x_paths, exMask_paths, qMRI_params);
end

% Get the subject & qMRI parameter -specific TSPOON path
TSPOON_sub_dir = fullfile(ds_dir, 'derivatives', 'AJ-TSPOON', sprintf('sub-%03d', subject), 'anat');
TSPOON_path = spm_select('ExtFPList', TSPOON_sub_dir, sprintf('^merged_.*%s\\.nii$', qMRI_param));

if isempty(TSPOON_path)
    dir = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives';
    x_paths = cell(1, 2);
    x_paths{1} = spm_select('ExtFPList', TSPOON_sub_dir, sprintf('^.*_GM.*%s\\.nii$', qMRI_param));
    x_paths{2} = spm_select('ExtFPList', TSPOON_sub_dir, sprintf('^.*_WM.*%s\\.nii$', qMRI_param));
    exMask_paths = cell(1, 2);
    exMask_paths{1} = spm_select('ExtFPList', dir, '^atlas-GM_space-MNI_mask\\.nii$');
    exMask_paths{2} = spm_select('ExtFPList', dir, '^atlas-WM_space-MNI_mask\\.nii$');
    TSPOON_path = aj_TS2One(x_paths, exMask_paths, qMRI_params);
end

% Get the population GM/WM masks (winner takes all)
deriv_dir = fullfile(ds_dir, 'derivatives');
WTAM_paths = cell(1, 2);
WTAM_paths{1} = spm_select('ExtFPList', deriv_dir, '^atlas-GM_space-MNI_mask\.nii$');
WTAM_paths{2} = spm_select('ExtFPList', deriv_dir, '^atlas-WM_space-MNI_mask\.nii$');

% Get 1D profiles for the first subplot
subplot1_paths = cell(1, 3);
subplot1_paths{1} = MNI_path;
subplot1_paths{2} = TWS_path;
subplot1_paths{3} = TSPOON_path;

[profiles_1, mni_coords] = aj_get_profile(subplot1_paths, coord);

% Get 1D profiles for the second subplot
subplot2_paths = cell(1, 2);
subplot2_paths{1} = WTAM_paths{1};
subplot2_paths{2} = WTAM_paths{2};

[profiles_2, ~] = aj_get_profile(subplot2_paths, coord);

% Group profiles
x = mni_coords{1}; % Use MNI coordinates
y1 = [profiles_1{1}(:), profiles_1{2}(:), profiles_1{3}(:)];
y2 = [profiles_2{1}(:), profiles_2{2}(:)];

% Plot y1 and y2
figureHandle = figure;

% Upper plot (y1)
subplot(2, 1, 1);
plot(x, y1, 'LineWidth', 1.5);
hold on;

xlabel('Length in MNI Space [mm]', 'FontSize', text_size);
if qMRI_param_num == 1
    qMRI = 'MTsat';
    ylabel('Intensity [AU]', 'FontSize', text_size);
elseif qMRI_param_num == 2
    qMRI = 'PD';
    ylabel('Intensity [%]', 'FontSize', text_size);
elseif qMRI_param_num == 3
    qMRI = 'R1';
    ylabel('Intensity [s^{-1}]', 'FontSize', text_size);
elseif qMRI_param_num == 4
    qMRI = 'R2*';
    ylabel('Intensity [s^{-1}]', 'FontSize', text_size);
end
grid on;
legend({'Not Smoothed', 'TWS', 'TSPOON'}, 'Location', 'best', 'FontSize', text_size);
hold off;

% Lower plot (y2)
subplot(2, 1, 2);
plot(x, y2, 'LineWidth', 1.5);
hold on;

% Annotate anterior/posterior for y1 if varying along Y
if isnan(str2double(coord{2}))
%     xline(0, '--k', 'Ligne Médiane');
    xline(0, '--k');
    text(min(x), max(ylim), 'Anterior', 'HorizontalAlignment', 'left', 'FontSize', text_size);
    text(max(x), max(ylim), 'Posterior', 'HorizontalAlignment', 'right', 'FontSize', text_size);
end

xlabel('Length in MNI [mm]', 'FontSize', text_size);
ylabel('Binary Mask Value', 'FontSize', text_size);
grid on;
legend({'GM', 'WM'}, 'Location', 'best', 'FontSize', text_size);
hold off;

mainTitle = sprintf('Brain 1D Merged Profile: Subject %d, qMRI Parameter %s, Profile (%s,%s,%s)',...
    subject, qMRI, coord{1}, coord{2}, coord{3});
sgtitle(mainTitle, 'FontSize', text_size);
