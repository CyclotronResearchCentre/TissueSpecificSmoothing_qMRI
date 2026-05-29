% aj_compare_smoothing_via_voxelwise_log_likelihood.m

% ---------------------------------------------------------
%% ============== USER SETTINGS ==============
clear all; close all; clc;

addpath('C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12');

spm('Defaults','fMRI'); spm_jobman('initcfg');

combination_names = {
    'MTsat_WM'
    'MTsat_WM'
    'PDmap_GM'
    'PDmap_WM'
    'R1map_GM'
    'R1map_WM'
    'R2starmap_GM'
    'R2starmap_WM'
};

basepaths = {
    'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-SUSAN_GLM_rft1\'
    'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TSPOON_GLM_rft1\'
    'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS_GLM_rft1\'
};

% basepaths = {
%     'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\SnPM-SUSAN_mreg-age\'
%     'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\SnPM-TSPOON_mreg-age\'
%     'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\SnPM-TWS_mreg-age\'
% };

out_path = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\LL_AIC_BIC_results\parametric_rft1';

SPMmat_names = 'SPM.mat';
beta_names = {'beta_0001.nii','beta_0002.nii','beta_0003.nii','beta_0004.nii','beta_0005.nii'};
ResMS_name = 'ResMS.nii';

% FOR SnPM
% SPMmat_names = 'SnPM.mat';
% beta_names = {'beta_0001.img','beta_0002.img','beta_0003.img','beta_0004.img','beta_0005.img'};
% ResMS_name = 'ResMS.img';

% Chunking to save memory
use_chunking = true;
chunk_vox = 60e3;  % voxels per block, depends on computer RAM

% Numerical floor for sigma2
sigma2_floor = 1e-8;

% Verbose
verbose = true;

for i = 1:numel(combination_names)
    % ---------- Determine tissue mask based on combination_names{i}  ----------
    if contains(combination_names{i}, 'WM')
        mask_fn = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\mask_WTA_WM_mask.nii';
    elseif contains(combination_names{i}, 'GM')
        mask_fn = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\mask_WTA_GM_mask.nii';
    else
        error('Neither WM nor GM found in combination name: %s', combination_names{i});
    end
    
    % ---------- Determine the output folder based on combination_names{i}  ----------
    outdir = fullfile(out_path, combination_names{i});
    if ~exist(outdir,'dir'), mkdir(outdir); end
    
    % ---------- Load SPM.mat files and basic checks ----------
    M = numel(basepaths);
    SPM_all = cell(M,1);
    scan_lists = cell(M,1);   % SPM.xY.P: list of subject images used in group GLM
    X_all = cell(M,1);        % design matrices
    SPMmat_paths = cell(M,1);
    beta_files_all = cell(M,1);
    ResMS_files = cell(M,1);

    for m = 1:M
        SPMmat_paths{m} = fullfile(basepaths{m}, combination_names{i}, SPMmat_names);
        beta_files_all{m} = fullfile(basepaths{m}, combination_names{i}, beta_names);
        ResMS_files{m}    = fullfile(basepaths{m}, combination_names{i}, ResMS_name);
    end

    for m = 1:M
        if ~exist(SPMmat_paths{m}, 'file'), error('SPM.mat not found: %s', SPMmat_paths{m}); end
        tmp = load(SPMmat_paths{m});
        if ~isfield(tmp,'SPM'), error('No SPM in %s', SPMmat_paths{m}); end
        SPM_all{m} = tmp.SPM;

        % get subject image list
        if isfield(tmp.SPM,'xY') && isfield(tmp.SPM.xY,'P')
            scan_lists{m} = tmp.SPM.xY.P(:);
        else
            error('SPM.xY.P not present in %s', SPMmat_paths{m});
        end

        % design matrix
        if isfield(tmp.SPM,'xX') && isfield(tmp.SPM.xX,'X')
            X_all{m} = tmp.SPM.xX.X;
        else
            error('SPM.xX.X not present in %s', SPMmat_paths{m});
        end

        if verbose, fprintf('Method %d: %d betas, %d scans\n', m, numel(beta_files_all{m}), numel(scan_lists{m})); end
    end

    % FOR SnPM
%     for m = 1:M
%         if ~exist(SPMmat_paths{m}, 'file')
%             error('SnPM.mat not found: %s', SPMmat_paths{m}); 
%         end
%         tmp = load(SPMmat_paths{m});
% 
%         % SnPM structure is directly in tmp (not tmp.SnPM)
%         if ~isfield(tmp,'V') || isempty(tmp.V)
%             error('No V field (image list) in %s', SPMmat_paths{m});
%         end
%         SnPM_all{m} = tmp;
% 
%         % get subject image list from V.fname (SnPM V structure contains image paths)
%         scan_lists{m} = {tmp.V.fname}';
%         scan_lists{m} = scan_lists{m}(~cellfun('isempty', scan_lists{m})); % clean empty
% 
%         % design matrix: reconstruct from CfgFile if needed, or use stored design
%         % SnPM doesn't store X directly, but design info is in CfgFile
%         if isfield(tmp,'CfgFile') && ischar(tmp.CfgFile)
%             cfg_tmp = load(tmp.CfgFile);
%             if isfield(cfg_tmp,'P')  % original scan list from config
%                 orig_scans{m} = cfg_tmp.P(:);
%             end
%             if isfield(cfg_tmp,'GX') % full design matrix often in GX
%                 X_all{m} = cfg_tmp.GX;
%             elseif isfield(cfg_tmp,'HCBGnames') % covariates info
%                 % reconstruct minimal design info if needed
%                 X_all{m} = [ones(138,1), 1:138]'; % intercept + subject index as fallback
%                 warning('GX not found, using minimal design in %s', SPMmat_paths{m});
%             end
%         else
%             % fallback: create dummy design matching scan count
%             X_all{m} = [ones(numel(scan_lists{m}),1)];
%             warning('No CfgFile, using intercept-only design in %s', SPMmat_paths{m});
%         end
% 
%         % 3. BETAS : chercher beta_*.img dans le même répertoire que SnPM.mat
%         snpm_dir = fileparts(SPMmat_paths{m});
%         beta_files = dir(fullfile(snpm_dir, 'beta_*.img'));
%         beta_files_all{m} = fullfile(snpm_dir, {beta_files.name}');
% 
%         % Vérifier qu'on a bien 4 betas (1 intérêt + 3 non-intérêt)
% %         expected_nbetas = size(cfg.GX,2);
% %         if numel(beta_files_all{m}) ~= expected_nbetas
% %             warning('Method %d: expected %d betas, found %d', ...
% %                 m, expected_nbetas, numel(beta_files_all{m}));
% %         end
% 
%         if verbose
%             fprintf('Method %d: %d images, %d betas, DF=%d\n', ...
%                 m, numel(scan_lists{m}), numel(beta_files_all{m}), tmp.df);
%         end
%     end

    old_root = 'C:\Users\lucad\Documents\smoothing\data\qMRI_AgingCallaghan\derivatives\';
    new_root = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\';

    for m = 1:M
        scan_lists{m} = cellfun(@(p) strrep(p, old_root, new_root), ...
                                scan_lists{m}, 'UniformOutput', false);
    end

    % Check number of subjects consistent across methods
    Ns = cellfun(@numel, scan_lists);
    if any(Ns ~= Ns(1)), error('Different number of subject images between methods: %s', mat2str(Ns)); end
    Nsub = Ns(1);
    if verbose, fprintf('Detected Nsub = %d subjects\n', Nsub); end

    % Check beta counts vs design columns
    p_all = cellfun(@numel, beta_files_all);
    if any(p_all ~= size(X_all{1},2))
        warning('Number of betas differs from design columns. Will use min(p,cols) per method.');
    end

    % ---------- Build mask / voxel list ----------
    if ~isempty(mask_fn) && exist(mask_fn,'file')
        Vmask = spm_vol(mask_fn);
        mask_vol = spm_read_vols(Vmask) > 0;
        mask_idx = find(mask_vol);
        if verbose, fprintf('Using provided mask: %d voxels\n', numel(mask_idx)); end
    else
        % derive intersection mask from first method's first N images (nonzero)
        Vref = spm_vol(scan_lists{1}{1});
        mask_vol = true(Vref.dim);
        for s = 1:Nsub
            img = spm_read_vols(spm_vol(scan_lists{1}{s}));
            mask_vol = mask_vol & (~isnan(img) & img~=0);
        end
        mask_idx = find(mask_vol);
        if verbose, fprintf('Derived mask from images: %d voxels\n', numel(mask_idx)); end
    end

    nvox_total = numel(mask_idx);

    % Chunking
    if use_chunking
        nblocks = ceil(nvox_total / chunk_vox);
        blocks = cell(nblocks,1);
        for b = 1:nblocks
            starti = (b-1)*chunk_vox + 1;
            endi = min(b*chunk_vox, nvox_total);
            blocks{b} = mask_idx(starti:endi);
        end
    else
        blocks = {mask_idx};
        nblocks = 1;
        disp('WARNINGGGG');
    end
    if verbose, fprintf('Processing %d voxels in %d blocks\n', nvox_total, nblocks); end

    % ---------- Prepare output volumes (templates) ----------
    Vtempl = spm_vol(scan_lists{1}{1}); % template header
    % prepare filenames
    LL_maps = cell(M,1); 
    
    BIC_maps = cell(M,1);
    for m = 1:M
        LL_maps{m} = fullfile(outdir, sprintf('LL_method%d.nii', m));
        % create empty volumes (zeros)
        Vo = Vtempl; Vo.fname = LL_maps{m}; Vo.dt = [16 0]; spm_write_vol(Vo, zeros(Vtempl.dim));
    end
    % best method maps
    best_LL_fn = fullfile(outdir,'best_method_by_LL.nii');
    Vo = Vtempl; Vo.fname = best_LL_fn; spm_write_vol(Vo, zeros(Vtempl.dim));

    % ---------- MAIN loop: per-block computations ----------
    % We'll store per-method per-voxel scalar LL_sum, AIC, BIC in arrays per block,
    % then write into the output volumes at the end of processing each block.

    for b = 1:nblocks
        voxlist = blocks{b};
        nb = numel(voxlist);
        if verbose, fprintf('Block %d/%d: voxels %d\n', b, nblocks, nb); end

        % Preallocate per-method results for this block
        LL_block = zeros(nb, M);

        % For each method:
        for m = 1:M
            if verbose, fprintf(' Method %d: loading beta images and ResMS\n', m); end

            % Load betas at voxlist: beta images are beta_files_all{m}{k}
            beta_files = beta_files_all{m};
            p = numel(beta_files);
            % limit p to number of cols in X if mismatch
            p_used = min(p, size(X_all{m},2));
            beta_vals = zeros(nb, p_used);
            for k = 1:p_used
                Vk = spm_vol(beta_files{k});
                imgk = spm_read_vols(Vk);
                beta_vals(:,k) = imgk(voxlist);
            end

            % Read ResMS image for method m
            Vres = spm_vol(ResMS_files{m});
            resms_vol = spm_read_vols(Vres);
            sigma2_block = resms_vol(voxlist); % voxel-wise residual variance (mean squares)
            sigma2_block = sigma2_block * (df / Nsub);
            % enforce numeric floor
            sigma2_block(isnan(sigma2_block)) = sigma2_floor;
            sigma2_block(sigma2_block < sigma2_floor) = sigma2_floor;

            % Read subject images Y_s at these voxels (from SPM.xY.P)
            Y_block = zeros(nb, Nsub);
            for s = 1:Nsub
                Vs = spm_vol(scan_lists{m}{s});
                imgs = spm_read_vols(Vs);
                Y_block(:,s) = imgs(voxlist);
            end

            % Reconstruct predicted values mu_s(v) = X(s,:) * beta(v,:)'
            % For each voxel, mu(:,v) = X * beta_voxel'
            % We'll compute matrix multiplication: mu = X * B' where B is nb x p_used
            % So mu' = B * X' ? We need mu as nb x Nsub
            % B (nb x p_used), X_all{m} is Nsub x p_all
            Xm = X_all{m}(:, 1:p_used); % Nsub x p_used
            mu_block = (Xm * beta_vals')'; % (p_used dims) => result nb x Nsub

            % Now compute log-likelihood per voxel:
            % ll(v) = sum_s -0.5*( log(2*pi*sigma2(v)) + (Y_s(v)-mu_s(v))^2 / sigma2(v) )
            % compute squared residuals per voxel across subjects:
            resid2 = (Y_block - mu_block).^2; % nb x Nsub
            % divide by sigma2 (nb x 1) broadcast:
            term = bsxfun(@rdivide, sum(resid2, 2), sigma2_block); % nb x 1 (sum over s)
            % compute log-term: sum_s log(2*pi*sigma2) = Nsub * log(2*pi*sigma2)
            ll_block = -0.5 * ( Nsub .* log(2*pi*sigma2_block) + term ); % nb x 1

            % Store LL
            LL_block(:, m) = ll_block;

            % AIC and BIC: k = number of parameters in model (p_used), n = number of observations (Nsub)
%             k = p_used;
%             n = Nsub;
        end % method loop

        % Decide best method per voxel according to metrics
        % For LL: choose argmax(LL)
        [~, best_by_LL_idx] = max(LL_block, [], 2); % 1..M per voxel

        % === Write block results back into NIfTI volumes ===
        for m = 1:M
            Vout = spm_vol(LL_maps{m});
            cur = spm_read_vols(Vout);
            cur(voxlist) = LL_block(:, m);
            spm_write_vol(Vout, cur);
        end

        % write best method maps (store index of best method)
        Vout = spm_vol(best_LL_fn); cur = spm_read_vols(Vout); cur(voxlist) = best_by_LL_idx; spm_write_vol(Vout, cur);

        if verbose, fprintf(' Block %d written (LL/AIC/BIC & best maps updated)\n', b); end
    end % block

    % Save summary table: global sums of LL per method (sum across voxels in mask)
    LL_total = zeros(M,1);
    for m = 1:M
        V = spm_vol(LL_maps{m});
        img = spm_read_vols(V);
        LL_total(m) = sum(img(mask_idx));
    end
    
    fprintf('Done. Outputs saved in %s\n', outdir);
    fprintf('LL_total per method: %s\n', mat2str(LL_total'));
end
