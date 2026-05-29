%% aj_get_voxelwise_log_likelihood
% =========================================================================
% This script compares several spatial smoothing approaches using voxelwise
% Gaussian log-likelihood (LL) estimation derived from GLM residuals.
%
% For each qMRI contrast and tissue class combination (e.g. MTsat_GM),
% the script:
%
%   1. Loads the subject-level GLM information for each smoothing method
%      (design matrix, beta images, residual variance maps, subject scans).
%
%   2. Reconstructs voxelwise predicted values from:
%           Y_hat = X * beta
%
%   3. Computes the voxelwise Gaussian log-likelihood:
%
%      LL(v) = -0.5 * [ N * log(2*pi*sigma²(v))
%                       + sum(residual²(v))/sigma²(v) ]
%
%   4. Determines, for each voxel, which smoothing strategy maximizes
%      the likelihood.
%
%   5. Saves:
%        - voxelwise LL maps for each smoothing method
%        - a voxelwise "best method" map
%        - summary LL values across the mask
%
% -------------------------------------------------------------------------
% INPUT DATA
% -------------------------------------------------------------------------
% The script supports:
%
%   - SPM parametric GLM analyses
%   - SnPM non-parametric analyses
%
% depending on the value of "analysis_type".
%
% -------------------------------------------------------------------------
% REQUIRED FILES
% -------------------------------------------------------------------------
% For each smoothing method and each combination:
%
%   SPM:
%       - SPM.mat
%       - beta_XXXX.nii
%       - ResMS.nii
%
%   SnPM:
%       - SnPM.mat
%       - beta_XXXX.img
%       - ResMS.img
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% For each combination:
%
%   LL_method1.nii
%   LL_method2.nii
%   LL_method3.nii
%   best_method_by_LL.nii
%
% where:
%   1 = SUSAN
%   2 = TSPOON
%   3 = TWS
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Antoine Jacquemin
% GIGA-CRC In Vivo Imaging
% University of Liège, Belgium
% =========================================================================

function aj_get_voxelwise_log_likelihood(cfg)
%% ========================================================================
% INITIALIZATION
% =========================================================================

addpath(cfg.spm_path);

spm('Defaults','fMRI');
spm_jobman('initcfg');

analysis_type = cfg.analysis_type;
combination_names = cfg.combination_names;
basepaths        = cfg.basepaths;
out_path         = fullfile(cfg.out_root, 'results_vwLL', cfg.analysis_type);

SPMmat_name = cfg.SPMmat_name;
beta_names  = cfg.beta_names;
ResMS_name  = cfg.ResMS_name;

method_names = cfg.method_names;

use_chunking = cfg.chunking.use;
chunk_vox = cfg.chunking.chunk_vox;
sigma2_floor = cfg.sigma2_floor;

verbose = cfg.verbose;

%% ========================================================================
% MAIN LOOP OVER qMRI COMBINATIONS
% =========================================================================
for i = 1:numel(combination_names)

    combination = combination_names{i};

    fprintf('\n====================================================\n');
    fprintf('Processing: %s\n', combination);
    fprintf('====================================================\n');

    %% --------------------------------------------------------------------
    % Select tissue mask
    % ---------------------------------------------------------------------
    if contains(combination, 'WM')

        mask_fn = ...
            'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\mask_WTA_WM_mask.nii';

    elseif contains(combination, 'GM')

        mask_fn = ...
            'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\mask_WTA_GM_mask.nii';

    else
        error('Neither GM nor WM found in "%s".', combination);
    end

    %% --------------------------------------------------------------------
    % Create output directory
    % ---------------------------------------------------------------------
    outdir = fullfile(out_path, combination);

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    %% --------------------------------------------------------------------
    % Initialize containers
    % ---------------------------------------------------------------------
    M = numel(basepaths);

    scan_lists     = cell(M,1);
    X_all          = cell(M,1);
    beta_files_all = cell(M,1);
    ResMS_files    = cell(M,1);
    SPMmat_paths   = cell(M,1);

    %% --------------------------------------------------------------------
    % Define paths
    % ---------------------------------------------------------------------
    for m = 1:M

        SPMmat_paths{m} = ...
            fullfile(basepaths{m}, combination, SPMmat_name);

        beta_files_all{m} = ...
            fullfile(basepaths{m}, combination, beta_names);

        ResMS_files{m} = ...
            fullfile(basepaths{m}, combination, ResMS_name);
    end

    %% --------------------------------------------------------------------
    % Load GLM information
    % ---------------------------------------------------------------------
    for m = 1:M

        if ~exist(SPMmat_paths{m}, 'file')
            error('File not found:\n%s', SPMmat_paths{m});
        end

        tmp = load(SPMmat_paths{m});

        switch lower(analysis_type)

            % ==============================================================
            % PARAMETRIC SPM
            % ==============================================================
            case {'spm_rft1_age_linear', 'spm_rft1_age_nonlinear'}

                if ~isfield(tmp, 'SPM')
                    error('No SPM structure found.');
                end

                SPM = tmp.SPM;

                scan_lists{m} = SPM.xY.P(:);
                X_all{m}      = SPM.xX.X;

            % ==============================================================
            % SNPM
            % ==============================================================
            case 'snpm'

                if ~isfield(tmp, 'V')
                    error('No image list found in SnPM structure.');
                end

                scan_lists{m} = {tmp.V.fname}';
                scan_lists{m} = ...
                    scan_lists{m}(~cellfun('isempty', scan_lists{m}));

                if isfield(tmp, 'CfgFile')

                    cfg_tmp = load(tmp.CfgFile);

                    if isfield(cfg_tmp, 'GX')

                        X_all{m} = cfg_tmp.GX;

                    else

                        X_all{m} = ...
                            ones(numel(scan_lists{m}),1);

                        warning('GX not found. Using intercept-only model.');
                    end

                else

                    X_all{m} = ones(numel(scan_lists{m}),1);

                    warning('CfgFile missing. Using intercept-only model.');
                end
        end

        if verbose
            fprintf('Method %d: %d scans | %d betas\n', ...
                m, ...
                numel(scan_lists{m}), ...
                numel(beta_files_all{m}));
        end
    end

    %% ------------------------------------------------------------------------
    % Robust path remapping: keep everything after "derivatives"
    % -------------------------------------------------------------------------

    new_root = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\';

    for m = 1:M

        scan_lists{m} = cellfun(@(p) local_rebase_path(p, new_root), ...
                                 scan_lists{m}, ...
                                 'UniformOutput', false);
    end

    %% --------------------------------------------------------------------
    % Check consistency across methods
    % ---------------------------------------------------------------------
    Ns = cellfun(@numel, scan_lists);

    if any(Ns ~= Ns(1))
        error('Different number of subject images between methods.');
    end

    Nsub = Ns(1);

    fprintf('Detected %d subjects.\n', Nsub);

    %% --------------------------------------------------------------------
    % Create voxel mask
    % ---------------------------------------------------------------------
    Vmask   = spm_vol(mask_fn);
    maskVol = spm_read_vols(Vmask) > 0;
    maskIdx = find(maskVol);

    fprintf('Mask contains %d voxels.\n', numel(maskIdx));

    %% --------------------------------------------------------------------
    % Chunking
    % ---------------------------------------------------------------------
    nvox_total = numel(maskIdx);

    if use_chunking

        nblocks = ceil(nvox_total / chunk_vox);
        blocks  = cell(nblocks,1);

        for b = 1:nblocks

            starti = (b-1)*chunk_vox + 1;
            endi   = min(b*chunk_vox, nvox_total);

            blocks{b} = maskIdx(starti:endi);
        end

    else

        blocks  = {maskIdx};
        nblocks = 1;
    end

    fprintf('Processing %d blocks.\n', nblocks);

    %% ------------------------------------------------------------------------
    % Prepare output volumes
    % -------------------------------------------------------------------------
    Vtempl = spm_vol(scan_lists{1}{1});

    LL_maps = cell(M,1);

    % -------------------------------------------------------------------------
    % Smoothing method names
    % -------------------------------------------------------------------------
    method_names = {
        'SUSAN'
        'TSPOON'
        'TWS'
    };

    for m = 1:M

        % ---------------------------------------------------------------------
        % Output filename:
        %   LL_method1_SUSAN.nii
        %   LL_method2_TSPOON.nii
        %   LL_method3_TWS.nii
        % ---------------------------------------------------------------------
        LL_maps{m} = fullfile( ...
            outdir, ...
            sprintf('LL_method%d_%s.nii', ...
            m, method_names{m}));

        % Create empty output volume
        Vo = Vtempl;
        Vo.fname = LL_maps{m};
        Vo.dt = [16 0];

        spm_write_vol(Vo, zeros(Vtempl.dim));
    end

    %% ------------------------------------------------------------------------
    % Best-method map
    %
    % Voxel values:
    %   1 = SUSAN
    %   2 = TSPOON
    %   3 = TWS
    % -------------------------------------------------------------------------
    best_LL_fn = fullfile(outdir, 'best_method_by_LL.nii');

    Vo = Vtempl;
    Vo.fname = best_LL_fn;

    spm_write_vol(Vo, zeros(Vtempl.dim));

    %% ====================================================================
    % MAIN VOXELWISE COMPUTATION
    % =====================================================================
    for b = 1:nblocks

        voxlist = blocks{b};
        nb      = numel(voxlist);

        fprintf('Block %d / %d\n', b, nblocks);

        LL_block = zeros(nb, M);

        %% ----------------------------------------------------------------
        % Loop across smoothing methods
        % -----------------------------------------------------------------
        for m = 1:M

            fprintf('Computing method %d...\n', m);

            % --------------------------------------------------------------
            % Load beta maps
            % --------------------------------------------------------------
            beta_files = beta_files_all{m};

            p      = numel(beta_files);
            p_used = min(p, size(X_all{m},2));

            beta_vals = zeros(nb, p_used);

            for k = 1:p_used

                Vk   = spm_vol(beta_files{k});
                imgk = spm_read_vols(Vk);

                beta_vals(:,k) = imgk(voxlist);
            end

            % --------------------------------------------------------------
            % Residual variance
            % --------------------------------------------------------------
            Vres      = spm_vol(ResMS_files{m});
            resms_vol = spm_read_vols(Vres);

            sigma2_block = resms_vol(voxlist);

            sigma2_block(isnan(sigma2_block)) = sigma2_floor;
            sigma2_block(sigma2_block < sigma2_floor) = sigma2_floor;

            % --------------------------------------------------------------
            % Load subject images
            % --------------------------------------------------------------
            Y_block = zeros(nb, Nsub);

            for s = 1:Nsub

                Vs   = spm_vol(scan_lists{m}{s});
                imgs = spm_read_vols(Vs);

                Y_block(:,s) = imgs(voxlist);
            end

            % --------------------------------------------------------------
            % Predicted values
            % --------------------------------------------------------------
            Xm = X_all{m}(:, 1:p_used);

            mu_block = (Xm * beta_vals')';

            % --------------------------------------------------------------
            % Log-likelihood computation
            % --------------------------------------------------------------
            resid2 = (Y_block - mu_block).^2;

            term = bsxfun( ...
                @rdivide, ...
                sum(resid2,2), ...
                sigma2_block);

            ll_block = -0.5 * ( ...
                Nsub .* log(2*pi*sigma2_block) + term);

            LL_block(:,m) = ll_block;
        end

        %% ----------------------------------------------------------------
        % Determine best smoothing method
        % -----------------------------------------------------------------
        [~, best_by_LL_idx] = max(LL_block, [], 2);

        %% ----------------------------------------------------------------
        % Write LL maps
        % -----------------------------------------------------------------
        for m = 1:M

            Vout = spm_vol(LL_maps{m});

            cur = spm_read_vols(Vout);

            cur(voxlist) = LL_block(:,m);

            spm_write_vol(Vout, cur);
        end

        %% ----------------------------------------------------------------
        % Write best-method map
        % -----------------------------------------------------------------
        Vout = spm_vol(best_LL_fn);

        cur = spm_read_vols(Vout);

        cur(voxlist) = best_by_LL_idx;

        spm_write_vol(Vout, cur);

        fprintf('Block %d completed.\n', b);
    end

    %% ------------------------------------------------------------------------
    % Global LL summary
    % -------------------------------------------------------------------------
    LL_total = zeros(M,1);

    for m = 1:M

        V   = spm_vol(LL_maps{m});
        img = spm_read_vols(V);

        LL_total(m) = sum(img(maskIdx));
    end

    %% ------------------------------------------------------------------------
    % Display summary
    % -------------------------------------------------------------------------
    fprintf('\nDone.\n');
    fprintf('Outputs saved in:\n%s\n', outdir);

    fprintf('\nTotal log-likelihood:\n');

    for m = 1:M
        fprintf('  Method %d (%s): %.3f\n', ...
            m, method_names{m}, LL_total(m));
    end
end
end

%% ------------------------------------------------------------------------
% Helper function
% -------------------------------------------------------------------------
function new_path = local_rebase_path(old_path, new_root)

    % Find "derivatives" keyword in path
    token = 'derivatives';

    idx = strfind(lower(old_path), lower(token));

    if isempty(idx)
        error('Path does not contain "derivatives": %s', old_path);
    end

    % Keep everything AFTER "derivatives"
    cut_pos = idx(1) + length(token);

    relative_part = old_path(cut_pos+1:end);

    % Rebuild clean path
    new_path = fullfile(new_root, relative_part);

    % Fix mixed separators just in case
    new_path = strrep(new_path, '/', filesep);
    new_path = strrep(new_path, '\', filesep);
end