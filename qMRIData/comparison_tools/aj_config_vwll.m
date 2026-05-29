function cfg = aj_config_vwll()
% =========================================================================
% CONFIGURATION FILE FOR VOXELWISE LOG-LIKELIHOOD ANALYSIS
% =========================================================================
% This function centralizes all user-defined parameters for the
% voxelwise smoothing comparison pipeline.
%
% OUTPUT
%   cfg : structure containing all paths, options, and constants
% =========================================================================

%% ------------------------------------------------------------------------
% GENERAL SETTINGS
% -------------------------------------------------------------------------
cfg.analysis_type = 'snpm';

cfg.verbose = true;

%% ------------------------------------------------------------------------
% PATHS
% -------------------------------------------------------------------------
cfg.spm_path = ...
    'C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12';

cfg.out_root = ...
    'E:\Master_Thesis\Data\BIDS_AgingData\derivatives';

%% ------------------------------------------------------------------------
% COMBINATIONS (qMRI + tissue class)
% -------------------------------------------------------------------------
cfg.combination_names = {
    'MTsat_GM'
    'MTsat_WM'
    'PDmap_GM'
    'PDmap_WM'
    'R1map_GM'
    'R1map_WM'
    'R2starmap_GM'
    'R2starmap_WM'
};

%% ------------------------------------------------------------------------
% METHOD LABELS
% -------------------------------------------------------------------------
cfg.method_names = {
    'SUSAN'
    'TSPOON'
    'TWS'
};

%% ------------------------------------------------------------------------
% MASKS
% -------------------------------------------------------------------------
cfg.mask.GM = fullfile(cfg.out_root, 'mask_WTA_GM_mask.nii');
cfg.mask.WM = fullfile(cfg.out_root, 'mask_WTA_WM_mask.nii');

%% ------------------------------------------------------------------------
% ANALYSIS-SPECIFIC SETTINGS
% -------------------------------------------------------------------------
switch lower(cfg.analysis_type)

    case 'snpm'

        cfg.basepaths = {
            fullfile(cfg.out_root, 'SnPM-SUSAN_mreg-age')
            fullfile(cfg.out_root, 'SnPM-TSPOON_mreg-age')
            fullfile(cfg.out_root, 'SnPM-TWS_mreg-age')
        };

        cfg.SPMmat_name = 'SnPM.mat';

        cfg.beta_names = {
            'beta_0001.img'
            'beta_0002.img'
            'beta_0003.img'
            'beta_0004.img'
            'beta_0005.img'
        };

        cfg.ResMS_name = 'ResMS.img';

    case 'spm_rft1_age_linear'

        cfg.basepaths = {
            fullfile(cfg.out_root, 'AJ-SUSAN_GLM_rft1_age_linear')
            fullfile(cfg.out_root, 'AJ-TSPOON_GLM_rft1_age_linear')
            fullfile(cfg.out_root, 'AJ-TWS_GLM_rft1_age_linear')
        };

        cfg.SPMmat_name = 'SPM.mat';

        cfg.beta_names = {
            'beta_0001.nii'
            'beta_0002.nii'
            'beta_0003.nii'
            'beta_0004.nii'
            'beta_0005.nii'
        };

        cfg.ResMS_name = 'ResMS.nii';

    otherwise

        cfg.basepaths = {
            fullfile(cfg.out_root, 'AJ-SUSAN_GLM_rft1_age_nonlinear')
            fullfile(cfg.out_root, 'AJ-TSPOON_GLM_rft1_age_nonlinear')
            fullfile(cfg.out_root, 'AJ-TWS_GLM_rft1_age_nonlinear')
        };

        cfg.SPMmat_name = 'SPM.mat';

        cfg.beta_names = {
            'beta_0001.nii'
            'beta_0002.nii'
            'beta_0003.nii'
            'beta_0004.nii'
            'beta_0005.nii'
        };

        cfg.ResMS_name = 'ResMS.nii';
end

%% ------------------------------------------------------------------------
% RUNTIME OPTIONS
% -------------------------------------------------------------------------
cfg.chunking.use = true;
cfg.chunking.chunk_vox = 60e3;

cfg.sigma2_floor = 1e-8;

end