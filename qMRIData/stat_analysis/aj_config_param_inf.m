function cfg = aj_config_param_inf()
% =========================================================================
% aj_config_param
% =========================================================================
% Configuration file for parametric qMRI GLM analyses.
%
% This file centralizes:
%   - dataset paths
%   - smoothing settings
%   - SPM parameters
%   - qMRI contrasts
%   - tissue classes
%   - statistical inference options
%
% =========================================================================

%% ========================================================================
% PATHS
% =========================================================================

cfg.dataset_dir = ...
    'E:\Master_Thesis\Data\BIDS_AgingData';

cfg.spm_path = ...
    'C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12';

%% ========================================================================
% SMOOTHING METHODS
% =========================================================================

cfg.smoothing_methods = {
    'TWsmoot'
    'TWS'
    'TSPOON'
    'SUSAN'
};

% -------------------------------------------------------------------------
% Selected smoothing method
%
% 1 = TWsmoot
% 2 = TWS
% 3 = TSPOON
% 4 = SUSAN
% -------------------------------------------------------------------------

cfg.smoothing_index = 4;

cfg.smoothing_name = ...
    cfg.smoothing_methods{cfg.smoothing_index};

%% ========================================================================
% SMOOTHING DIRECTORY
% =========================================================================

switch lower(cfg.smoothing_name)

    case 'twsmoot'

        cfg.smoothing_dir = fullfile( ...
            cfg.dataset_dir, ...
            'derivatives', ...
            'TWsmoot');

    case 'susan'

        cfg.smoothing_dir = fullfile( ...
            cfg.dataset_dir, ...
            'derivatives', ...
            'SUSAN');

    otherwise

        cfg.smoothing_dir = fullfile( ...
            cfg.dataset_dir, ...
            'derivatives', ...
            sprintf('AJ-%s', cfg.smoothing_name));
end

%% ========================================================================
% SPM BATCH
% =========================================================================

cfg.jobfile = {fullfile(pwd, 'batch_spm_ageNonLinear_job.m')};

%% ========================================================================
% RANDOM FIELD THEORY PARAMETERS
% =========================================================================

cfg.rft_hypothesis = 2;

cfg.rft_names = {
    'rft0'
    'rft1'
};

cfg.rft_name = ...
    cfg.rft_names{cfg.rft_hypothesis};

% -------------------------------------------------------------------------
% SPM DEFAULTS
% -------------------------------------------------------------------------

cfg.rft_nonstationary = 1;

cfg.ufp = 0.5;

%% ========================================================================
% qMRI METRICS
% =========================================================================

cfg.qmetrics = {
    'MTsat'
    'PDmap'
    'R1map'
    'R2starmap'
};

%% ========================================================================
% TISSUE CLASSES
% =========================================================================

cfg.tissues = {
    'GM'
    'WM'
};

end