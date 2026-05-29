function cfg = aj_config_stat_tools()
% =========================================================================
% aj_config_stats
% =========================================================================
% Central configuration file for:
%
%   1. Cluster extraction analyses
%   2. Similarity metric computations
%
% This function centralizes:
%   - dataset locations
%   - smoothing-method directories
%   - qMRI combinations
%   - output folders
%   - Excel export settings
%
% The goal is to:
%   - avoid duplicated hardcoded paths
%   - standardize naming conventions
%   - simplify reproducibility
%   - facilitate future extensions
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% cfg : structure containing all analysis settings
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Antoine Jacquemin
% GIGA-CRC In Vivo Imaging
% University of Liège, Belgium
% =========================================================================

%% ========================================================================
% GENERAL PATHS
% =========================================================================

cfg.derivatives_root = ...
    'E:\Master_Thesis\Data\BIDS_AgingData\derivatives';

%% ========================================================================
% SPM TOOLBOX
% =========================================================================

cfg.spm_path = ...
    'C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12';

%% ========================================================================
% qMRI COMBINATIONS
% =========================================================================

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

%% ========================================================================
% ANALYSIS TYPE
% =========================================================================
% Available:
%   'snpm'
%   'rft0'
%   'rft1_age_linear'
%   'rft1_age_nonlinear'
% =========================================================================

cfg.analysis = 'rft0';

%% ========================================================================
% SMOOTHING METHODS
% =========================================================================
% Available:
%   'SUSAN'
%   'TSPOON'
%   'TWS'
% =========================================================================

cfg.method1 = 'SUSAN';
cfg.method2 = 'TWS';

%% ========================================================================
% FLAGS
% =========================================================================

cfg.flag = struct();

cfg.flag.saveExcel = true;

%% ========================================================================
% BASE DIRECTORIES
% =========================================================================

switch lower(cfg.analysis)

    % =====================================================================
    % SnPM
    % =====================================================================
    case 'snpm'

        cfg.base_dirs.SUSAN = fullfile( ...
            cfg.derivatives_root, ...
            'SnPM-SUSAN_mreg-age');

        cfg.base_dirs.TSPOON = fullfile( ...
            cfg.derivatives_root, ...
            'SnPM-TSPOON_mreg-age');

        cfg.base_dirs.TWS = fullfile( ...
            cfg.derivatives_root, ...
            'SnPM-TWS_mreg-age');

        cfg.output_dir = fullfile( ...
            cfg.derivatives_root, ...
            'results_simMetrics', ...
            'snpm');

    % =====================================================================
    % SPM RFT0
    % =====================================================================
    case 'rft0'

        cfg.base_dirs.SUSAN = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-SUSAN_GLM_rft0_age_linear');

        cfg.base_dirs.TSPOON = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-TSPOON_GLM_rft0_age_linear');

        cfg.base_dirs.TWS = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-TWS_GLM_rft0_age_linear');

        cfg.output_dir = fullfile( ...
            cfg.derivatives_root, ...
            'results_simMetrics', ...
            'rft0');

    % =====================================================================
    % SPM RFT1 LINEAR
    % =====================================================================
    case 'rft1_age_linear'

        cfg.base_dirs.SUSAN = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-SUSAN_GLM_rft1_age_linear');

        cfg.base_dirs.TSPOON = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-TSPOON_GLM_rft1_age_linear');

        cfg.base_dirs.TWS = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-TWS_GLM_rft1_age_linear');

        cfg.output_dir = fullfile( ...
            cfg.derivatives_root, ...
            'results_simMetrics', ...
            'rft1_age_linear');

    % =====================================================================
    % SPM RFT1 NONLINEAR
    % =====================================================================
    otherwise

        cfg.base_dirs.SUSAN = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-SUSAN_GLM_rft1_age_nonlinear');

        cfg.base_dirs.TSPOON = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-TSPOON_GLM_rft1_age_nonlinear');

        cfg.base_dirs.TWS = fullfile( ...
            cfg.derivatives_root, ...
            'AJ-TWS_GLM_rft1_age_nonlinear');

        cfg.output_dir = fullfile( ...
            cfg.derivatives_root, ...
            'results_simMetrics', ...
            'rft1_age_nonlinear');
end

%% ========================================================================
% AUTOMATIC METHOD PATHS
% =========================================================================

cfg.base_dir1 = cfg.base_dirs.(cfg.method1);
cfg.base_dir2 = cfg.base_dirs.(cfg.method2);

%% ========================================================================
% AUTOMATIC EXCEL OUTPUTS
% =========================================================================

if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

% -------------------------------------------------------------------------
% Similarity metrics Excel output
% -------------------------------------------------------------------------
cfg.excel.simMetrics = fullfile( ...
    cfg.output_dir, ...
    sprintf('SimMetrics_%s_vs_%s.xlsx', ...
    cfg.method1, ...
    cfg.method2));

% -------------------------------------------------------------------------
% Cluster data Excel output
% -------------------------------------------------------------------------
cfg.excel.clusterData = fullfile( ...
    cfg.base_dir1, ...
    'clusterdata.xlsx');

end