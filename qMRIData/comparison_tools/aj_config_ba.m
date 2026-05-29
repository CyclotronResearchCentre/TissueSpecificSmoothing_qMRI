function cfg = aj_config_ba()
% =========================================================================
% AJ_CONFIG_BA
% =========================================================================
% Centralized configuration file for the Bland–Altman analysis pipeline.
%
% This function defines:
%   - SPM installation path
%   - qMRI metric / tissue combinations
%   - smoothing methods to compare
%   - analysis options
%   - input/output directories
%   - anatomical masks
%
% The configuration structure is intended to standardize the Bland–Altman
% workflow and avoid hard-coded parameters within analysis scripts.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% cfg : structure containing all user-defined settings and paths.
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Antoine Jacquemin
% GIGA-CRC In Vivo Imaging
% University of Liège, Belgium
% =========================================================================

%% ========================================================================
% SPM PATH
% =========================================================================
% Path to the SPM12 installation directory.
% -------------------------------------------------------------------------
cfg.spm_path = ...
'C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12';

%% ========================================================================
% qMRI COMBINATIONS
% =========================================================================
% qMRI contrasts and tissue classes included in the analysis.
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

%% ========================================================================
% SMOOTHING METHODS
% =========================================================================
% Available methods:
%   - 'SUSAN'
%   - 'TSPOON'
%   - 'TWS'
%
% The Bland–Altman difference is computed as:
%       method1 - method2
% -------------------------------------------------------------------------
cfg.method1_name = 'TWS';
cfg.method2_name = 'TSPOON';

%% ========================================================================
% ANALYSIS FLAGS
% =========================================================================
% analysis:
%   - 'snpm'
%   - 'spm_rft0'
%   - 'spm_rft1'
%
% thresholded:
%   Apply threshold masks to statistical maps.
%
% useLL:
%   Use voxelwise log-likelihood (LL) maps instead of statistical maps.
% -------------------------------------------------------------------------
cfg.flag = struct();

cfg.flag.analysis    = 'snpm';
cfg.flag.drawPlot    = false;
cfg.flag.savePlot    = true;
cfg.flag.saveExcel   = true;
cfg.flag.paperReady  = true;

cfg.flag.thresholded = true;
cfg.flag.useLL       = false;

%% ========================================================================
% ROOT DIRECTORIES
% =========================================================================
% Define root directories and anatomical masks depending on the selected
% analysis type.
% -------------------------------------------------------------------------

% Base derivatives directory
cfg.derivatives_dir = ...
'E:\Master_Thesis\Data\BIDS_AgingData\derivatives';

%% ========================================================================
% INPUT DATA TYPE
% =========================================================================
% Two input modes are supported:
%
%   1. Statistical maps
%      → derived from SPM or SnPM analyses
%
%   2. Voxelwise log-likelihood maps
%      → LL_methodX.nii
%
% When LL maps are used:
%   - thresholding is automatically disabled
%   - analysis type is forced to 'spm_rft1'
% -------------------------------------------------------------------------
if cfg.flag.useLL
    cfg.flag.thresholded = false;

    % LL maps are generated from the voxelwise likelihood framework
    cfg.flag.analysis = 'spm_rft1';

    cfg.root_dir = fullfile( ...
        cfg.derivatives_dir, ...
        'results_LL');
    
    % Anatomical masks
    cfg.flag.maskGM = fullfile(cfg.derivatives_dir, 'mask_WTA_GM_mask.nii');
    cfg.flag.maskWM = fullfile(cfg.derivatives_dir, 'mask_WTA_WM_mask.nii');
else
    cfg.root_dir = cfg.derivatives_dir;
    
    % Anatomical masks
    cfg.flag.maskGM = '';
    cfg.flag.maskWM = '';
end

%% ========================================================================
% CONSISTENCY CHECKS
% =========================================================================

valid_methods = {'SUSAN','TSPOON','TWS'};

if ~ismember(cfg.method1_name, valid_methods)
    error('Invalid method1_name: %s', cfg.method1_name);
end

if ~ismember(cfg.method2_name, valid_methods)
    error('Invalid method2_name: %s', cfg.method2_name);
end

valid_analyses = {'snpm','spm_rft0','spm_rft1'};

if ~ismember(lower(cfg.flag.analysis), valid_analyses)
    error('Invalid analysis type: %s', cfg.flag.analysis);
end

end
