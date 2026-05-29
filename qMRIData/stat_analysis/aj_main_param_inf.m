%% aj_main_param_inf
% =========================================================================
% aj_main_param_inf
% =========================================================================
% Main execution script for parametric voxelwise GLM analyses performed on
% smoothed quantitative MRI (qMRI) maps derived from the AgingData dataset.
%
% The pipeline reproduces the statistical framework described in:
%
%   Callaghan et al. (2014)
%
% and applies several spatial smoothing strategies prior to group-level
% statistical inference.
%
% -------------------------------------------------------------------------
% PIPELINE OVERVIEW
% -------------------------------------------------------------------------
%
% 1. Initialize SPM environment
%
% 2. Generate covariate/regressor files:
%       - age
%       - total intracranial volume (TIV)
%       - scanner
%       - sex
%
% 3. Select smoothed qMRI maps for:
%       - MTsat
%       - PDmap
%       - R1map
%       - R2starmap
%
% 4. Run voxelwise GLM analyses separately for:
%       - gray matter (GM)
%       - white matter (WM)
%
% 5. Apply either:
%       - stationary RFT inference
%       - non-stationary RFT inference
%
% -------------------------------------------------------------------------
% SUPPORTED SMOOTHING METHODS
% -------------------------------------------------------------------------
%
%   1 = TWsmoot
%   2 = TWS
%   3 = TSPOON
%   4 = SUSAN
%
% -------------------------------------------------------------------------
% REFERENCES
% -------------------------------------------------------------------------
%
% Callaghan MF et al. (2014)
% Neurobiology of Aging
% http://dx.doi.org/10.1016/j.neurobiolaging.2014.02.008
%
% AgingData dataset:
% https://github.com/CyclotronResearchCentre/BIDS_AgingData
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Antoine Jacquemin
% GIGA-CRC In Vivo Imaging
% University of Liège, Belgium
% =========================================================================

%% ========================================================================
% INITIALIZATION
% =========================================================================
clear;
close all;
clc;

%% ========================================================================
% LOAD CONFIGURATION
% =========================================================================
cfg = aj_config_param_inf();

%% ========================================================================
% INITIALIZE SPM
% =========================================================================
addpath(cfg.spm_path);

spm_jobman('initcfg');
spm('defaults', 'fmri');

%% ========================================================================
% SPM DEFAULT PARAMETERS
% =========================================================================

% -------------------------------------------------------------------------
% Non-stationary random field correction
% -------------------------------------------------------------------------
spm_get_defaults('stats.rft.nonstat', cfg.rft_nonstationary);

current_rft = spm_get_defaults('stats.rft.nonstat');

if current_rft ~= cfg.rft_nonstationary

    warning( ...
        'SPM parameter "stats.rft.nonstat" was set to %f instead of %f.', ...
        current_rft, ...
        cfg.rft_nonstationary);
end

% -------------------------------------------------------------------------
% Intrinsic smoothness estimation parameter
% Recommended for quantitative MRI data
% -------------------------------------------------------------------------
spm_get_defaults('stats.fmri.ufp', cfg.ufp);

current_ufp = spm_get_defaults('stats.fmri.ufp');

if current_ufp ~= cfg.ufp

    warning( ...
        'SPM parameter "stats.fmri.ufp" was set to %f instead of %f.', ...
        current_ufp, ...
        cfg.ufp);
end

%% ========================================================================
% PREPARE REGRESSOR FILES
% =========================================================================
% Generate covariate files if they do not already exist.
% =========================================================================

aj_get_regfile(cfg.dataset_dir);

%% ========================================================================
% TSPOON FILENAME STANDARDIZATION
% =========================================================================
% Rename GTSPOON files if required for compatibility.
% =========================================================================

if strcmpi(cfg.smoothing_name, 'TSPOON')

    aj_get_gtspoon_rename(cfg.smoothing_dir);
end

%% ========================================================================
% START PARALLEL POOL
% =========================================================================

if isempty(gcp('nocreate'))
    parpool;
end

%% ========================================================================
% MAIN PROCESSING LOOP
% =========================================================================

for iMetric = 1:numel(cfg.qmetrics)

    metric = cfg.qmetrics{iMetric};

    for iTissue = 1:numel(cfg.tissues)

        tissue = cfg.tissues{iTissue};

        fprintf('\n====================================================\n');
        fprintf('Processing: %s | %s\n', metric, tissue);
        fprintf('====================================================\n');

        %% ----------------------------------------------------------------
        % Output directory
        % -----------------------------------------------------------------
        outdir = fullfile( ...
            cfg.dataset_dir, ...
            'derivatives', ...
            sprintf( ...
            '%s_GLM_%s_age_nonlinear', ...
            cfg.smoothing_name, ...
            cfg.rft_name), ...
            sprintf('%s_%s', metric, tissue));

        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end

        %% ----------------------------------------------------------------
        % Input qMRI maps
        % -----------------------------------------------------------------
        scans = aj_get_smoofiles_selection( ...
            cfg.smoothing_index, ...
            cfg.smoothing_dir, ...
            cfg.smoothing_methods, ...
            metric, ...
            tissue);

        if isempty(scans) || all(cellfun(@isempty, scans))

            warning( ...
                'No files found for %s | %s.', ...
                metric, ...
                tissue);

            continue;
        end

        %% ----------------------------------------------------------------
        % Explicit tissue mask
        % -----------------------------------------------------------------
        mask = cellstr( ...
            spm_select( ...
            'ExtFPList', ...
            fullfile(cfg.dataset_dir, 'derivatives'), ...
            sprintf('^atlas-.*%s_space-MNI_mask.*\\.nii$', tissue)));

        if isempty(mask)

            warning('No mask found for tissue: %s', tissue);

            continue;
        end

        %% ----------------------------------------------------------------
        % Run SPM batch
        % -----------------------------------------------------------------
        inputs = cell(3,1);

        inputs{1} = cellstr(outdir);
        inputs{2} = cellstr(scans);
        inputs{3} = mask;

        spm_jobman('run', cfg.jobfile, inputs{:});
    end
end

%% ========================================================================
% CLOSE PARALLEL POOL
% =========================================================================

try
    delete(gcp('nocreate'));

catch
    warning('Parallel pool already closed or unavailable.');
end

fprintf('\nParametric GLM pipeline completed.\n');