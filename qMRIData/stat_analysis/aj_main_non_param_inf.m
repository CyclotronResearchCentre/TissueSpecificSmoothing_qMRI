% =========================================================================
% SnPM MULTIPLE REGRESSION PIPELINE (Callaghan framework)
% =========================================================================
% This script performs a voxelwise non-parametric multiple regression
% analysis using SnPM, based on the methodology of Callaghan et al. (2014).
%
% The pipeline includes:
%   1. Model specification (SnPM MREG design)
%   2. Model computation (permutation-based inference setup)
%   3. Statistical inference (positive and negative contrasts)
%
% The analysis is performed across:
%   - Multiple qMRI metrics
%   - Tissue classes (GM / WM)
%   - Smoothing strategies
%
% =========================================================================
% REFERENCES
% Callaghan, M.F. et al. (2014)
% Neurobiology of Aging
% http://dx.doi.org/10.1016/j.neurobiolaging.2014.02.008
%
% BIDS Aging Dataset:
% https://github.com/CyclotronResearchCentre/BIDS_AgingData
% =========================================================================

%% ========================================================================
% INITIALIZATION
% =========================================================================
clear; clc; close all;

ds_dir = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives';

addpath("C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12");

spm_jobman('initcfg');
spm('defaults', 'FMRI');

%% ========================================================================
% USER PARAMETERS
% =========================================================================
smoo_approaches = {'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'};
qmetrics        = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
tcs             = {'GM', 'WM'};

%% ========================================================================
% PARALLEL SETUP
% =========================================================================
if isempty(gcp('nocreate'))
    parpool;
end

%% ========================================================================
% JOB FILES
% =========================================================================
job_mreg     = fullfile(pwd, 'batch_snpm_mreg_job.m');
job_compute  = fullfile(pwd, 'batch_snpm_compute_job.m');
job_infer_pos = fullfile(pwd, 'batch_snpm_inference_ar-pos_job.m');
job_infer_neg = fullfile(pwd, 'batch_snpm_inference_ar-neg_job.m');

%% ========================================================================
% 1. MODEL SPECIFICATION (SnPM MREG)
% =========================================================================
for smoo = 1:length(smoo_approaches)

    smoo_name = smoo_approaches{smoo};
    smoo_dir  = fullfile(ds_dir, smoo_name);

    % Optional preprocessing step (dataset-specific correction)
    if smoo == 3
        aj_get_gtspoon_rename(smoo_dir);
    end

    for q = 1:length(qmetrics)

        for tc = 1:length(tcs)

            outdir = fullfile( ...
                ds_dir, ...
                sprintf('TEST-SnPM-%s_mreg-age', smoo_name), ...
                sprintf('%s_%s', qmetrics{q}, tcs{tc}));

            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end

            scans = aj_get_smoofiles_selection(smoo, smoo_dir, smoo_approaches, ...
                                        qmetrics{q}, tcs{tc});

            if isempty(scans) || all(cellfun(@isempty, scans))
                warning('Missing input images: %s %s', tcs{tc}, qmetrics{q});
                continue;
            end

            mask = spm_select('ExtFPList', ds_dir, ...
                sprintf('^atlas-.*%s_space-MNI_mask.*\\.nii$', tcs{tc}));

            if isempty(mask)
                warning('Missing mask for %s', tcs{tc});
                continue;
            end

            inputs = {
                cellstr(outdir)
                scans
                cellstr(mask)
            };

            spm_jobman('run', job_mreg, inputs{:});
        end
    end
end

%% ========================================================================
% 2. MODEL COMPUTATION
% =========================================================================
for smoo = 1:length(smoo_approaches)

    for q = 1:length(qmetrics)

        for tc = 1:length(tcs)

            cfgfile = fullfile( ...
                ds_dir, ...
                sprintf('SnPM-%s_mreg-age', smoo_approaches{smoo}), ...
                sprintf('%s_%s', qmetrics{q}, tcs{tc}), ...
                'SnPMcfg.mat');

            if ~exist(cfgfile, 'file')
                warning('Missing SnPMcfg: %s', cfgfile);
                continue;
            end

            spm_jobman('run', job_compute, {cfgfile});
        end
    end
end

%% ========================================================================
% 3. STATISTICAL INFERENCE (POSITIVE & NEGATIVE CONTRASTS)
% =========================================================================
for smoo = 1:length(smoo_approaches)

    for q = 1:length(qmetrics)

        for tc = 1:length(tcs)

            snpm_file = fullfile( ...
                ds_dir, ...
                sprintf('SnPM-%s_mreg-age', smoo_approaches{smoo}), ...
                sprintf('%s_%s', qmetrics{q}, tcs{tc}), ...
                'SnPM.mat');

            if ~exist(snpm_file, 'file')
                warning('Missing SnPM result: %s', snpm_file);
                continue;
            end

            % ---------------------------------------------------------
            % Run BOTH contrasts explicitly (robust and reproducible)
            % ---------------------------------------------------------
            job_list = {job_infer_pos, job_infer_neg};

            for j = 1:2
                spm_jobman('run', job_list{j}, {snpm_file});
            end
        end
    end
end

%% ========================================================================
% CLEANUP
% =========================================================================
try
    delete(gcp('nocreate'));
catch
    warning('Parallel pool already closed or does not exist.');
end
