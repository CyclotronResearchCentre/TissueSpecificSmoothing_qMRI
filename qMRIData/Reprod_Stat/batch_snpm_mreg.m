% Script to:
% - run a multiple regression (mreg) designed from Callaghan et al. (2014)
% - use SnPM "MultiSub: Simple Regression; 1 covariate of interest"
% - use QUANTITAVE data
%--------------------------------------------------------------------------
% REFERENCES
% M.F. Callaghan and al. (2014) http://dx.doi.org/10.1016/j.neurobiolaging.2014.02.008
% https://hackmd.io/u_vOEzA8TzS1yGj52V6Txg
% https://github.com/CyclotronResearchCentre/BIDS_AgingData
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% User inputs
clear all; clc; close all;
% ds_dir = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives';
ds_dir = 'C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\data_sample\derivatives';

smoo_approaches = {'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'};
qmetrics = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
tcs = {'GM', 'WM'};

%% Script
% Initialize SPM configuration
addpath("C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12");
spm_jobman('initcfg');
spm('defaults', 'FMRI');

% Select the jobfile for GLM
script_dir = 'C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\TissueSpecificSmoothing\qMRIData\Reprod_Stat';
jobfile = {fullfile(script_dir, 'batch_snpm_mreg_job.m')};

% if isempty(gcp('nocreate')), parpool; end % Start a default parallel pool if not already open

for smoo = 3:3 %1:length(smoo_approaches)
    smoo_dir = fullfile(ds_dir, smoo_approaches{smoo});
    if smoo==3, aj_rename_gtspoon(smoo_dir); end % Function to rename GTSPOON files if needed
    
    for qmetric = 1:length(qmetrics)
        inputs  = cell(3,1);

        for tc = 1:length(tcs)
            % OUTPUT Folder
            inputs{1} = cellstr(fullfile(ds_dir,sprintf('SnPM-%s_mreg-age',smoo_approaches{smoo}),sprintf('%s_%s',qmetrics{qmetric},tcs{tc})));
            if ~exist(char(inputs{1}), 'dir'), mkdir(char(inputs{1})); end % Create output directory if it doesn't exist

            % INPUT Arguments: all nifti files for each qMRI metric and each TC for all subjects
            inputs{2} = aj_select_smoofiles(smoo, smoo_dir, smoo_approaches, qmetrics{qmetric}, tcs{tc});
            if isempty(inputs{2}) || all(cellfun(@isempty, inputs{2})), warning('No INPUT files found for %s %s', tcs{tc}, qmetrics{qmetric}); continue; end

            % Select explicit masks defining GM and WM voxels
            inputs{3} = cellstr(spm_select('ExtFPList',ds_dir,sprintf('^atlas-.*%s_space-MNI_mask.*\\.nii$',tcs{tc})));
            if isempty(inputs{3}), warning('No mask file found for %s', tcs{tc}); continue; end
            
            spm_jobman('run', jobfile, inputs{:});
        end
    end
end

jobfile = {fullfile(script_dir, 'batch_snpm_compute_job.m')};

for smoo = 3:3 %1:length(smoo_approaches)    
    for qmetric = 1:length(qmetrics)
        inputs  = cell(1,1);
        for tc = 1:length(tcs)
            % Compute
            inputs{1} = cellstr(fullfile(ds_dir,sprintf('SnPM-%s_mreg-age',smoo_approaches{smoo}),sprintf('%s_%s',qmetrics{qmetric},tcs{tc}), 'SnPMcfg.mat'));
            if isempty(inputs{1}) || all(cellfun(@isempty, inputs{1})), warning('No mreg file found for %s %s', tcs{tc}, qmetrics{qmetric}); continue; end
            
            spm_jobman('run', jobfile, inputs{:});
        end
    end
end

jobfile = {fullfile(script_dir, 'batch_snpm_inference_job.m')};

for smoo = 3:3 %1:length(smoo_approaches)    
    for qmetric = 1:length(qmetrics)
        inputs  = cell(1,1);
        for tc = 1:length(tcs)
            % Inefrence at voxel level
            inputs{1} = cellstr(fullfile(ds_dir,sprintf('SnPM-%s_mreg-age',smoo_approaches{smoo}),sprintf('%s_%s',qmetrics{qmetric},tcs{tc}), 'SnPM.mat'));
            if isempty(inputs{1}) || all(cellfun(@isempty, inputs{1})), warning('No results file from computed model found for %s %s', tcs{tc}, qmetrics{qmetric}); continue; end
            
            spm_jobman('run', jobfile, inputs{:});
        end
    end
end

% try, delete(gcp('nocreate')); catch, warning('Parallel pool already closed or does not exist.'); end
