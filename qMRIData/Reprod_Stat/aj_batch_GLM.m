% Script to process smoothed quantitative maps from published Callaghan
% dataset: AgingData (qMRI).
%
% First the script creates the regressors files (age, TIV, scanner and
% gender) before applying the GLM processing. This processing follows the
% processing described in the article. 
%
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
close all; clear all; clc;

% Choose between 'Perso PC', 'Desktop'
flags.users = 0;

% Choose your smoothing approach between 'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'
smoo = 4;

% Choose your random field theory value: 1 stationary hypothesis vs 2 non
% stationary hypothesis
hyp = 2;

% Quantitative MRI data ? 1 for yes and 0 for no
qdata = 0;

% USER WARNING: default spm values have to be changed in BOTH scripts (this
% one and in the corresponding jobfile)

if flags.users == 0
    addpath('C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\TissueSpecificSmoothing\qMRIData\Reprod_Stat');
end

% Get user paths, set up spm12 env
paths = aj_default_GLM(flags);

%% Prepare SPM default values for GLM
% Set spm stats.fmri.ufp default value to 0.5 for QUANTITATIVE MRI data
if qdata == 1
    spm_get_defaults('stats.fmri.ufp',0.5);
    current_ufp = spm_get_defaults('stats.fmri.ufp');
    if current_ufp ~= 0.5, warning('stats.fmri.ufp a été défini à %f au lieu de 0.5.', current_ufp); end
else
    current_ufp = spm_get_defaults('stats.fmri.ufp');
    if current_ufp ~= 0.001, warning('stats.fmri.ufp is equal to %f instead of 0.001 (spm default one).', current_ufp); end
    warning('Input GLM preparation is well done only for quantitative data.');
end

% Set spm stats.rft.nonstat default value to 1 if assuming non stationary smoothing
if hyp == 2
    spm_get_defaults('stats.rft.nonstat',1);
    current_rft = spm_get_defaults('stats.rft.nonstat');
    if current_rft ~= 1, warning('stats.rft.nonstat a été défini à %f au lieu de 1.', current_rft); end
else
    current_rft = spm_get_defaults('stats.rft.nonstat');
    if current_rft ~= 0, warning('stats.rft.nonstat is equal to %f instead of 0 (spm default one).', current_rft); end
end

%% Prepare input QUANTITATIVE smoothed data for GLM
% Set up the smoothing approaches, metrics and TC names lists
smoo_approachs = {'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'};
qmetrics = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
TCs = {'GM', 'WM'};

if smoo==1
    smoo_foldername = smoo_approachs{smoo};
else
    smoo_foldername = sprintf('AJ-%s', smoo_approachs{smoo});
    
    if smoo==4
        tmp_unzip_dir = fullfile(paths.ds_dir,'derivatives','AJ-SUSAN_unzipped');
        if ~exist(tmp_unzip_dir,'dir'), mkdir(tmp_unzip_dir); end
    end
end

smoo_dir = fullfile(paths.ds_dir,'derivatives', smoo_foldername);

if smoo==3
    % Function to rename GTSPOON files if needed
    aj_rename_gtspoon(smoo_dir);
end

%% GLM from Callaghan et al. (2014)
rft_hyps = {'rft0', 'rft1'};

% Create the regressor files otherwise prodive the existing ones
aj_compute_regfile(paths.ds_dir);

% Initialize SPM configuration
spm_jobman('initcfg');
spm('defaults', 'fmri');

% Select the jobfile for GLM
jobfile = {fullfile(paths.script_dir, 'aj_batch_GLM_job.m')};

% Start parallel pool if not already open
% if isempty(gcp('nocreate'))
%     parpool; % Create a default parallel pool
% end

for i = 4:4 %1:length(qmetrics)
    inputs  = cell(3,1);
    
    for ii = 1:1 %1:length(TCs)
        % OUTPUT Folder
        inputs{1} = cellstr(fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM_%s',smoo_approachs{smoo}, rft_hyps{hyp}),sprintf('%s_%s',qmetrics{i},TCs{ii})));
        
        % Create output directory if it doesn't exist
        if ~exist(char(inputs{1}), 'dir'), mkdir(char(inputs{1})); end
        
        % INPUT Arguments: all nifti files for each qMRI metric and each TC for all subjects
        inputs{2} = aj_select_smoofiles(tmp_unzip_dir, smoo, smoo_dir, smoo_approachs, qmetrics{i}, TCs{ii});
        if isempty(inputs{2}) || all(cellfun(@isempty, inputs{2})), warning('No file found for %s %s', TCs{ii}, qmetrics{i}); continue; end

        % Select explicit masks defining GM and WM voxels
        inputs{3} = cellstr(spm_select('ExtFPList',fullfile(paths.ds_dir,'derivatives'),sprintf('^atlas-.*%s_space-MNI_mask.*\\.nii$',TCs{ii})));
        if isempty(inputs{3}), warning('No mask file found for %s', TCs{ii}); continue; end

        spm_jobman('run', jobfile, inputs{:});
    end
end

% try
%     delete(gcp('nocreate'));
% catch
%     warning('Parallel pool already closed or does not exist.');
% end

%% IN WORK : Apply a FWE p-value of 0.05 on the GLM results & save it
smoo_approachs = {'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'};
rft_hyps = {'rft0', 'rft1'};
GLM_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM_%s',smoo_approachs{smoo}, rft_hyps{hyp}));
if ~exist(GLM_dir, 'dir'), warning('No GLM folder at %s', GLM_dir); end

aj_automate_spm_results(GLM_dir);
