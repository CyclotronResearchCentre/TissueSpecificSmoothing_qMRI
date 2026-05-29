% Main script to apply the different smoothing methods on 3D data
%
% This script apply Gaussian Smoothing (GS), Tissue-Weigthed Smoothing
% (TWS) and Tissue-SPecific smOOthing CompeNsated (TSPOON) to two datasets
% (fMRI and qMRI). Configurations and options are defined in the default
% fonction. Dataset path has to be in the BIDS standard.
%
%--------------------------------------------------------------------------
% fMRI Smoothing: Smoothing is applied to all the contrast images resulting
% from the First-Level Analysis GLM described in the article by Wakeman and
% al. (2015) in Section REFERENCES. Only the smoothing preprocessing step
% has been removed and the creation of the modulated warped tissue-class
% maps have been added.
% Here is a wrap up of resulting contrast files:
% ess0001 F-contrast named Canonical HRF effects of interest, which tests
% if there is any significant activation across the defined regressors.
% con_0002 is the contrast for Faces > Scrambled Faces.
% con_0003 is the contrast for Famous.
% con_0004 is the contrast for Unfamiliar.
% con_0005 is the contrast for Scrambled.
%--------------------------------------------------------------------------
% qMRI Smoothing: Smoothing is applied to the quantitative maps for all the
% four qMRI parameters (MTsat, PD, R1 and R2*) and use the tissue segmented
% maps.
%
%--------------------------------------------------------------------------
% REFERENCES
% BIDS: Gorgolewski, K.J., et al. (2016). BIDS: The brain imaging data structure. A standard for organizing and describing outputs of neuroimaging experiments. Scientific Data, 3, 160044.
% fMRI dataset: Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal human neuroimaging dataset. Sci. Data 2:150001 doi: 10.1038/sdata.2015.1
% fMRI dataset: Henson, R.N., Wakeman, D.G., Litvak, V. & Friston, K.J. (2011). A Parametric Empirical Bayesian framework for the EEG/MEG inverse problem: generative models for multisubject and multimodal integration. Frontiers in Human Neuroscience, 5, 76, 1-16.
% fMRI dataset: Chapter 42 of the SPM12 manual (http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf)
% qMRI dataset: M.F. Callaghan and al. (2014) http://dx.doi.org/10.1016/j.neurobiolaging.2014.02.008
% qMRI dataset: https://hackmd.io/u_vOEzA8TzS1yGj52V6Txg
% qMRI dataset: https://github.com/CyclotronResearchCentre/BIDS_AgingData
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% fMRI Dataset: Public OpenNeuro Dataset "ds000117 v6"
% Cleaning environment & setting up SPM path
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');
[param, flag] = aj_smooth_default();

% Paths to access to the data
work_bids_root = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117';
param.outDerivName = '1stat';
stat_path = fullfile(work_bids_root, 'derivatives', param.outDerivName);
preproc_path = fullfile(work_bids_root, 'derivatives', 'preprocessing');

% Use BIDS to get data information
BIDS_stat = spm_BIDS(stat_path);
BIDS_preproc = spm_BIDS(preproc_path);
nsub = length(BIDS_stat.subjects); % stat is the most limiting for number of subjects

% Looking for constrat MPMs and modulated warped TC for each subject
conMPM_paths = cell(nsub,1);
mwTC_paths = cell(nsub,1);
for i = 1:nsub
    conMPM_paths{i} = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'func'), '^.*con.*\.nii$');
    mwTC_paths{i} = spm_select('FPListRec', fullfile(BIDS_preproc.subjects(i).path, 'anat'), '^.*mwc.*\.nii$'); % well sorted
end

% Initialize cell arrays to store the smoothed results for all subjects
gs_imgaussfilt3_paths = cell(nsub,1);
gs_spm_paths = cell(nsub,1);
tws_paths = cell(nsub,1);
smwTC_paths = cell(nsub,1);
tspoon_paths = cell(nsub,1);
sMask_paths = cell(nsub,1);

% Start parallel pool if not already open
if isempty(gcp('nocreate'))
    parpool; % Create a default parallel pool
end

% Call smoothing functions for each subject with loaded NIfTI data
parfor i = 1:nsub
    fprintf('Executing smoothing for subject %d...\n', i);
    
    [gs_imgaussfilt3_paths{i}, gs_spm_paths{i},...  % Gaussian results
        tws_paths{i}, smwTC_paths{i},...            % TWS results
        tspoon_paths{i}, sMask_paths{i}] = ...      % TSPOON results
        aj_smoothing(conMPM_paths{i}, mwTC_paths{i}, param, flag, 3); % 3D so dim = 3
end

delete(gcp('nocreate'));

%% qMRI Dataset: Aging Data from M. Callaghan
% Cleaning environment & setting up SPM path
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

[param, flag] = aj_smooth_default();
param.fwhm_gs = 3;          % Kernel width for Gaussian Smoothing (GS)
param.fwhm_tws = 3;         % Kernel width for Tissue-Weighted Smoothing (TWS)
param.fwhm_tspoon = 3;      % Kernel width for Tissue-SPecific smOOthing compeNsated (TSPOON)

% Paths to access to the data
ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
param.outDerivName = 'SPM12_dartel';
BIDS_warped_data = fullfile(ds_dir, 'derivatives', param.outDerivName);

% Start parallel pool if not already open
if isempty(gcp('nocreate'))
    parpool; % Create a default parallel pool
end

% Use BIDS to get data information
BIDS_stat = spm_BIDS(BIDS_warped_data);
nsub = length(BIDS_stat.subjects);

% Set up the MPM names list
MPMs_listname = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
nMPMnames = length(MPMs_listname);

% Looking for warped MPMs and TC segmentation maps for each subject
wMPM_paths = cell(nsub,1);
TCseg_paths = cell(nsub,1);
parfor i = 1:nsub
    CSF_path_i = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'anat'), '^*CSF_probseg*\.nii$');
    GM_path_i = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'anat'), '^*GM_probseg*\.nii$');
    WM_path_i = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'anat'), '^*WM_probseg*\.nii$');
    TCseg_paths{i} = char(GM_path_i, WM_path_i, CSF_path_i); % well sorted
    
    for ii = 1:nMPMnames
        % Get the file list for current subject and MPM
        file_list = spm_select('FPListRec', fullfile(BIDS_stat.subjects(i).path,'anat'), ['^*' MPMs_listname{ii} '*\.nii$']);
        if ~isempty(file_list)
            % Ensure wMPM_paths{i} is a character array
            if isempty(wMPM_paths{i})
                wMPM_paths{i} = file_list; % Initialize with first file list
            else
                % Concatenate character arrays vertically
                wMPM_paths{i} = char(wMPM_paths{i}, file_list); 
            end
        end
    end
end

% Initialize cell arrays to store the smoothed results for all subjects
gs_imgaussfilt3_paths = cell(nsub,1);
gs_spm_paths = cell(nsub,1);
tws_paths = cell(nsub,1);
smwTC_paths = cell(nsub,1);
tspoon_paths = cell(nsub,1);
sMask_paths = cell(nsub,1);

% Parallel loop for smoothing
parfor i = 1:nsub
    fprintf('Executing smoothing for subject %d...\n', i);
    
    [gs_imgaussfilt3_paths{i}, gs_spm_paths{i},...  % Gaussian results
        tws_paths{i}, smwTC_paths{i},...            % TWS results
        tspoon_paths{i}, sMask_paths{i}] = ...      % TSPOON results
        aj_smoothing(wMPM_paths{i}, TCseg_paths{i}, param, flag, 3); % 3D so dim = 3
end

delete(gcp('nocreate'));
