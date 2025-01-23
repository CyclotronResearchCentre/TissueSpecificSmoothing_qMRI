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
%% Create regressors files
% Cleaning environment & setting up SPM path
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

% Paths to access to the data
ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
tsv_path = fullfile(ds_dir,'participants.tsv');

% Read the .tsv file
opts = detectImportOptions(tsv_path, 'FileType', 'text');
data = readtable(tsv_path, opts);

% Convert categorical variables (sex, scanner) to boolean
sex = double(strcmp(data.sex, 'M')); % Male=1, Female=0
scanner = double(strcmp(data.scanner, 'trio')); % Trio=1, Quatro=0

age = data.age;
TIV = data.TIV;         % Total Intracranial Volume

save(fullfile(ds_dir, 'reg_age.mat'), 'age');
save(fullfile(ds_dir, 'reg_sex.mat'), 'sex');
save(fullfile(ds_dir, 'reg_TIV.mat'), 'TIV');
save(fullfile(ds_dir, 'reg_scanner.mat'), 'scanner');

%% TSPOON Renaming
clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

param.smoothName = 'AJ-TSPOON';
ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
bids_dir = fullfile(ds_dir, 'derivatives', param.smoothName);

% Get list of GM files
GM_files = spm_select('FPListRec', bids_dir, sprintf('^tspoon1_sub*'));
nfiles = size(GM_files, 1);

% Rename each file
for i = 1:nfiles
    original_file = strtrim(GM_files(i, :));
    original_basename = spm_file(original_file, 'basename');
    new_basename = regexprep(original_basename, 'tspoon1_', 'TSPOON_GM_');
    new_file = spm_file(original_file, 'basename', new_basename);
    movefile(original_file, new_file);
    
    fprintf('Renamed: %s -> %s\n', original_file, new_file);
end

% Get list of WM files
WM_files = spm_select('FPListRec', bids_dir, sprintf('^tspoon2_sub*'));
nfiles = size(WM_files, 1);

% Rename each file
for i = 1:nfiles
    original_file = strtrim(WM_files(i, :));
    original_basename = spm_file(original_file, 'basename');
    new_basename = regexprep(original_basename, 'tspoon2_', 'TSPOON_WM_');
    new_file = spm_file(original_file, 'basename', new_basename);
    movefile(original_file, new_file);
    
    fprintf('Renamed: %s -> %s\n', original_file, new_file);
end

%% GLM on smoothed data
% Cleaning environment & setting up SPM path
clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

% Choose 'TSPOON', 'TWS' or 'TWsmoot'
smoothing_method = {'TWsmoot', 'TWS', 'TSPOON'};
meth = 3;

% Paths to access to script, smoothed data and data
script_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\agingdata\reprod_article';
ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
if meth==1
    param.smoothName = smoothing_method{meth};
else
    param.smoothName = sprintf('AJ-%s', smoothing_method{meth});
end
bids_dir = fullfile(ds_dir,'derivatives',param.smoothName);

% Set up the metrics and TC names lists
metrics_names = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
nMetricsNames = length(metrics_names);
TC_names = {'GM', 'WM'};
nTCNames = length(TC_names);
% Initialize SPM configuration
spm_jobman('initcfg'); 
spm('defaults', 'fmri');

% Start parallel pool if not already open
if isempty(gcp('nocreate'))
    parpool; % Create a default parallel pool
end

parfor i = 1:nMetricsNames
    smoothing_method = {'TWsmoot', 'TWS', 'TSPOON'};
    metrics_names = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
    TC_names = {'GM', 'WM'};
    inputs  = cell(3,1);
    jobfile = cellstr(fullfile(script_dir,'aj_batch_GLM_job.m'));
    
    for ii = 1:nTCNames
        % Output file
        inputs{1} = cellstr(fullfile(ds_dir,'derivatives',sprintf('AJ-%s_GLM',smoothing_method{meth}),sprintf('%s_%s',metrics_names{i},TC_names{ii})));
        
        % Select files for each qMRI metric and each TC for all subjects
        if meth==1
            inputs{2} = cellstr(spm_select('FPListRec',bids_dir,sprintf('^.*%ssmo_%s.*\\.nii$',TC_names{ii},metrics_names{i})));
        else
            inputs{2} = cellstr(spm_select('FPListRec',bids_dir,sprintf('^%s_%s.*%s.*\\.nii$',smoothing_method{meth},TC_names{ii},metrics_names{i}))); 
        end
         
        % Select explicit masks defining GM and WM voxels
        inputs{3} = cellstr(spm_select('ExtFPList',fullfile(ds_dir,'derivatives'),sprintf('^atlas-.*%s_space-MNI_mask.*\\.nii$',TC_names{ii})));
        
        spm_jobman('run', jobfile, inputs{:});
    end
end

delete(gcp('nocreate'));

%% CORRECTED GLM on smoothed data
% WARNING: have to remove % in job file to do spm_get_defaults('stats.fmri.ufp',0.5);

% New script to run the initial batch for a specific qMRI parameter and
% tissue class with spm.stats.fmri.ufp = 0.5 instead of 0.001 -> test if
% the generated mask.nii during model estimation corresponds better.

% Cleaning environment & setting up SPM environment
clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');
spm_get_defaults('stats.fmri.ufp',0.5);

% Choose 'TSPOON', 'TWS' or 'TWsmoot'
smoothing_method = {'TWsmoot', 'TWS', 'TSPOON'};
meth = 2;

% Paths to access to script, smoothed data and data
script_dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\agingdata';
ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
if meth==1
    param.smoothName = smoothing_method{meth};
else
    param.smoothName = sprintf('AJ-%s', smoothing_method{meth});
end
bids_dir = fullfile(ds_dir,'derivatives',param.smoothName);

% Set up the metrics and TC names lists
metrics_names = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
nMetricsNames = length(metrics_names);
TC_names = {'GM', 'WM'};
nTCNames = length(TC_names);

% Initialize SPM configuration
spm_jobman('initcfg'); 
spm('defaults', 'fmri');

% Start parallel pool if not already open
% if isempty(gcp('nocreate'))
%     parpool; % Create a default parallel pool
% end

for i = 4:4 %1:nMetricsNames
    smoothing_method = {'TWsmoot', 'TWS', 'TSPOON'};
    metrics_names = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
    TC_names = {'GM', 'WM'};
    inputs  = cell(3,1);
    jobfile = cellstr(fullfile(script_dir,'aj_batch_GLM_job.m'));
    
    for ii = 1:1 %1:nTCNames
        % Output file
        inputs{1} = cellstr(fullfile(ds_dir,'derivatives',sprintf('AJ-%s_corrGLM',smoothing_method{meth}),sprintf('%s_%s',metrics_names{i},TC_names{ii})));
        % Create output directory if it doesn't exist
        if ~exist(char(inputs{1}), 'dir')
            mkdir(char(inputs{1}));
        end
        
        % Select files for each qMRI metric and each TC for all subjects
        if meth==1
            inputs{2} = cellstr(spm_select('FPListRec',bids_dir,sprintf('^.*%ssmo_%s.*\\.nii$',TC_names{ii},metrics_names{i})));
        else
            inputs{2} = cellstr(spm_select('FPListRec',bids_dir,sprintf('^%s_%s.*%s.*\\.nii$',smoothing_method{meth},TC_names{ii},metrics_names{i}))); 
        end
         
        % Select explicit masks defining GM and WM voxels
        inputs{3} = cellstr(spm_select('ExtFPList',fullfile(ds_dir,'derivatives'),sprintf('^atlas-.*%s_space-MNI_mask.*\\.nii$',TC_names{ii})));
        
        spm_jobman('run', jobfile, inputs{:});
        % WARNING: have to remove % in job file to do spm_get_defaults('stats.fmri.ufp',0.5);
    end
end

% delete(gcp('nocreate'));