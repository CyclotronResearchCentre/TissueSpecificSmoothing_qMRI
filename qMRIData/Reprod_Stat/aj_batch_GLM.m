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
smoo_dir = fullfile(ds_dir, 'derivatives', param.smoothName);

% Get list of GM files
GM_files = spm_select('FPListRec', smoo_dir, sprintf('^tspoon1_sub*'));
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
WM_files = spm_select('FPListRec', smoo_dir, sprintf('^tspoon2_sub*'));
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
% % Cleaning environment & setting up SPM path
% clear;clc;
% % addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');
% addpath("C:\Users\aj\Documents\toolbox\spm12");
% 
% % Choose 'TSPOON', 'TWS' or 'TWsmoot'
% smoothing_method = {'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'};
% meth = 4;
% 
% % Paths to access to script, smoothed data and data
% ds_dir = 'C:\Users\aj\Documents\SMOOTHING\data\qMRI_AgingCallaghan';
% if meth==1
%     param.smoothName = smoothing_method{meth};
% else
%     param.smoothName = sprintf('AJ-%s', smoothing_method{meth});
% end
% bids_dir = fullfile(ds_dir,'derivatives',param.smoothName);
% 
% if meth==4
%     tmp_unzip_dir = fullfile(ds_dir,'derivatives','AJ-SUSAN_unzipped');
%     if ~exist(tmp_unzip_dir,'dir'), mkdir(tmp_unzip_dir); end
% end
% 
% % Set up the metrics and TC names lists
% metrics_names = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
% nMetricsNames = length(metrics_names);
% TC_names = {'GM', 'WM'};
% nTCNames = length(TC_names);
% 
% % Initialize SPM configuration
% spm_jobman('initcfg'); 
% spm('defaults', 'fmri');
% 
% % Start parallel pool if not already open
% if isempty(gcp('nocreate'))
%     parpool; % Create a default parallel pool
% end
% 
% parfor i = 4:4 %1:nMetricsNames
%     smoothing_method = {'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'};
%     metrics_names = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
%     TC_names = {'GM', 'WM'};
%     inputs  = cell(3,1);
%     jobfile = cellstr(fullfile(pwd,'aj_batch_GLM_job.m'));
% 
%     for ii = 2:2 %1:nTCNames
%         % Output file
%         inputs{1} = cellstr(fullfile(ds_dir,'derivatives',sprintf('AJ-%s_GLM',smoothing_method{meth}),sprintf('%s_%s',metrics_names{i},TC_names{ii})));
% 
%         % Select files for each qMRI metric and each TC for all subjects
%         if meth==1
%             inputs{2} = cellstr(spm_select('FPListRec',bids_dir,sprintf('^.*%ssmo_%s.*\\.nii$',TC_names{ii},metrics_names{i})));
%         else
%             inputs{2} = cellstr(spm_select('FPListRec',bids_dir,sprintf('^%s_%s.*%s.*\\.nii$',smoothing_method{meth},TC_names{ii},metrics_names{i}))); 
%         end
% 
%         % Select explicit masks defining GM and WM voxels
%         inputs{3} = cellstr(spm_select('ExtFPList',fullfile(ds_dir,'derivatives'),sprintf('^atlas-.*%s_space-MNI_mask.*\\.nii$',TC_names{ii})));
% 
%         spm_jobman('run', jobfile, inputs{:});
%     end
% end
% 
% delete(gcp('nocreate'));

%% GLM for QUANTITATIVE smoothed data
% Cleaning environment & setting up SPM environment
clear;clc;
addpath("C:\Users\lucad\Documents\smoothing\repo\spm12");
spm_get_defaults('stats.fmri.ufp',0.5);
spm_get_defaults('stats.rft.nonstat',1);

% Choose 'TSPOON', 'TWS' or 'TWsmoot'
smoothing_method = {'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'};
meth = 4;

% Paths to access to script, smoothed data and data
ds_dir = 'C:\Users\lucad\Documents\smoothing\data\qMRI_AgingCallaghan';
if meth==1
    param.smoothName = smoothing_method{meth};
else
    param.smoothName = sprintf('AJ-%s', smoothing_method{meth});
end
smoo_dir = fullfile(ds_dir,'derivatives',param.smoothName);

if meth==4
    tmp_unzip_dir = fullfile(ds_dir,'derivatives','AJ-SUSAN_unzipped');
    if ~exist(tmp_unzip_dir,'dir'), mkdir(tmp_unzip_dir); end
end

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

for i = 1:nMetricsNames
    inputs  = cell(3,1);
    jobfile = {fullfile('C:\Users\aj\Documents\SMOOTHING\TissueSpecificSmoothing\qMRIData\Reprod_Stat', 'aj_batch_GLM_job.m')};
    
    for ii = 1:nTCNames
        % Output file
        inputs{1} = cellstr(fullfile(ds_dir,'derivatives',sprintf('AJ-%s_GLM',smoothing_method{meth}),sprintf('%s_%s',metrics_names{i},TC_names{ii})));
        % Create output directory if it doesn't exist
        if ~exist(char(inputs{1}), 'dir')
            mkdir(char(inputs{1}));
        end
        
        % Select files for each qMRI metric and each TC for all subjects
        if meth==1
            inputs{2} = cellstr(spm_select('FPListRec',smoo_dir,sprintf('^.*%ssmo_%s.*\\.nii$',TC_names{ii},metrics_names{i})));
        elseif meth==4
            pattern = sprintf('^susan_Mask_%s.*%s.*\\.nii(\\.gz)?$', TC_names{ii}, metrics_names{i});
            file_list = cellstr(spm_select('FPListRec', smoo_dir, pattern));

            if isempty(file_list), warning('No file found for %s %s', TC_names{ii}, metrics_names{i}); continue; end

            % Unzip files if needed
            unzipped_files = cell(size(file_list));

            for f = 1:numel(file_list)
                [~,base,ext] = fileparts(file_list{f});

                if strcmp(ext,'.gz')
                    % full filename was .nii.gz → remove .gz
                    nii_out = fullfile(tmp_unzip_dir, base);
                    if ~exist(nii_out,'file')
                        fprintf('Unzipping %s\n', file_list{f});
                        gunzip(file_list{f}, tmp_unzip_dir);
                    end
                    unzipped_files{f} = nii_out;

                else
                    % already .nii
                    unzipped_files{f} = file_list{f};
                end
            end
            
            inputs{2} = unzipped_files;

        else
            inputs{2} = cellstr(spm_select('FPListRec',smoo_dir,sprintf('^%s_%s.*%s.*\\.nii$',smoothing_method{meth},TC_names{ii},metrics_names{i}))); 
        end

        % Select explicit masks defining GM and WM voxels
        inputs{3} = cellstr(spm_select('ExtFPList',fullfile(ds_dir,'derivatives'),sprintf('^atlas-.*%s_space-MNI_mask.*\\.nii$',TC_names{ii})));
        if isempty(inputs{3})
            warning('No mask file found for %s', TC_names{ii});
            continue;
        end

        spm_jobman('run', jobfile, inputs{:});
        % WARNING: have to remove % in job file to do spm_get_defaults('stats.fmri.ufp',0.5);
    end
end

try
    delete(gcp('nocreate'));
catch
    warning('Parallel pool already closed or does not exist.');
end
