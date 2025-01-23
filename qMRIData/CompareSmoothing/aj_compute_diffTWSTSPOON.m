% Main script to compute the differences induced by the applied smoothing
% approach. First, the differences between TWS and TSPOON quantitative maps
% are computed for all qMRI parameters and tissues. Then a One-Sample
% T-Test (second level GLM) is applied to these substracted maps. A F-Test
% is applied to the intercept to get the regions affected by smoothing.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Create difference data: TWS - TSPOON 
clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

out_dir = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-Diff_TWSTSPOON';

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
smoothing_methods = {'TWS', 'TSPOON'};
TWS_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s',smoothing_methods{1}));
TSPOON_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s',smoothing_methods{2}));

qMRI_params = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
TCs = {'GM', 'WM'};

% Combining parameters using combvec (Cartesian product)
[TCs_idx, qMRI_idx] = ndgrid(1:numel(TCs), 1:numel(qMRI_params));

% Creating TWS and TSPOON patterns
TWS_patterns = arrayfun(@(i, j) sprintf('^%s_%s.*%s\\.nii$', ...
	smoothing_methods{1}, TCs{i}, qMRI_params{j}), ...
	TCs_idx(:), qMRI_idx(:), 'UniformOutput', false);

TSPOON_patterns = arrayfun(@(i, j) sprintf('^%s_%s.*%s\\.nii$', ...
	smoothing_methods{2}, TCs{i}, qMRI_params{j}), ...
	TCs_idx(:), qMRI_idx(:), 'UniformOutput', false);

TWS_paths = cellfun(@(pattern) spm_select('FPListRec', TWS_dir, pattern), ...
	TWS_patterns, 'UniformOutput', false);
TSPOON_paths = cellfun(@(pattern) spm_select('FPListRec', TSPOON_dir, pattern), ...
	TSPOON_patterns, 'UniformOutput', false);

% Do the job 
for i = 1:numel(TWS_paths)
    TWS_files = cellstr(TWS_paths{i});
    TSPOON_files = cellstr(TSPOON_paths{i});

    if numel(TWS_files) ~= numel(TSPOON_files)
        error('The number of TWS and TSPOON files does not match for the pattern %d.', i);
    end

    for ii = 1:numel(TWS_files)
        TWS_nii = spm_vol(TWS_files{ii});
        TSPOON_nii = spm_vol(TSPOON_files{ii});
        TWS_data = spm_read_vols(TWS_nii);
        TSPOON_data = spm_read_vols(TSPOON_nii);
        
        % Calculating the difference
        diff_data = TWS_data - TSPOON_data;
        
        % Build the new file name
        [~, file_name, ext] = fileparts(TSPOON_files{ii});
        new_name = regexprep(file_name, sprintf('^%s_', smoothing_methods{2}), 'diff_');
        sub_out_dir = fullfile(out_dir, sprintf('sub-%03d', ii));
        if ~exist(sub_out_dir, 'dir')
            mkdir(sub_out_dir);
        end
        output_file = fullfile(sub_out_dir, [new_name, ext]);
        
        diff_nii = TWS_nii;
        diff_nii.fname = output_file;
        spm_write_vol(diff_nii, diff_data);
    end
end

disp('All files have been processed and saved.');

%% Compute GLM of TWS-TSPOON difference
% Initialize SPM configuration
spm_jobman('initcfg'); 
spm('defaults', 'fmri');

qMRI_params = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
TCs = {'GM', 'WM'};

% Combining parameters using combvec (Cartesian product)
[TCs_idx, qMRI_idx] = ndgrid(1:numel(TCs), 1:numel(qMRI_params));

% Creating searching and output patterns
search_patterns = arrayfun(@(i, j) sprintf('^diff_%s.*%s\\.nii$', ...
	TCs{i}, qMRI_params{j}), ...
	TCs_idx(:), qMRI_idx(:), 'UniformOutput', false);
out_patterns = arrayfun(@(i, j) sprintf('diff_%s_%s', ...
	TCs{i}, qMRI_params{j}), ...
	TCs_idx(:), qMRI_idx(:), 'UniformOutput', false);
                    
ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
input_dir = fullfile(ds_dir,'derivatives','AJ-Diff_TWSTSPOON');
output_dir = fullfile(ds_dir,'derivatives','AJ-Diff_TWSTSPOON_GLM');

for i = 1:numel(search_patterns)
    jobfile = {'C:\Users\antoi\Documents\master_thesis\MATLAB\agingdata\compare_smoo\aj_GLM_diffTWSTSPOONGen_job.m'};
    
    inputs  = cell(2,1);
    % Output filename
    inputs{1} = cellstr(fullfile(output_dir,out_patterns{i}));
    % Input files
    inputs{2} = cellstr(spm_select('FPListRec', input_dir, search_patterns{i}));
    
    spm_jobman('run', jobfile, inputs{:});
end
