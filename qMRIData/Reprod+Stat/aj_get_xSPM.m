function GLM_xSPM = aj_get_xSPM(GLM_dir)
% This function processes GLM results stored in a specified directory, 
% generating the `xSPM` structures needed for subsequent statistical 
% analysis. These structures are used to apply thresholds, handle 
% contrasts and manage other parameters for statistical maps.
% The function loads and processes all the T-statistics (`spmT_*.nii`) and 
% F-statistics (`spmF_*.nii`) files found in the given directory and 
% prepares the `xSPM` structure for each contrast.
%
% INPUTS
% GLM_dir:  A string specifying the path to the GLM directory that contains 
%           the statistical maps (e.g., `spmT_*.nii` or `spmF_*.nii`) and 
%           the associated `SPM.mat` file.
%
% OUTPUT
% GLM_xSPM: A cell array containing the `xSPM` structure for each contrast, 
%           which is used in the analysis and visualization of statistical 
%           maps in SPM.
%
% REFERENCE
% SPM documentation https://www.fil.ion.ucl.ac.uk/spm/doc/
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Dealing with inputs
if nargin<1
    error('aj_get_xSPM: Not enough inputs.');
end

if ~exist(GLM_dir, 'dir')
    error('aj_get_xSPM: the GLM directory does not exist.');
end

% Store the initial directory
initial_dir = pwd;

% Find all spmT_000x.nii and spmF_000x.nii files
file_list = dir(fullfile(GLM_dir, 'spmT_*.nii'));
file_list = [file_list; dir(fullfile(GLM_dir, 'spmF_*.nii'))];
nContrasts = length(file_list);

% Load the SPM.mat file associated with the results folder
load(fullfile(GLM_dir, 'SPM.mat'), 'SPM'); % Ensure SPM.mat exists

%% Do the job
GLM_xSPM = cell(1,nContrasts);
for i = 1:nContrasts    
    % Set up xSPM structure
    xSPM = struct();
    xSPM.swd = GLM_dir;             % Working directory
    xSPM.title = file_list(i).name; % Title from the file name
    xSPM.Ic = i;                    % Contrast index
    xSPM.Im = [];                   % Mask
    xSPM.Ex = 0;                    % Mask exclusion flag
    xSPM.thresDesc = 'FWE';         % Threshold description
    xSPM.u = 0.05;                  % Threshold (FWE-corrected p-value)
    xSPM.k = 0;                     % Extent threshold
    
    % Call spm_getSPM to populate the xSPM structure
    try
        [~, xSPM] = spm_getSPM(xSPM);
        GLM_xSPM{i} = xSPM; % Store the xSPM structure
    catch ME
        warning('Error processing contrast %d: %s', i, ME.message);
    end
end

% Return to the initial directory at the end of the script
cd(initial_dir);

end
