function Y_thr = aj_get_thrMatrix(GLM_dir, flagBinary)
%--------------------------------------------------------------------------
% Function to get a (binary) thresholded 3D matrix from a NIfTI file path.
% This function loads statistical maps (e.g., spmT_*.nii or spmF_*.nii files) 
% from a specified GLM directory and thresholds them based on a specified 
% statistical threshold. The thresholding can be done either to create a 
% binary mask (using a fixed threshold) or to zero out values below the threshold 
% while keeping the rest unchanged.
%
% The function uses the Family-Wise Error (FWE) corrected threshold obtained 
% from the SPM structure (`SPM.mat`) and applies this threshold to the 
% statistical maps (e.g., T- or F-statistics).
%
% INPUTS:
% GLM_dir:  A string specifying the path to the GLM directory that contains 
%           the statistical maps (e.g., spmT_000x.nii files).
% flagBinary: A boolean flag (0 or 1) indicating the type of thresholding:
%   - 1 (binary thresholding): Any voxel with a value above the threshold 
%     is set to 1, while below-threshold voxels are set to 0.
%   - 0 (non-binary thresholding): Voxels below the threshold are set to 0, 
%     while other values are retained.
%
% OUTPUT
% Y_thr:    A cell array of thresholded 3D matrices (one per contrast). If 
%           `flagBinary` is set to 1, the output is a binary mask; if set 
%           to 0, the output retains non-zero values above the threshold.
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

% Load SPM.mat to set up the contrast
load(fullfile(GLM_dir, 'SPM.mat'), 'SPM');

%% Do the job
Y_thr = cell(nContrasts,1);
for i = 1:nContrasts 
    % Initialize the xSPM structure
    xSPM = struct();
    xSPM.swd = GLM_dir;
    xSPM.Ic = i;                % Contrast index
    xSPM.Im = [];               % Mask
    xSPM.Ex = 0;                % Inclusive mask flag
    xSPM.thresDesc = 'FWE';     % FWE-corrected
    xSPM.u = 0.05;              % p-value threshold
    xSPM.k = 0;                 % Cluster extent threshold

    % Generate xSPM structure
    [~, xSPM] = spm_getSPM(xSPM);

    % Compute the FWE threshold using spm_uc
    % Extract required parameters for spm_uc
    R = xSPM.R;      % Resel counts (search space)
    df = xSPM.df;    % Degrees of freedom
    STAT = xSPM.STAT;
    n = 1;           % Number of conjoint tests
    p_thresh = 0.05; % FWE-corrected p-value

    % Compute threshold
    u_fwe = spm_uc(p_thresh, df, STAT, R, n);

    % Threshold the statistical map
    Y = spm_read_vols(xSPM.Vspm);
    if flagBinary
        Y_thr{i} = Y > u_fwe;
    else
        Z = Y;
        Z(Z < u_fwe) = 0;
        Y_thr{i} = Z;
    end

    % Return to the initial directory at the end of the script
    cd(initial_dir);
end
end