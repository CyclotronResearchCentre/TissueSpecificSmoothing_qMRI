%--------------------------------------------------------------------------
% Main script to create population masks:
% - ICV mask over the population
% - WTA mask over the population
% - ICV-WTA mask over the population
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Create mask mean ICV
close all; clear;clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

% Load the NIfTI file
V = spm_vol('D:\Master_Thesis\Data\BIDS_AgingData\derivatives\atlas-MTsat_space-MNI_res-high_desc-meanICV.nii');
data = spm_read_vols(V);

% Binarize the data
binary_data = data ~= 0;

% Create a new NIfTI structure for saving
V.fname = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\atlas-MTsat_space-MNI_res-high_desc-maskMeanICV.nii';
spm_write_vol(V, binary_data);

%% Compute the "winner takes all" masks over all the population
ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
MNI_dir = fullfile(ds_dir,'derivatives','SPM12_dartel');

TCs = {'GM', 'WM', 'CSF'};

GM_density_pattern = sprintf('^.*label-%s_probseg\\.nii$', TCs{1});
WM_density_pattern = sprintf('^.*label-%s_probseg\\.nii$', TCs{2});
CSF_density_pattern = sprintf('^.*label-%s_probseg\\.nii$', TCs{3});

GM_density_paths = cellstr(spm_select('FPListRec', MNI_dir, GM_density_pattern));
WM_density_paths = cellstr(spm_select('FPListRec', MNI_dir, WM_density_pattern));
CSF_density_paths = cellstr(spm_select('FPListRec', MNI_dir, CSF_density_pattern));

pth_out = fullfile(ds_dir,'derivatives');

WTAM_paths = aj_compute_WTAM(GM_density_paths, WM_density_paths, CSF_density_paths, pth_out);

%% Apply maskMeanICV to WTAMasks
pth_out = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives';

maskMeanICV_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\atlas-MTsat_space-MNI_res-high_desc-maskMeanICV.nii';
WTA_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\WTA_CSF_mask.nii';

Vi = char(maskMeanICV_path, WTA_path);
mask_WTA_path = spm_file(WTA_path,'prefix','masked_','path',pth_out);
f = '(i1.*i2)';
% Flags for image calculation
ic_flag = struct(... % type is set below based on that of input image
    'interp', -4);   % 4th order spline interpolation
spm_imcalc(Vi, mask_WTA_path, f, ic_flag);