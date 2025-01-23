function WTAM_paths = aj_compute_WTAM(GM_paths, WM_paths, CSF_paths, pth_out)
%--------------------------------------------------------------------------
% Function to create Winner-Takes-All (WTA) masks.
%
% INPUT
% GM_paths:     Path list of the GM probability maps
% WM_paths:     Path list of the WM probability maps
% CSF_paths:    Path list of the CSF probability maps
% pth_out:      Output path to save the resulting WTA masks
%
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
nsub = numel(GM_paths);

% Initialize matrices to store the average maps
V = spm_vol(GM_paths{1}); % Load the header of the first file
[dim1, dim2, dim3] = deal(V.dim(1), V.dim(2), V.dim(3)); % Dimensions

GM_mean = zeros(dim1, dim2, dim3);
WM_mean = zeros(dim1, dim2, dim3);
CSF_mean = zeros(dim1, dim2, dim3);

% Load and sum probability maps for each subject
for i = 1:nsub
    GM = spm_read_vols(spm_vol(GM_paths{i}));
    WM = spm_read_vols(spm_vol(WM_paths{i}));
    CSF = spm_read_vols(spm_vol(CSF_paths{i}));
    
    GM_mean = GM_mean + GM;
    WM_mean = WM_mean + WM;
    CSF_mean = CSF_mean + CSF;
end

% Divide by the number of subjects to get the average probability maps
GM_mean = GM_mean / nsub;
WM_mean = WM_mean / nsub;
CSF_mean = CSF_mean / nsub;

% Apply Winner Takes All (WTA) logic on the average maps
[~, max_idx] = max(cat(4, GM_mean, WM_mean, CSF_mean), [], 4);

% Generate separate binary masks for GM, WM, and CSF
GM_mask = (max_idx == 1); % GM
WM_mask = (max_idx == 2); % WM
CSF_mask = (max_idx == 3); % CSF

% Save the GM, WM, and CSF masks as separate NIfTI files
WTAM_paths = cell(1,3);

% Use the header from the original files for proper metadata
V.fname = fullfile(pth_out, 'WTA_GM_mask.nii');
V.dt = [2, 0]; % 2 for binary/int16
spm_write_vol(V, GM_mask);
WTAM_paths{1} = V.fname;

V.fname = fullfile(pth_out, 'WTA_WM_mask.nii');
spm_write_vol(V, WM_mask);
WTAM_paths{2} = V.fname;

V.fname = fullfile(pth_out, 'WTA_CSF_mask.nii');
spm_write_vol(V, CSF_mask);
WTAM_paths{3} = V.fname;

fprintf('WTA masks saved successfully for GM, WM, and CSF in %s.\n', pth_out);

end