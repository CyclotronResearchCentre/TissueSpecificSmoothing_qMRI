% Toy script to better understand TWS and TSPOON smoothing
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Histogramm of diff_GM vs. smwc1
clear;clc;
diff_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-Diff_TWSTSPOON_GLM\diff_GM_MTsat\spmF_0001.nii';
dir = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS';
GMseg_path = fullfile(dir, 'mean_sGMprobseg.nii');

diff = spm_read_vols(spm_vol(diff_path));
GMseg = spm_read_vols(spm_vol(GMseg_path));

valid_idx = ~isnan(diff) & diff ~= 0 & ~isnan(GMseg) & GMseg ~= 0;
GMseg_values = GMseg(valid_idx);
reverse_GMseg_values = 1./GMseg_values;
diff_values = diff(valid_idx);

figure;
scatter(reverse_GMseg_values, diff_values, 5, 'filled');
xlabel('GMseg Values ^{-1}');
ylabel('F Values'); % beta 1 and mask (not >5)
title('Scatter Plot of GMseg vs Diff');
grid on;

%% In Silico Simulation
close all;clear;clc;

GM_d = [1 1 1 1 1 1 1 1 0.75 0.75 0.6 0.5 0.25 0 0 0 0 0 0 0 0 0 0 0];

n = length(GM_d);
TC_d = zeros(3,n);
TC_d(1,:) = GM_d;
TC_d(2,:) = ones(1,n) - GM_d; % WM_d
TC_d(3,:) = zeros(1,n); % CSF_d

GM_s = zeros(1,n);
GM_s(:) = 20;
WM_s = zeros(1,n);
WM_s(:) = 10;

signal = TC_d(1,:).*GM_s + TC_d(2,:).*WM_s;

% Apply standard smoothing
% Y = filtfilt(B, A, X): The length of the input X must be more than three
% times the filter order, defined as max(length(B)-1,length(A)-1).
sm_kern = 4;
wg = gausswin(sm_kern);
wg = wg/sum(wg); % normalize
wg2 = gausswin(2*sm_kern); % double width smoothing
wg2 = wg2/sum(wg2); % normalize

gs_signal = filtfilt(wg,1,signal);
gs_TC_d = filtfilt(wg,1,TC_d')';

% Apply tissue-weighted smoothing
% Appplying the smoothing as implemented for VBQ,
% assuming the TPMs are like the tissue probability but smoothed with a
% kernel twice the size, for simplicity.
TWS_signal = zeros(2,n);
for ii=1:2
    tmp1 = signal .* TC_d(ii,:) .* (filtfilt(wg2,1,TC_d(ii,:))>.05); % Like the TPM masking
    tmp2 = gs_TC_d(ii,:) .* (gs_TC_d(ii,:)>.05); % masking from smoothed tissue
    TWS_signal(ii,:) = filtfilt(wg,1,tmp1) ./ tmp2;
end

% Explicit mask
% majority and above 20%
exMask = [  ...
    TC_d(1,:)>TC_d(2,:) & ... % GM>WM
    TC_d(1,:)>TC_d(3,:) & ... % GM>CSF
    TC_d(1,:)>.2 ; ...                % GM>.2
    TC_d(2,:)>TC_d(1,:) & ... % WM>GM
    TC_d(2,:)>TC_d(3,:) & ... % WM>CSF
    TC_d(2,:)>.2 ] ;                  % WM>.2

% Apply TSPOON smoothing
gs_exMask = filtfilt(wg,1,double(exMask)')';
TSPOON_signal = zeros(2,n);
for ii=1:2
    TSPOON_signal(ii,:) = filtfilt(wg,1,exMask(ii,:).*signal) ./ gs_exMask(ii,:);
end

figure;
hold on;
plot(signal, '-k', 'LineWidth', 1.5, 'DisplayName', 'signal');
plot(TWS_signal(1,:), '-b', 'LineWidth', 1.5, 'DisplayName', 'GM TWS signal');
plot(TWS_signal(2,:), '-c', 'LineWidth', 1.5, 'DisplayName', 'WM TWS signal');
plot(TSPOON_signal(1,:), '-r', 'LineWidth', 1.5, 'DisplayName', 'GM TSPOON signal');
plot(TSPOON_signal(2,:), '-y', 'LineWidth', 1.5, 'DisplayName', 'WM TSPOON signal');
legend('Location', 'best');
xlabel('Voxel Index');
ylabel('Intensity');
title('TWS + TSPOON');
ylim([0 25]);
grid on;
hold off;

%% Compute mean of smoothed TCprobseg for all TWS/TSPOON subjects
clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

dir = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS';
GMprobseg_paths = spm_select('FPListRec',dir,'^gs_sub-.*CSF_probseg.nii$');
data = zeros(181,217,181,138);
for i = 1:138
    data(:,:,:,i) = spm_read_vols(spm_vol(GMprobseg_paths(i,:)));
end
mean_data = mean(data, 4);

% Save the mean as a nifti file
mean_vol = spm_vol(GMprobseg_paths(1,:));
mean_vol.fname = fullfile(dir, 'mean_sCSFprobseg.nii');
spm_write_vol(mean_vol, mean_data);

%% Build binary thresholded mean of smoothed TCprobseg for all TWS subjects
clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

nTC = 3;
gs_mwc_paths = cell(1,nTC);
gs_mwc_paths{1} = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS\mean_sGMprobseg.nii';
gs_mwc_paths{2} = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS\mean_sWMprobseg.nii';
gs_mwc_paths{3} = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS\mean_sCSFprobseg.nii';

pth_out = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS';

for i = 1:nTC
    nii = spm_vol(gs_mwc_paths{i});
    data = spm_read_vols(nii);
    
    binary_data = double(data > 0.05 & ~isnan(data) & ~isinf(data));
        
    exMask_info = spm_vol(char(gs_mwc_paths{i})); % Use original info for output header
    exMask_info.fname = spm_file(exMask_info.fname,'prefix','thr_','path',pth_out);
    spm_write_vol(exMask_info, binary_data);
end

%% Estimation of the spatial areas difference
clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

% TWS_mask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS\thr_mean_sGMprobseg.nii';
% TSPOON_mask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TSPOON\thr_gs_Mask_mean_sGMprobseg.nii';

TWS_data_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS\sub-138\anat\TWS_GMw_sub-S138_space-MNI_MTsat.nii';
TSPOON_data_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TSPOON\sub-138\anat\TSPOON_GM_sub-S138_space-MNI_MTsat.nii';

TWS_data = spm_read_vols(spm_vol(TWS_data_path));
TSPOON_data = spm_read_vols(spm_vol(TSPOON_data_path));

TWS_mask = ~isnan(TWS_data);
TSPOON_mask = ~isnan(TSPOON_data);

volume_TWS = sum(TWS_mask(:));
volume_TSPOON = sum(TSPOON_mask(:));
fprintf('Volume couvert par TWS : %d voxels\n', volume_TWS);
fprintf('Volume couvert par TSPOON : %d voxels\n', volume_TSPOON);
volume_diff_TWS = sum(TWS_mask(:) & ~TSPOON_mask(:));
volume_diff_TSPOON = sum(TSPOON_mask(:) & ~TWS_mask(:));
fprintf('Volume non couvert par TSPOON : %d voxels, soit %2f pourcent du volume couvert.\n', volume_diff_TWS, volume_diff_TWS/volume_TWS*100);
fprintf('Volume non couvert par TWS : %d voxels, soit %2f pourcent du volume couvert.\n', volume_diff_TSPOON, volume_diff_TSPOON/volume_TSPOON*100);

mask_intersection = TWS_mask & TSPOON_mask;
TWS_common = TWS_data(mask_intersection);
TSPOON_common = TSPOON_data(mask_intersection);

flag.drawPlot = 1;
flag.savePlot = 0;
aj_BlandAltman(TWS_common, TSPOON_common, flag);

mean_diff = mean(TWS_common - TSPOON_common);
% relative_diff = mean((TWS_common - TSPOON_common) ./ TSPOON_common) * 100;
rmse = sqrt(mean((TWS_common - TSPOON_common).^2));
corr_coeff = corr(TWS_common, TSPOON_common);

fprintf('Différence moyenne : %.2f\n', mean_diff);
% fprintf('Différence relative moyenne : %.2f %%\n', relative_diff);
fprintf('RMSE : %.2f\n', rmse);
fprintf('Coefficient de corrélation : %.2f\n', corr_coeff);

%% Smoothing Quality
clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

TWS_data_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS\sub-138\anat\TWS_GMw_sub-S138_space-MNI_MTsat.nii';
TSPOON_data_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TSPOON\sub-138\anat\TSPOON_GM_sub-S138_space-MNI_MTsat.nii';

TWS_data = spm_read_vols(spm_vol(TWS_data_path));
TSPOON_data = spm_read_vols(spm_vol(TSPOON_data_path));

% NaN -> 0
TWS_data(isnan(TWS_data)) = 0;
TSPOON_data(isnan(TSPOON_data)) = 0;

% Assessment of spatial detail preservation
% Spatial Gradient: Measures the intensity of local changes. Excessive 
% preservation of detail could mean insufficient smoothing. The higher the 
% gradient, the more the image preserves fine spatial variations.
[Gx_TWS, Gy_TWS, Gz_TWS] = gradient(TWS_data);
TWS_gradient = sqrt(Gx_TWS.^2 + Gy_TWS.^2 + Gz_TWS.^2);

[Gx_TSPOON, Gy_TSPOON, Gz_TSPOON] = gradient(TSPOON_data);
TSPOON_gradient = sqrt(Gx_TSPOON.^2 + Gy_TSPOON.^2 + Gz_TSPOON.^2);

TWS_gradient_mean = mean(TWS_gradient(:));
TSPOON_gradient_mean = mean(TSPOON_gradient(:));

% Measure of global smoothness.
% Calculation of variance in local windows: A lower variance corresponds to
% a smoother image. Global average of local variances to quantify the 
% degree of smoothing.
win = ones(3, 3, 3);
TWS_local_variance = convn((TWS_data - mean(TWS_data(:))).^2, win, 'same');
TSPOON_local_variance = convn((TSPOON_data - mean(TSPOON_data(:))).^2, win, 'same');

TWS_smoothness = mean(TWS_local_variance(:));
TSPOON_smoothness = mean(TSPOON_local_variance(:));

% Image Contrast vs. Noise
% A measure of the signal-to-noise ratio (SNR). A higher SNR value 
% indicates better image quality.
TWS_active_mask = TWS_data ~= 0;
TWS_noise_mask = TWS_data == 0;
TWS_signal_mean = mean(TWS_data(TWS_active_mask));
TWS_noise_std = std(TWS_data(TWS_noise_mask));
TWS_SNR = TWS_signal_mean / TWS_noise_std;

TSPOON_active_mask = TSPOON_data ~= 0;
TSPOON_noise_mask = TSPOON_data == 0;
TSPOON_signal_mean = mean(TSPOON_data(TSPOON_active_mask));
TSPOON_noise_std = std(TSPOON_data(TSPOON_noise_mask));
TSPOON_SNR = TSPOON_signal_mean / TSPOON_noise_std;

% Inter-smoothing correlation (TWS vs TSPOON)
% Calculate the spatial correlation between the two smoothed images. A weak
% correlation may indicate that the two methods produce very different 
% results. Use a Pearson correlation coefficient.
TWS_vector = TWS_data(:);
TSPOON_vector = TSPOON_data(:);
correlation = corr(TWS_vector, TSPOON_vector);

% Intra-regional homogeneity. Compare variances or mean differences within
% defined anatomical regions.
TWS_region_variance = var(TWS_data);
TSPOON_region_variance = var(TSPOON_data);
