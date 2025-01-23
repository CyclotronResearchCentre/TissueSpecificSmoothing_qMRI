% Toy script to check if the smoothed maps are well defined into the GM and
% WM masks. These maps are smoothed according to three different smoothing
% approaches: TWsmoot (originally smoothed using an old version of TWS, TWS
% (from hMRI toolbox) and TSPOON (generalized version).
%
% Conclusion: The old version of TWS (from VBQ toolbox) smooths maps based
% on subject-specific GM/WM masks. It doesn't use the TPM from SPM which
% are used as additional condition in the new version of TWS (greater than
% 5% of TPM)
%
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Check if all MPM TC-specific contrasts are inside the GLM TC-specific mask spatial area
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

smoothing_method = {'TWsmoot', 'AJ-TWS', 'AJ-TSPOON'};
meth = 1;
sub = 35; % [1, 138]
qMRI_param = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
param1 = 1;
param2 = 4;
TC_list = {'GM', 'WM'};
TC = 1;

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
sub_dir = fullfile(ds_dir, 'derivatives', smoothing_method{meth}, sprintf('sub-S%03d',sub), 'anat');
file1_path = spm_select('ExtFPList', sub_dir, sprintf('^.*%s.*%s\\.nii$', TC_list{TC},qMRI_param{param1}));
file2_path = spm_select('ExtFPList', sub_dir, sprintf('^.*%s.*%s\\.nii$', TC_list{TC},qMRI_param{param2}));

v1 = spm_read_vols(spm_vol(file1_path));
v2 = spm_read_vols(spm_vol(file2_path));

useful1 = (v1 ~= 0) & ~isnan(v1) & ~isinf(v1);
countUsefulVoxels = sum(useful1(:));
disp([smoothing_method{meth} '-' TC_list{TC} '-' qMRI_param{param1} ': Number of useful values: ', num2str(countUsefulVoxels)]);

useful2 = (v2 ~= 0) & ~isnan(v2) & ~isinf(v2);
countUsefulVoxels = sum(useful2(:));
disp([smoothing_method{meth} '-' TC_list{TC} '-' qMRI_param{param2} ': Number of useful values: ', num2str(countUsefulVoxels)]);

diffMask_12 = useful1 & ~useful2;
diffMask_21 = useful2 & ~useful1;

% Display the results
disp(['Subject:' num2str(sub) ', TC:' TC_list{TC} '- Number of differing values: ', num2str(nnz(diffMask_21) + nnz(diffMask_12))]);
disp(['Number of ' qMRI_param{param1} '!=0 and ' qMRI_param{param2} '=0: ', num2str(nnz(diffMask_12))]);
disp(['Number of ' qMRI_param{param2} '!=0 and ' qMRI_param{param1} '=0: ', num2str(nnz(diffMask_21))]);

%% Check if all MPM TC-specific contrasts are inside the GLM TC-specific mask spatial area
% close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
smoothing_method = {'TWsmoot', 'TWS', 'TSPOON'};
meth1 = 2;
meth2 = 3;

Callaghan_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s_GLM',smoothing_method{meth1}));
ourTWS_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s_GLM',smoothing_method{meth2}));

smoothing_method = {'origTWS', 'newTWS', 'TSPOON'};
contrasts = {'spmT_0001.nii', 'spmT_0002.nii', 'spmF_0003.nii'};
contrast_patterns = strcat('^', contrasts, '$');

Callaghan_lists = cellfun(@(pattern) spm_select('FPListRec', Callaghan_dir, pattern), ...
                              contrast_patterns, 'UniformOutput', false);
ourTWS_lists = cellfun(@(pattern) spm_select('FPListRec', ourTWS_dir, pattern), ...
                           contrast_patterns, 'UniformOutput', false);
Callaghan_list = vertcat(Callaghan_lists{:});
ourTWS_list = vertcat(ourTWS_lists{:});

% up to 24 -> odd numbers are GM
pathCall = Callaghan_list(8,:);
pathOur = ourTWS_list(8,:);
pathMask='D:\Master_Thesis\Data\BIDS_AgingData\derivatives\atlas-WM_space-MNI_mask.nii';

meth1Info = spm_vol(pathCall);
meth1Vol = spm_read_vols(meth1Info);
nVoxel_1 = meth1Info.dim(1) * meth1Info.dim(2) * meth1Info.dim(3);
nZero_1 = nVoxel_1 - nnz(meth1Vol);
tpm = isnan(meth1Vol);
nNan_1 = numel(tpm(tpm==1));
tpm = isinf(meth1Vol);
nInf_1 = numel(tpm(tpm==1));
nValues_1 = nVoxel_1-nZero_1-nNan_1-nInf_1;

meth2Info = spm_vol(pathOur);
meth2Vol = spm_read_vols(meth2Info);
nVoxel_2 = meth2Info.dim(1) * meth2Info.dim(2) * meth2Info.dim(3);
nZero_2 = nVoxel_2 - nnz(meth2Vol);
tpm = isnan(meth2Vol);
nNan_2 = numel(tpm(tpm==1));
tpm = isinf(meth2Vol);
nInf_2 = numel(tpm(tpm==1));
nValues_2 = nVoxel_2-nZero_2-nNan_2-nInf_2;

maskInfo = spm_vol(pathMask);
maskVol = spm_read_vols(maskInfo);
nVoxel_mask = maskInfo.dim(1) * maskInfo.dim(2) * maskInfo.dim(3);
nZero_mask = nVoxel_mask - nnz(maskVol);
tpm = isnan(maskVol);
nNan_mask = numel(tpm(tpm==1));
tpm = isinf(maskVol);
nInf_mask = numel(tpm(tpm==1));
nValues_mask = nVoxel_mask-nZero_mask-nNan_mask-nInf_mask;

% Calculate logical indices for mask and non-valid values
GMVol = meth1Vol; % GMCallVol GMOurVol
outOfMask = maskVol == 0; % Voxels outside the mask
nonZeroValues = GMVol ~= 0; % Non-zero values
nonNaNValues = ~isnan(GMVol); % Non-NaN values
nonInfValues = ~isinf(GMVol); % Non-Infinity values
validValuesOutsideMask = outOfMask & nonZeroValues & nonNaNValues & nonInfValues;
countValidOutsideMask = sum(validValuesOutsideMask(:));
disp([smoothing_method{meth1} ': Number of valid values outside the mask: ', num2str(countValidOutsideMask)]);

GMVol = meth2Vol; % GMCallVol GMOurVol
outOfMask = maskVol == 0; % Voxels outside the mask
nonZeroValues = GMVol ~= 0; % Non-zero values
nonNaNValues = ~isnan(GMVol); % Non-NaN values
nonInfValues = ~isinf(GMVol); % Non-Infinity values
validValuesOutsideMask = outOfMask & nonZeroValues & nonNaNValues & nonInfValues;
countValidOutsideMask = sum(validValuesOutsideMask(:));
disp([smoothing_method{meth2} ': Number of valid values outside the mask: ', num2str(countValidOutsideMask)]);

% Find the positions where GMOurVol has valid values but GMCallVol does not
% Logical masks for valid values
validMeth2 = (meth2Vol ~= 0) & ~isnan(meth2Vol) & ~isinf(meth2Vol);
validMeth1 = (meth1Vol ~= 0) & ~isnan(meth1Vol) & ~isinf(meth1Vol);

% Find the positions where meth2Vol is valid, but meth1Vol is not
diffMask_21 = validMeth2 & ~validMeth1;
linearIndicesDiff_21 = find(diffMask_21);
[subX_21, subY_21, subZ_21] = ind2sub(size(diffMask_21), linearIndicesDiff_21);

% Find the positions where meth1Vol is valid, but meth2Vol is not
diffMask_12 = validMeth1 & ~validMeth2;
linearIndicesDiff_12 = find(diffMask_12);
[subX_12, subY_12, subZ_12] = ind2sub(size(diffMask_12), linearIndicesDiff_12);

% Display the results
disp(['Number of differing values: ', num2str(nnz(diffMask_21) + nnz(diffMask_12))]);
disp(['Number of ' smoothing_method{meth2} '=0 and ' smoothing_method{meth1} '!=0: ', num2str(nnz(diffMask_12))]);
disp(['Number of ' smoothing_method{meth2} '!=0 and ' smoothing_method{meth1} '=0: ', num2str(nnz(diffMask_21))]);


%% Count the number of voxels with a value among all subjects (MPM-specific)
% Sum the binary qMRI MPM specific smoothed maps for all the 138 subjects
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
smoothing_method = {'TWsmoot', 'TWS', 'TSPOON'};
nS = length(smoothing_method);

% Loop through smoothing methods
binaryAllSubMatrix = cell(nS,1);
missingContributors = cell(nS,1);

for meth = 1:3 %1:nS
    if meth==1
        sub_dir = fullfile(ds_dir,'derivatives',sprintf('%s',smoothing_method{meth}));
        file_paths = spm_select('FPListRec', sub_dir, sprintf('^sub.*GMsmo_MTsat\\.nii$'));
    elseif meth==2
        sub_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s',smoothing_method{meth}));
        file_paths = spm_select('FPListRec', sub_dir, sprintf('^TWS_GMw.*MTsat\\.nii$'));
    elseif meth==3
        sub_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s',smoothing_method{meth}));
        file_paths = spm_select('FPListRec', sub_dir, sprintf('^TSPOON_GM.*MTsat\\.nii$'));
    end

    % Initialize variables
    sumBinaryMatrix = [];
    contributionTracker = []; % 4D matrix to track contributions

    for i = 1:size(file_paths, 1)
        % Read the volume
        matrix = spm_read_vols(spm_vol(file_paths(i, :)));

        % Binarize the matrix
        if meth == 1
            binaryMatrix = double(matrix ~= 0);
        else
            binaryMatrix = double(~isnan(matrix));
        end

        % Initialize matrices during the first iteration
        if i == 1
            sumBinaryMatrix = zeros(size(binaryMatrix));
            contributionTracker = zeros([size(binaryMatrix), size(file_paths, 1)]);
        end

        % Accumulate the binary matrices
        sumBinaryMatrix = sumBinaryMatrix + binaryMatrix;

        % Track contributions for this subject
        contributionTracker(:,:,:,i) = binaryMatrix;
    end

    % Binarize the summed matrix for voxels with value 138
    binaryAllSubMatrix{meth} = double(sumBinaryMatrix == 138);

    % Find indices of voxels with value 137, 136, or 135
    targetMask = (sumBinaryMatrix == 137); % | (sumBinaryMatrix == 136) | (sumBinaryMatrix == 135);
    [subX, subY, subZ] = ind2sub(size(sumBinaryMatrix), find(targetMask));

    % Identify missing contributors for these voxels
    missingContributors{meth} = cell(numel(subX), 1); % Each entry corresponds to a voxel
    for v = 1:numel(subX)
        x = subX(v); y = subY(v); z = subZ(v);

        % Find which subjects did NOT contribute at this voxel
        contributingSubjects = squeeze(contributionTracker(x, y, z, :));
        missingContributors{meth}{v} = find(contributingSubjects == 0); % List missing contributors
    end

    % Display information
    disp([smoothing_method{meth} ': Number of non zero elements: ' num2str(nnz(binaryAllSubMatrix{meth}))...
        ' over ' num2str(size(binaryAllSubMatrix{meth},1)*size(binaryAllSubMatrix{meth},2)*size(binaryAllSubMatrix{meth},3))...
        ' elements.']);
    disp([smoothing_method{meth} ': Missing contributors identified for target voxels.']);
end

%% Plot + subject missing counts (MPM-specific)
% Loop through the missingContributors for each smoothing method
subject_redundancy = cell(nS,1);
for meth = 3:3 %1:nS
    plotM = binaryAllSubMatrix{meth};
    sliceIndex = round(size(plotM,3)/2);
    slice = plotM(:, :, sliceIndex);
    figure;
    imagesc(slice);
    colormap(gray);
    colorbar;
    axis image;
    title(['Slice ', num2str(sliceIndex)]);
    xlabel('X-axis');
    ylabel('Y-axis');

    % Initialize a map to track redundancy for missing subjects
    missingCounts = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    
    % For each voxel in the current method's missingContributors
    for v = 1:numel(missingContributors{meth})
        % Get the list of missing subjects for this voxel
        missingSubjects = missingContributors{meth}{v};
        
        % Update the count for each missing subject
        for subj = missingSubjects
            if isKey(missingCounts, subj)
                missingCounts(subj) = missingCounts(subj) + 1;
            else
                missingCounts(subj) = 1;
            end
        end
    end
    
    % Store results
    subject_redundancy{meth} = missingCounts;
    
    % Convert the results into a sorted list
    missingSubjectsList = keys(missingCounts);
    redundancyCounts = values(missingCounts);

    % Combine into a matrix for easier display
    results = [cell2mat(missingSubjectsList(:)), cell2mat(redundancyCounts(:))];
    results = sortrows(results, -2); % Sort by redundancy count in descending order

    % Display the results
    disp('Subject redundancy count (Subject ID, Redundancy):');
    disp(results);
    % a=values(subject_redundancy{3}); a{138}
end

%% Count the number of voxels with a value among all subjects (gs_TC)
% Sum the binary qMRI MPM specific smoothed maps for all the 138 subjects
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
smoothing_method = {'TWS', 'TSPOON'};
nS = length(smoothing_method);
TC_names = {'GM', 'WM'};
nTC = length(TC_names);

for meth = 2:2 %1:nS
    if meth==1
        sub_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s',smoothing_method{meth}));
        GM_paths = spm_select('FPListRec', sub_dir, sprintf('^gs_sub.*%s_probseg\\.nii$',TC_names{1}));
        WM_paths = spm_select('FPListRec', sub_dir, sprintf('^gs_sub.*%s_probseg\\.nii$',TC_names{2}));
    elseif meth==2
        sub_dir = fullfile(ds_dir,'derivatives',sprintf('AJ-%s',smoothing_method{meth}));
        GM_paths = spm_select('FPListRec', sub_dir, sprintf('^gs_sub.*%s_probseg\\.nii$',TC_names{1}));
        WM_paths = spm_select('FPListRec', sub_dir, sprintf('^gs_sub.*%s_probseg\\.nii$',TC_names{2}));
    
%         % Combine paths into a cell array
%         fn_smwTC = {GM_paths; WM_paths};
%         
%         % Options for mask creation
%         opts = struct('minTCp', 0.2, 'noOvl', true, 'outPth', fullfile(sub_dir, 'hmri_masks'));
% 
%         % Call the function
%         [fn_maskTC, fn_meanTC] = hmri_proc_crtMask(fn_smwTC, opts);
%         
%         % Display output filenames
%         disp('Mask filenames:');
%         disp(fn_maskTC);
%         disp('Mean tissue class filenames:');
%         disp(fn_meanTC);
    end
end

GM_binaryAllSubMatrix = cell(nS,1);
WM_binaryAllSubMatrix = cell(nS,1);
for meth = 2:2 %1:nS
    % Initialize a matrix to hold the sum of binary matrices
    GM_sumBinaryMatrix = [];
    WM_sumBinaryMatrix = [];

    % Loop through each file
    for i = 1:size(GM_paths, 1)
        GM_matrix = spm_read_vols(spm_vol(GM_paths(i, :)));
        WM_matrix = spm_read_vols(spm_vol(WM_paths(i, :)));

        % Binarize the matrix
        GM_binaryMatrix = double(GM_matrix > 0.05);
        WM_binaryMatrix = double(WM_matrix > 0.05);

        % Initialize sumBinaryMatrix during the first iteration
        if i == 1
            GM_sumBinaryMatrix = zeros(size(GM_binaryMatrix));
            WM_sumBinaryMatrix = zeros(size(WM_binaryMatrix));
        end

        % Accumulate the binary matrices
        GM_sumBinaryMatrix = GM_sumBinaryMatrix + GM_binaryMatrix;
        WM_sumBinaryMatrix = WM_sumBinaryMatrix + WM_binaryMatrix;
    end

    % Binarize the summed matrix: 1 if value equals 138, otherwise 0
    GM_binaryAllSubMatrix{i} = double(GM_sumBinaryMatrix == 138);
    WM_binaryAllSubMatrix{i} = double(WM_sumBinaryMatrix == 138);
    GM_linearIndices138 = find(GM_binaryAllSubMatrix{i});
    [GM_subX, GM_subY, GM_subZ] = ind2sub(size(GM_binaryAllSubMatrix{i}), GM_linearIndices138);
    WM_linearIndices138 = find(WM_binaryAllSubMatrix{i});
    [WM_subX, WM_subY, WM_subZ] = ind2sub(size(WM_binaryAllSubMatrix{i}), WM_linearIndices138);

    disp([smoothing_method{meth} '-GM: Number of non zero elements: ' num2str(nnz(GM_binaryAllSubMatrix{i}))...
        ' over ' num2str(size(GM_binaryAllSubMatrix{i},1)*size(GM_binaryAllSubMatrix{i},2)*size(GM_binaryAllSubMatrix{i},3))...
        ' elements.']);
    disp([smoothing_method{meth} '-WM: Number of non zero elements: ' num2str(nnz(WM_binaryAllSubMatrix{i}))...
        ' over ' num2str(size(WM_binaryAllSubMatrix{i},1)*size(WM_binaryAllSubMatrix{i},2)*size(WM_binaryAllSubMatrix{i},3))...
        ' elements.']);
end

%% Compare TSPOON hmri masks with Callaghan mean masks
clear;clc;
% both pairs are already binerazised (mask)
GM_origMask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\atlas-GM_space-MNI_mask.nii';
WM_origMask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\atlas-WM_space-MNI_mask.nii';
GM_hmriMask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TSPOON\hmri_masks\mask_gs_sub-S001_MTsat_space-MNI_desc-mod_label-GM_probseg.nii';
WM_hmriMask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TSPOON\hmri_masks\mask_gs_sub-S001_MTsat_space-MNI_desc-mod_label-WM_probseg.nii';

GM_origMask = spm_read_vols(spm_vol(GM_origMask_path));
WM_origMask = spm_read_vols(spm_vol(WM_origMask_path));
GM_hmriMask = spm_read_vols(spm_vol(GM_hmriMask_path));
WM_hmriMask = spm_read_vols(spm_vol(WM_hmriMask_path));

GMdiffMask_21 = GM_origMask & ~GM_hmriMask;
GMlinearIndicesDiff_21 = find(GMdiffMask_21);
[GMsubX_21, GMsubY_21, GMsubZ_21] = ind2sub(size(GMdiffMask_21), GMlinearIndicesDiff_21);

GMdiffMask_12 = GM_hmriMask & ~GM_origMask;
GMlinearIndicesDiff_12 = find(GMdiffMask_12);
[GMsubX_12, GMsubY_12, GMsubZ_12] = ind2sub(size(GMdiffMask_12), GMlinearIndicesDiff_12);

disp(['Number of differing values: ', num2str(nnz(GMdiffMask_21) + nnz(GMdiffMask_12))]);
disp(['Number of GM_hmriMask & ~GM_origMask: ', num2str(nnz(GMdiffMask_12))]);
disp(['Number of GM_origMask & ~GM_hmriMask: ', num2str(nnz(GMdiffMask_21))]);

WMdiffMask_21 = WM_origMask & ~WM_hmriMask;
WMlinearIndicesDiff_21 = find(WMdiffMask_21);
[WMsubX_21, WMsubY_21, WMsubZ_21] = ind2sub(size(WMdiffMask_21), WMlinearIndicesDiff_21);

WMdiffMask_12 = WM_hmriMask & ~WM_origMask;
WMlinearIndicesDiff_12 = find(WMdiffMask_12);
[WMsubX_12, WMsubY_12, WMsubZ_12] = ind2sub(size(WMdiffMask_12), WMlinearIndicesDiff_12);

disp(['Number of differing values: ', num2str(nnz(WMdiffMask_21) + nnz(WMdiffMask_12))]);
disp(['Number of GM_hmriMask & ~GM_origMask: ', num2str(nnz(WMdiffMask_12))]);
disp(['Number of GM_origMask & ~GM_hmriMask: ', num2str(nnz(WMdiffMask_21))]);

%% DONE! Create hmri masks based on NON smoothed data
close all; clear; clc;
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\spm12');

ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
TC_names = {'GM', 'WM'};
nTC = length(TC_names);

sub_dir = fullfile(ds_dir,'derivatives','SPM12_dartel');
GM_paths = spm_select('FPListRec', sub_dir, sprintf('^sub.*%s_probseg\\.nii$',TC_names{1}));
WM_paths = spm_select('FPListRec', sub_dir, sprintf('^sub.*%s_probseg\\.nii$',TC_names{2}));

% Combine paths into a cell array
fn_smwTC = {GM_paths; WM_paths};

% Options for mask creation
opts = struct('minTCp', 0.2, 'noOvl', true, 'outPth', fullfile(sub_dir, 'hmri_masks'));

% Call the function
[fn_maskTC, fn_meanTC] = hmri_proc_crtMask(fn_smwTC, opts);

% Display output filenames
disp('Mask filenames:');
disp(fn_maskTC);
disp('Mean tissue class filenames:');
disp(fn_meanTC);

%% Compare created hmri masks on NON smoothed data with original ones
clear;clc;
% both pairs are already binerazised (mask)
GM_origMask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\atlas-GM_space-MNI_mask.nii';
WM_origMask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\atlas-WM_space-MNI_mask.nii';
GM_hmriMask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\SPM12_dartel\hmri_masks\mask_sub-S001_MTsat_space-MNI_desc-mod_label-GM_probseg.nii';
WM_hmriMask_path = 'D:\Master_Thesis\Data\BIDS_AgingData\derivatives\SPM12_dartel\hmri_masks\mask_sub-S001_MTsat_space-MNI_desc-mod_label-WM_probseg.nii';

GM_origMask = spm_read_vols(spm_vol(GM_origMask_path));
WM_origMask = spm_read_vols(spm_vol(WM_origMask_path));
GM_hmriMask = spm_read_vols(spm_vol(GM_hmriMask_path));
WM_hmriMask = spm_read_vols(spm_vol(WM_hmriMask_path));

GMdiffMask_21 = GM_origMask & ~GM_hmriMask;
GMlinearIndicesDiff_21 = find(GMdiffMask_21);
[GMsubX_21, GMsubY_21, GMsubZ_21] = ind2sub(size(GMdiffMask_21), GMlinearIndicesDiff_21);

GMdiffMask_12 = GM_hmriMask & ~GM_origMask;
GMlinearIndicesDiff_12 = find(GMdiffMask_12);
[GMsubX_12, GMsubY_12, GMsubZ_12] = ind2sub(size(GMdiffMask_12), GMlinearIndicesDiff_12);

disp(['Number of differing values: ', num2str(nnz(GMdiffMask_21) + nnz(GMdiffMask_12))]);
disp(['Number of GM_hmriMask & ~GM_origMask: ', num2str(nnz(GMdiffMask_12))]);
disp(['Number of GM_origMask & ~GM_hmriMask: ', num2str(nnz(GMdiffMask_21))]);

WMdiffMask_21 = WM_origMask & ~WM_hmriMask;
WMlinearIndicesDiff_21 = find(WMdiffMask_21);
[WMsubX_21, WMsubY_21, WMsubZ_21] = ind2sub(size(WMdiffMask_21), WMlinearIndicesDiff_21);

WMdiffMask_12 = WM_hmriMask & ~WM_origMask;
WMlinearIndicesDiff_12 = find(WMdiffMask_12);
[WMsubX_12, WMsubY_12, WMsubZ_12] = ind2sub(size(WMdiffMask_12), WMlinearIndicesDiff_12);

disp(['Number of differing values: ', num2str(nnz(WMdiffMask_21) + nnz(WMdiffMask_12))]);
disp(['Number of GM_hmriMask & ~GM_origMask: ', num2str(nnz(WMdiffMask_12))]);
disp(['Number of GM_origMask & ~GM_hmriMask: ', num2str(nnz(WMdiffMask_21))]);