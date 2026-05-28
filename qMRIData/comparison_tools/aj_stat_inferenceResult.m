
%% User settings
clear; clc;
addpath("C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12");

combination_names = {
    'MTsat_GM'
    'MTsat_WM'
    'PDmap_GM'
    'PDmap_WM'
    'R1map_GM'
    'R1map_WM'
    'R2starmap_GM'
    'R2starmap_WM'
};

%==========================================================================
% base_dir = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-SUSAN_GLM_rft0_age_linear';
% flag.saveExcel = 1;
% excel_path = fullfile(base_dir, 'clusterdata.xlsx');
% clusterdata = aj_get_clusterData(base_dir, combination_names, flag, excel_path);
%==========================================================================
% excel_path = fullfile('E:\Master_Thesis\Data\BIDS_AgingData\derivatives', ...
%     'qMRI_results', 'rft0', 'SimMetrics_SUSAN_vs_TWS.xlsx');
% base_dir1 = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-SUSAN_GLM_rft0_age_linear';
% base_dir2 = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS_GLM_rft0_age_linear';
% flag.saveExcel = 1;
% thrSimMetrics = aj_get_simMetrics(base_dir1, base_dir2, combination_names, flag, excel_path);
% %==========================================================================


method1_name = 'TWS'; % or TWS / TSPOON / SUSAN
method2_name = 'TSPOON'; % or TWS / TSPOON / SUSAN

flag = struct();
flag.analysis   = 'snpm';      % or 'spm_rft0' / 'spm_rft1'
flag.drawPlot   = 0;
flag.savePlot   = 1;
flag.saveExcel  = 1;
flag.paperReady = 1;

flag.thresholded = 1; % for thresholded statistical (non) parametric maps
flag.useLL = 0; % for log-likelyhood maps

if flag.useLL == 1
    flag.thresholded = 0;
    flag.analysis   = 'spm_rft1';
    root_dir = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\LL_AIC_BIC_results';
    mask_GM = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\mask_WTA_GM_mask.nii';
    mask_WM = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\mask_WTA_WM_mask.nii';
    flag.maskGM = mask_GM;
    flag.maskWM = mask_WM;
else
    root_dir = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives';
    flag.maskGM = '';
    flag.maskWM = '';
end

results = aj_BlandAltman_unified( ...
    root_dir, combination_names, ...
    method1_name, method2_name, flag);

%% Get Cluster Data for SnPM (loop + auto binarization)
% function results = aj_get_clusterData(base_dir, combination_names, flag, excel_path)
% 
% if nargin < 3
%     flag.saveExcel = 0;
% end
% 
% results = struct();
% table_rows = [];
% 
% for i = 1:length(combination_names)
%     
%     name = combination_names{i};
%     fprintf('\nProcessing: %s\n', name);
%     
%     pattern = fullfile(base_dir, name, '*FWE05_mask.nii'); % 'Tmap_thr_bin.nii'
%     d = dir(pattern);
% 
%     if ~isempty(d)
%         infile_bin = fullfile(d(1).folder, d(1).name);
%     end
%     
%     % --- Create binary map if needed ---
%     if ~exist(infile_bin, 'file')
%         
%         % --- Check thresholded map exists ---
%         infile_thr = fullfile(base_dir, name, 'Tmap_thr.nii');
%         
%         if ~exist(infile_thr, 'file')
%             warning('Thresholded file not found: %s', infile_thr);
%             continue;
%         end
%         
%         fprintf('Binary file not found → creating it...\n');
%         binarize_nifti(infile_thr, infile_bin)
%     end
%     
%     % --- Step 3: load binary image ---
%     V = spm_vol(infile_bin);
%     Y = spm_read_vols(V);
%     
%     num_significant_voxels = nnz(Y);
%     Y = double(Y);
%     
%     % --- Step 4: cluster analysis ---
%     [cluster_labels, num_clusters] = spm_bwlabel(Y, 6);
%     
%     unique_clusters = unique(cluster_labels(:));
%     unique_clusters(unique_clusters == 0) = [];
%     
%     cluster_sizes = zeros(length(unique_clusters), 1);
%     for ii = 1:length(unique_clusters)
%         cluster_sizes(ii) = sum(cluster_labels(:) == unique_clusters(ii));
%     end
%     
%     if isempty(cluster_sizes)
%         cluster_sizes = 0;
%     end
%     
%     % --- Step 5: metrics ---
%     data = struct();
%     data.nSigVox = num_significant_voxels;
%     data.nClusters = num_clusters;
%     data.clusterSizeMean = mean(cluster_sizes);
%     data.clusterSizeGeomMean = exp(mean(log(cluster_sizes + 1))) - 1;
%     data.clusterSizeMedian = median(cluster_sizes);
%     data.clusterSizeSTD = std(cluster_sizes);
%     prc = prctile(cluster_sizes, [10, 25, 50, 75, 90]);
%     
%     % Store struct
%     results.(name) = data;
%     
%     % --- Step 6: build table row ---
%     row = {name, ...
%            data.nSigVox, ...
%            data.nClusters, ...
%            data.clusterSizeMean, ...
%            data.clusterSizeGeomMean, ...
%            data.clusterSizeMedian, ...
%            data.clusterSizeSTD, ...
%            prc(1), prc(2), prc(3), prc(4), prc(5)};
%        
%     table_rows = [table_rows; row];
% end
% 
% % --- Step 7: convert to table ---
% if ~isempty(table_rows)
%     
%     T = cell2table(table_rows, 'VariableNames', { ...
%         'Combination', ...
%         'nSigVox', ...
%         'nClusters', ...
%         'MeanClusterSize', ...
%         'GeomMeanClusterSize', ...
%         'MedianClusterSize', ...
%         'STDClusterSize', ...
%         'P10', 'P25', 'P50', 'P75', 'P90'});
%     
%     % --- Step 8: save Excel ---
%     if flag.saveExcel
%         if nargin < 4 || isempty(excel_path)
%             error('You must provide excel_path when saveExcel = 1');
%         end
%         
%         writetable(T, excel_path);
%         fprintf('\nExcel file saved at:\n%s\n', excel_path);
%     end
%     
% else
%     warning('No data available to export.');
% end
% 
% end

% %% Binarize Nifti
% function aj_binarize_nifti(infile, outfile)
%     V = spm_vol(infile);
%     Y = spm_read_vols(V);
%     Ybin = double((Y ~= 0) & ~isnan(Y));
%     V.fname = outfile;
%     V.dt = [2 0];
%     spm_write_vol(V, Ybin);
%     fprintf('Binary image saved: %s\n', outfile);
% end
% 
% %% Similarity Metrics for SnPM binary maps
% function thrSimMetrics = aj_get_simMetrics(base_dir1, base_dir2, combination_names, flag, excel_path)
% 
% if nargin < 4
%     flag.saveExcel = 0;
% end
% 
% all_simMetrics = [];
% all_row_titles = {};
% 
% for i = 1:length(combination_names)
%     
%     name = combination_names{i};
%     
%     pattern1 = fullfile(base_dir1, name, '*FWE05_mask.nii'); % 'Tmap_thr_bin.nii'
%     d1 = dir(pattern1);
% 
%     if ~isempty(d1)
%         file1 = fullfile(d1(1).folder, d1(1).name);
%     end
%     
%     pattern2 = fullfile(base_dir2, name, '*FWE05_mask.nii'); % 'Tmap_thr_bin.nii'
%     d2 = dir(pattern2);
%     
%     if ~isempty(d2)
%         file2 = fullfile(d2(1).folder, d2(1).name);
%     end
%     
%     
%     fprintf('\nProcessing: %s\n and %s\n', file1, file2);
%     
%     % --- Check existence ---
%     if ~exist(file1, 'file') || ~exist(file2, 'file')
%         warning('Missing file for %s → skipping', name);
%         continue;
%     end
%     
%     % --- Load images ---
%     V1 = spm_vol(file1); Y1 = spm_read_vols(V1);
%     V2 = spm_vol(file2); Y2 = spm_read_vols(V2);
%     
%     % --- Binarize ---
%     Y1 = Y1 ~= 0;
%     Y2 = Y2 ~= 0;
%     
%     % --- Flatten ---
%     Y1 = Y1(:);
%     Y2 = Y2(:);
%     
%     % --- Remove NaNs ---
%     valid = ~isnan(Y1) & ~isnan(Y2);
%     Y1 = Y1(valid);
%     Y2 = Y2(valid);
%     
%     % --- Compute metrics ---
%     TP = sum(Y1 == 1 & Y2 == 1);
%     TN = sum(Y1 == 0 & Y2 == 0);
%     FP = sum(Y1 == 0 & Y2 == 1);
%     FN = sum(Y1 == 1 & Y2 == 0);
%     
%     % Safety (avoid division by zero)
%     if (TP + FP + FN) == 0
%         J = NaN;
%         D = NaN;
%     else
%         J = TP / (TP + FP + FN);
%         D = 2*TP / (2*TP + FP + FN);
%     end
%     
%     N = TP + TN + FP + FN;
%     po = (TP + TN) / N;
%     pe = ((TP+FP)*(TP+FN) + (FN+TN)*(FP+TN)) / N^2;
%     
%     if (1 - pe) == 0
%         K = NaN;
%     else
%         K = (po - pe) / (1 - pe);
%     end
%     
%     % --- Store ---
%     all_simMetrics = [all_simMetrics; J, D, K];
%     all_row_titles = [all_row_titles; {name}];
% end
% 
% % --- Create table only if data exist ---
% if isempty(all_simMetrics)
% 
%     warning('No similarity metrics were computed.');
%     
%     thrSimMetrics = table();
% 
% else
% 
%     % Ensure matrix format
%     all_simMetrics = reshape(all_simMetrics, [], 3);
% 
%     % --- Create table ---
%     thrSimMetrics = table( ...
%         all_row_titles, ...
%         all_simMetrics(:,1), ...
%         all_simMetrics(:,2), ...
%         all_simMetrics(:,3), ...
%         'VariableNames', {'Combination','Jaccard','Dice','CohenKappa'});
% 
%     disp(thrSimMetrics);
% 
% end
% 
% % --- Save Excel ---
% if flag.saveExcel
%     if nargin < 5 || isempty(excel_path)
%         error('You must provide excel_path when saveExcel = 1');
%     end
%     
%     writetable(thrSimMetrics, excel_path);
%     fprintf('\nExcel file saved at:\n%s\n', excel_path);
% end
% 
% end
