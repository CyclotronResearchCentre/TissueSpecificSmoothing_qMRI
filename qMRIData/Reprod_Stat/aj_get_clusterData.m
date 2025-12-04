function [GLM_data, GLM_clusterData_paths] = aj_get_clusterData(GLM_dir,pth_out,basename)
%--------------------------------------------------------------------------
% Function to extract cluster statistics and save them for GLM contrasts.
% This script processes GLM (General Linear Model) results, specifically 
% statistical maps from neuroimaging studies. It computes cluster-level 
% statistics, such as the number of significant voxels, number of clusters, 
% mean cluster size, geometric mean size, median size, standard deviation of 
% cluster sizes, and various percentiles of the cluster sizes. These metrics 
% are computed for each GLM contrast (e.g., for different statistical tests) 
% based on thresholding the statistical maps at a Family-Wise Error (FWE) 
% corrected p-value threshold.

% The function saves the resulting cluster statistics in a table format as 
% `.mat` files, and returns the statistics and paths to the saved files.

% INPUTS
% GLM_dir:      The directory where the GLM statistical maps (spmT_*.nii, 
%               spmF_*.nii) and SPM.mat file are located.
% pth_out:      The output directory where the cluster data tables will be 
%               saved.
% basename:     A string to use as a prefix for the output files.

% OUTPUTS
% GLM_data:     A cell array of structures containing cluster statistics 
%               for each GLM contrast.
% GLM_clusterData_paths: A cell array of paths to the saved cluster data 
%               `.mat` files.

% PROCESS
% The script performs the following steps:
% 1. It reads the GLM directory to find all spmT_*.nii and spmF_*.nii files, 
%    which correspond to statistical test results for each contrast.
% 2. The SPM.mat file is loaded to obtain the contrasts setup in the GLM model.
% 3. For each contrast, the script generates an xSPM structure and computes 
%    the threshold using the spm_uc function, which calculates the FWE-corrected 
%    threshold based on the degrees of freedom and the search space.
% 4. The statistical map is thresholded, and the binary mask of significant 
%    voxels is computed.
% 5. The script identifies clusters in the binary mask using spm_bwlabel 
%    and calculates cluster sizes.
% 6. Various statistical metrics (e.g., mean size, median size, percentiles) 
%    are calculated for the clusters and saved in a table.
% 7. The data is stored in a `.mat` file, and the path to this file is saved.
%
% REFERENCE
% L. Thurfjell & al. (1992) https://doi.org/10.1016/1049-9652(92)90083-A
% SPM documentation https://www.fil.ion.ucl.ac.uk/spm/doc/
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Dealing with inputs
if nargin<2
    error('aj_get_xSPM: Not enough inputs.');
end

if ~exist(GLM_dir, 'dir')
    error('aj_get_xSPM: the GLM directory does not exist.');
end

if ~exist(pth_out, 'dir')
    mkdir(pth_out);
end

% Store the initial directory
initial_dir = pwd;

% Find all spmT_000x.nii and spmF_000x.nii files
% sorted as : spmT_0001, spmT_0002, spmF_0003
file_list = dir(fullfile(GLM_dir, 'spmT_*.nii'));
file_list = [file_list; dir(fullfile(GLM_dir, 'spmF_*.nii'))];
nContrasts = length(file_list);

% Load SPM.mat to set up the contrast
load(fullfile(GLM_dir, 'SPM.mat'), 'SPM');

%% Do the job
GLM_data = cell(nContrasts,1);
GLM_clusterData_paths = cell(nContrasts,1);
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
    binary_mask = Y > u_fwe;

    % Return to the initial directory at the end of the script
    cd(initial_dir);

    % Count significant voxels
    num_significant_voxels = nnz(binary_mask);

    % Label clusters
    binary_mask = double(binary_mask); % BW must be numeric, real, full and double
    [cluster_labels, num_clusters] = spm_bwlabel(binary_mask, 6); % connectivity criterion = surface (6)

    % Calculate cluster sizes
    unique_clusters = unique(cluster_labels(:)); % Unique cluster labels
    unique_clusters(unique_clusters == 0) = []; % Remove background (label 0)
    cluster_sizes = zeros(length(unique_clusters), 1);
    for ii = 1:length(unique_clusters)
        cluster_sizes(ii) = sum(cluster_labels(:) == unique_clusters(ii));
    end

    % Calculate statistical metrics
    data = struct();
    data.nSigVox = num_significant_voxels;
    data.nClusters = num_clusters;
    data.mean_size = mean(cluster_sizes);
    data.geom_mean = exp(mean(log(cluster_sizes + 1))) - 1;
    data.median_size = median(cluster_sizes);
    data.std_size = std(cluster_sizes);
    data.percentiles = prctile(cluster_sizes, [10, 25, 50, 75, 90]);
    
    GLM_data{i} = data;
    
    % Convert the data into a table
    clusterDataTable = table( ...
        num_significant_voxels, ...
        num_clusters, ...
        data.mean_size, ...
        data.geom_mean, ...
        data.median_size, ...
        data.std_size, ...
        data.percentiles(1), data.percentiles(2), data.percentiles(3), data.percentiles(4), data.percentiles(5), ...
        'VariableNames', {'NumSigVoxels', 'NumClusters', 'MeanSize', 'GeomMean', 'MedianSize', 'StdSize', 'P10', 'P25', 'P50', 'P75', 'P90'} ...
    );

    % Save the table into a .mat file
    basename_i = sprintf('%s_%s',basename,strrep(xSPM.title,' ',''));
    pth_out_i = fullfile(pth_out,basename_i);
    save(pth_out_i, 'clusterDataTable');
    
    % Store restult paths
    if isempty(GLM_clusterData_paths{i})
        GLM_clusterData_paths{i} = pth_out_i;
    else
        GLM_clusterData_paths{i} = char(GLM_clusterData_paths{i}, pth_out_i);
    end
    
    if isequal(strrep(xSPM.title,' ',''), 'Increasewithage')
        % Log-transform cluster sizes
        log_cluster_sizes = log10(cluster_sizes);
        figure;
        histogram(log_cluster_sizes, 'BinWidth', 0.1);
        xlabel('Log10(Cluster Size)');
        ylabel('Frequency');
        title('Log-Scaled Distribution of Cluster Sizes');
        grid on;
        
        % Compute cumulative percentages
        sorted_cluster_sizes = sort(cluster_sizes);
        cum_percentages = (1:length(sorted_cluster_sizes)) / length(sorted_cluster_sizes) * 100;
        figure;
        plot(sorted_cluster_sizes, cum_percentages, '-o', 'LineWidth', 1.5);
        xlabel('Cluster Size (Number of Voxels)');
        ylabel('Cumulative Percentage (%)');
        title('Cumulative Distribution of Cluster Sizes');
        grid on;

        % Bar Chart of Binned Sizes
        edges = [0, 1, 10, 50, 100, 500];
        max_cluster_size = max(cluster_sizes);
        if max_cluster_size > edges(end)
            edges = [edges, max_cluster_size];
        end
%         bin_labels = {'0-1', '2-10', '11-50', '51-100', '101-500', sprintf('%d-%d', edges(end-1)+1, edges(end))};
        % Ensure bin_labels are unique
        if max_cluster_size > edges(end-1)
            bin_labels = {'0-1', '2-10', '11-50', '51-100', '101-500', sprintf('%d-%d', edges(end-1)+1, edges(end))};
        else
            bin_labels = {'0-1', '2-10', '11-50', '51-100', '101-500'};
        end
        bin_counts = histcounts(cluster_sizes, edges);
        figure;
        bar(1:length(bin_counts), bin_counts);
        xlabel('Cluster Size Ranges');
        ylabel('Number of Clusters');
        title('Distribution of Cluster Sizes by Range');
        grid on;
%         edges = [0, 1, 10, 50, 100, 500, max(cluster_sizes)];
%         bin_labels = {'0-1', '2-10', '11-50', '51-100', '101-500', '>500'};
%         bin_counts = histcounts(cluster_sizes, edges);
%         figure;
%         bar(categorical(bin_labels, bin_labels), bin_counts);
%         xlabel('Cluster Size Ranges');
%         ylabel('Number of Clusters');
%         title('Distribution of Cluster Sizes by Range');
%         grid on;
    end

end

fprintf('Cluster data are saved in: %s \n',pth_out);

%% Print results
% fprintf('Cluster size statistics:\n');
% fprintf('Number of significant voxels: %d voxels\n', num_significant_voxels);
% fprintf('Number of clusters: %d clusters\n', num_clusters);
% fprintf('Mean: %.2f voxels\n', mean_size);
% fprintf('Geometric Mean: %.2f voxels\n', geom_mean);
% fprintf('Median: %.2f voxels\n', median_size);
% fprintf('Standard Deviation: %.2f voxels\n', std_size);
% 
% fprintf('10th Percentile: %.2f voxels\n', percentiles(1));
% fprintf('25th Percentile: %.2f voxels\n', percentiles(2));
% fprintf('50th Percentile (Median): %.2f voxels\n', percentiles(3));
% fprintf('75th Percentile: %.2f voxels\n', percentiles(4));
% fprintf('90th Percentile: %.2f voxels\n', percentiles(5));

%% Plot clusters distribtution
% % Plot the histogram for filtered sizes
% threshold_size = 10;
% filtered_cluster_sizes = cluster_sizes(cluster_sizes >= threshold_size);
% figure;
% histogram(filtered_cluster_sizes, 'BinWidth', 25);
% xlabel('Cluster Size (Number of Voxels)');
% ylabel('Frequency');
% title(sprintf('Histogram of Cluster Sizes (<= %d Voxels)', threshold_size));
% grid on;
% 
% % Log-transform cluster sizes
% log_cluster_sizes = log10(cluster_sizes);
% figure;
% histogram(log_cluster_sizes, 'BinWidth', 0.1);
% xlabel('Log10(Cluster Size)');
% ylabel('Frequency');
% title('Log-Scaled Distribution of Cluster Sizes');
% grid on;
% 
% % Boxplot for Outlier Identification
% figure;
% boxplot(cluster_sizes, 'Orientation', 'horizontal');
% xlabel('Cluster Size (Number of Voxels)');
% title('Boxplot of Cluster Sizes');
% grid on;
% 
% % Compute cumulative percentages
% sorted_cluster_sizes = sort(cluster_sizes);
% cum_percentages = (1:length(sorted_cluster_sizes)) / length(sorted_cluster_sizes) * 100;
% figure;
% plot(sorted_cluster_sizes, cum_percentages, '-o', 'LineWidth', 1.5);
% xlabel('Cluster Size (Number of Voxels)');
% ylabel('Cumulative Percentage (%)');
% title('Cumulative Distribution of Cluster Sizes');
% grid on;
% 
% % Bar Chart of Binned Sizes
% edges = [0, 1, 10, 50, 100, 500, max(cluster_sizes)];
% bin_labels = {'0-1', '2-10', '11-50', '51-100', '101-500', '>500'};
% bin_counts = histcounts(cluster_sizes, edges);
% figure;
% bar(categorical(bin_labels, bin_labels), bin_counts);
% xlabel('Cluster Size Ranges');
% ylabel('Number of Clusters');
% title('Distribution of Cluster Sizes by Range');
% grid on;

end