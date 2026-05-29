%% Get Cluster Data from Thresholded SnPM/SPM Maps
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% This function extracts cluster-level descriptive statistics from
% thresholded statistical maps (e.g., SnPM or SPM voxelwise FWE-corrected
% masks). The function automatically searches for an existing binary mask
% (*FWE05_mask.nii). If no binary image is found, a binary mask is created
% from the thresholded statistical image (Tmap_thr.nii).
%
% For each qMRI combination, the function computes:
%   - Number of significant voxels
%   - Number of clusters
%   - Mean cluster size
%   - Geometric mean cluster size
%   - Median cluster size
%   - Standard deviation of cluster sizes
%   - Cluster size percentiles (10th, 25th, 50th, 75th, 90th)
%
% Cluster identification is performed using SPM's spm_bwlabel function
% with 6-connectivity (surface connectivity criterion).
%
% Optionally, results can be exported to an Excel spreadsheet.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% base_dir : char
%     Root directory containing statistical result folders.
%
% combination_names : cell array of char
%     Cell array containing the names of qMRI combinations/folders, e.g.:
%         {'MTsat_GM', 'MTsat_WM', 'R1map_GM', ...}
%
% flag : struct
%     Structure containing optional processing flags:
%
%     flag.saveExcel : logical (default = 0)
%         If true, saves the extracted metrics into an Excel file.
%
% excel_path : char (optional)
%     Full path to the Excel output file.
%     Required only if flag.saveExcel = 1.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% results : struct
%     Structure containing cluster statistics for each combination.
%
%     Example:
%         results.MTsat_GM.nSigVox
%         results.MTsat_GM.nClusters
%         results.MTsat_GM.clusterSizeMean
%         ...
%
% -------------------------------------------------------------------------
% REQUIREMENTS
% -------------------------------------------------------------------------
% - SPM12 must be added to the MATLAB path.
% - Helper function "aj_binarize_nifti.m" must be available.
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Antoine Jacquemin
%
% -------------------------------------------------------------------------
% EXAMPLE
% -------------------------------------------------------------------------
% flag.saveExcel = 1;
%
% excel_path = fullfile(outdir, 'clusterData.xlsx');
%
% results = aj_get_clusterData( ...
%     base_dir, ...
%     combination_names, ...
%     flag, ...
%     excel_path);
%
% -------------------------------------------------------------------------

function results = aj_get_clusterData(base_dir, combination_names, flag, excel_path)

%--------------------------------------------------------------------------
% Default flags
%--------------------------------------------------------------------------
if nargin < 3
    flag.saveExcel = 0;
end

% Initialize outputs
results = struct();
table_rows = [];

%--------------------------------------------------------------------------
% Loop across qMRI combinations
%--------------------------------------------------------------------------
for i = 1:length(combination_names)
    
    name = combination_names{i};
    fprintf('\nProcessing: %s\n', name);
    
    %----------------------------------------------------------------------
    % Search for an existing binary thresholded mask
    %----------------------------------------------------------------------
    pattern = fullfile(base_dir, name, '*FWE05_mask.nii');
    d = dir(pattern);

    if ~isempty(d)
        infile_bin = fullfile(d(1).folder, d(1).name);
    end
    
    %----------------------------------------------------------------------
    % Create binary image if it does not exist
    %----------------------------------------------------------------------
    if ~exist(infile_bin, 'file')
        
        % Define thresholded statistical map
        infile_thr = fullfile(base_dir, name, 'Tmap_thr.nii');
        
        % Check existence
        if ~exist(infile_thr, 'file')
            warning('Thresholded file not found: %s', infile_thr);
            continue;
        end
        
        fprintf('Binary file not found -> creating it...\n');
        aj_binarize_nifti(infile_thr, infile_bin)
    end
    
    %----------------------------------------------------------------------
    % Load binary image
    %----------------------------------------------------------------------
    V = spm_vol(infile_bin);
    Y = spm_read_vols(V);
    
    % Count significant voxels
    num_significant_voxels = nnz(Y);
    
    % Convert to double precision for spm_bwlabel
    Y = double(Y);
    
    %----------------------------------------------------------------------
    % Cluster analysis
    %----------------------------------------------------------------------
    % 6-connectivity corresponds to surface connectivity
    [cluster_labels, num_clusters] = spm_bwlabel(Y, 6);
    
    % Retrieve cluster labels
    unique_clusters = unique(cluster_labels(:));
    unique_clusters(unique_clusters == 0) = [];
    
    % Compute cluster sizes
    cluster_sizes = zeros(length(unique_clusters), 1);
    
    for ii = 1:length(unique_clusters)
        cluster_sizes(ii) = ...
            sum(cluster_labels(:) == unique_clusters(ii));
    end
    
    % Handle empty cases
    if isempty(cluster_sizes)
        cluster_sizes = 0;
    end
    
    %----------------------------------------------------------------------
    % Compute cluster metrics
    %----------------------------------------------------------------------
    data = struct();
    
    data.nSigVox            = num_significant_voxels;
    data.nClusters          = num_clusters;
    data.clusterSizeMean    = mean(cluster_sizes);
    data.clusterSizeGeomMean = ...
        exp(mean(log(cluster_sizes + 1))) - 1;
    data.clusterSizeMedian  = median(cluster_sizes);
    data.clusterSizeSTD     = std(cluster_sizes);
    
    % Cluster size percentiles
    prc = prctile(cluster_sizes, [10, 25, 50, 75, 90]);
    
    % Store metrics into output structure
    results.(name) = data;
    
    %----------------------------------------------------------------------
    % Build table row for Excel export
    %----------------------------------------------------------------------
    row = { ...
        name, ...
        data.nSigVox, ...
        data.nClusters, ...
        data.clusterSizeMean, ...
        data.clusterSizeGeomMean, ...
        data.clusterSizeMedian, ...
        data.clusterSizeSTD, ...
        prc(1), prc(2), prc(3), prc(4), prc(5)};
       
    table_rows = [table_rows; row];
end

%--------------------------------------------------------------------------
% Convert results to table
%--------------------------------------------------------------------------
if ~isempty(table_rows)
    
    T = cell2table(table_rows, ...
        'VariableNames', { ...
        'Combination', ...
        'nSigVox', ...
        'nClusters', ...
        'MeanClusterSize', ...
        'GeomMeanClusterSize', ...
        'MedianClusterSize', ...
        'STDClusterSize', ...
        'P10', 'P25', 'P50', 'P75', 'P90'});
    
    %----------------------------------------------------------------------
    % Save Excel file if requested
    %----------------------------------------------------------------------
    if flag.saveExcel
        
        if nargin < 4 || isempty(excel_path)
            error('You must provide excel_path when saveExcel = 1');
        end
        
        writetable(T, excel_path);
        
        fprintf('\nExcel file saved at:\n%s\n', excel_path);
    end
    
else
    warning('No data available to export.');
end

end

%% Binarize NIfTI Image
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% This function converts a NIfTI image into a binary mask.
%
% All voxels satisfying:
%       voxel ~= 0 AND voxel is not NaN
% are assigned a value of 1.
%
% All remaining voxels are assigned a value of 0.
%
% The resulting binary image is saved as a new NIfTI file using SPM.
%
% This helper function is particularly useful for:
%   - Thresholded statistical maps
%   - Significant voxel masks
%   - SnPM/SPM FWE-corrected images
%   - Cluster-based analyses
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% infile : char
%     Full path to the input NIfTI image.
%
% outfile : char
%     Full path to the output binary NIfTI image.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% A binary NIfTI image is written to disk:
%
%     1 = significant/non-zero valid voxel
%     0 = background, zero or NaN voxel
%
% The output datatype is uint8-compatible integer format (SPM datatype 2).
%
% -------------------------------------------------------------------------
% REQUIREMENTS
% -------------------------------------------------------------------------
% - SPM12 must be added to the MATLAB path.
% -------------------------------------------------------------------------


function aj_binarize_nifti(infile, outfile)

    %----------------------------------------------------------------------
    % Load input NIfTI image
    %----------------------------------------------------------------------
    V = spm_vol(infile);
    Y = spm_read_vols(V);

    %----------------------------------------------------------------------
    % Binarization
    %----------------------------------------------------------------------
    % Keep only non-zero and non-NaN voxels
    Ybin = double((Y ~= 0) & ~isnan(Y));

    %----------------------------------------------------------------------
    % Define output image properties
    %----------------------------------------------------------------------
    V.fname = outfile;

    % SPM datatype = 2 → unsigned 8-bit integer
    V.dt = [2 0];

    %----------------------------------------------------------------------
    % Save binary image
    %----------------------------------------------------------------------
    spm_write_vol(V, Ybin);

    fprintf('Binary image saved: %s\n', outfile);

end