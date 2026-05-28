%% Similarity Metrics for Thresholded SnPM/SPM Binary Maps
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% This function computes spatial similarity metrics between two sets of
% binary thresholded statistical maps (e.g., SnPM or SPM FWE-corrected
% masks).
%
% For each qMRI combination, the function automatically searches for
% binary statistical masks matching:
%
%     *FWE05_mask.nii
%
% and computes the following similarity measures:
%
%   - Jaccard Index
%   - Dice Coefficient
%   - Cohen's Kappa
%
% The comparison is performed voxelwise after:
%   - binarization of the images
%   - flattening into vectors
%   - removal of NaN voxels
%
% The function can optionally export the results into an Excel file.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% base_dir1 : char
%     Root directory of the first method/statistical pipeline.
%
% base_dir2 : char
%     Root directory of the second method/statistical pipeline.
%
% combination_names : cell array of char
%     Cell array containing qMRI combinations/folders, e.g.:
%
%         {'MTsat_GM', 'MTsat_WM', 'R1map_GM', ...}
%
% flag : struct
%     Structure containing optional processing flags.
%
%     flag.saveExcel : logical (default = 0)
%         If true, saves the similarity metrics into an Excel file.
%
% excel_path : char (optional)
%     Full path to the Excel output file.
%     Required only if flag.saveExcel = 1.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% thrSimMetrics : table
%     Table containing similarity metrics for each qMRI combination.
%
%     Columns:
%         - Combination
%         - Jaccard
%         - Dice
%         - CohenKappa
%
% -------------------------------------------------------------------------
% METRICS
% -------------------------------------------------------------------------
% Jaccard Index:
%
%                   TP
%     J = ---------------------
%         TP + FP + FN
%
% Dice Coefficient:
%
%                  2TP
%     D = ---------------------
%         2TP + FP + FN
%
% Cohen's Kappa:
%
%              po - pe
%     K = ----------------
%             1 - pe
%
% where:
%
%     TP = True Positives
%     TN = True Negatives
%     FP = False Positives
%     FN = False Negatives
%
% -------------------------------------------------------------------------
% REQUIREMENTS
% -------------------------------------------------------------------------
% - SPM12 must be added to the MATLAB path.
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
% excel_path = fullfile(outdir, 'SimilarityMetrics.xlsx');
%
% thrSimMetrics = aj_get_simMetrics( ...
%     base_dir1, ...
%     base_dir2, ...
%     combination_names, ...
%     flag, ...
%     excel_path);
%
% -------------------------------------------------------------------------

function thrSimMetrics = aj_get_simMetrics( ...
    base_dir1, ...
    base_dir2, ...
    combination_names, ...
    flag, ...
    excel_path)

%--------------------------------------------------------------------------
% Default flags
%--------------------------------------------------------------------------
if nargin < 4
    flag.saveExcel = 0;
end

% Initialize outputs
all_simMetrics = [];
all_row_titles = {};

%--------------------------------------------------------------------------
% Loop across qMRI combinations
%--------------------------------------------------------------------------
for i = 1:length(combination_names)
    
    name = combination_names{i};
    
    %----------------------------------------------------------------------
    % Search binary masks for method 1
    %----------------------------------------------------------------------
    pattern1 = fullfile(base_dir1, name, '*FWE05_mask.nii');
    d1 = dir(pattern1);

    if ~isempty(d1)
        file1 = fullfile(d1(1).folder, d1(1).name);
    end
    
    %----------------------------------------------------------------------
    % Search binary masks for method 2
    %----------------------------------------------------------------------
    pattern2 = fullfile(base_dir2, name, '*FWE05_mask.nii');
    d2 = dir(pattern2);
    
    if ~isempty(d2)
        file2 = fullfile(d2(1).folder, d2(1).name);
    end
    
    fprintf('\nProcessing:\n%s\n%s\n', file1, file2);
    
    %----------------------------------------------------------------------
    % Check file existence
    %----------------------------------------------------------------------
    if ~exist(file1, 'file') || ~exist(file2, 'file')
        warning('Missing file for %s -> skipping', name);
        continue;
    end
    
    %----------------------------------------------------------------------
    % Load images
    %----------------------------------------------------------------------
    V1 = spm_vol(file1);
    Y1 = spm_read_vols(V1);
    
    V2 = spm_vol(file2);
    Y2 = spm_read_vols(V2);
    
    %----------------------------------------------------------------------
    % Binarization
    %----------------------------------------------------------------------
    Y1 = Y1 ~= 0;
    Y2 = Y2 ~= 0;
    
    %----------------------------------------------------------------------
    % Flatten images into vectors
    %----------------------------------------------------------------------
    Y1 = Y1(:);
    Y2 = Y2(:);
    
    %----------------------------------------------------------------------
    % Remove NaN voxels
    %----------------------------------------------------------------------
    valid = ~isnan(Y1) & ~isnan(Y2);
    
    Y1 = Y1(valid);
    Y2 = Y2(valid);
    
    %----------------------------------------------------------------------
    % Compute confusion matrix terms
    %----------------------------------------------------------------------
    TP = sum(Y1 == 1 & Y2 == 1);
    TN = sum(Y1 == 0 & Y2 == 0);
    FP = sum(Y1 == 0 & Y2 == 1);
    FN = sum(Y1 == 1 & Y2 == 0);
    
    %----------------------------------------------------------------------
    % Compute Jaccard and Dice metrics
    %----------------------------------------------------------------------
    % Safety check to avoid division by zero
    if (TP + FP + FN) == 0
        
        J = NaN;
        D = NaN;
        
    else
        
        J = TP / (TP + FP + FN);
        
        D = 2*TP / (2*TP + FP + FN);
    end
    
    %----------------------------------------------------------------------
    % Compute Cohen's Kappa
    %----------------------------------------------------------------------
    N = TP + TN + FP + FN;
    
    po = (TP + TN) / N;
    
    pe = ((TP+FP)*(TP+FN) + (FN+TN)*(FP+TN)) / N^2;
    
    % Safety check
    if (1 - pe) == 0
        
        K = NaN;
        
    else
        
        K = (po - pe) / (1 - pe);
    end
    
    %----------------------------------------------------------------------
    % Store results
    %----------------------------------------------------------------------
    all_simMetrics = [all_simMetrics; J, D, K];
    
    all_row_titles = [all_row_titles; {name}];
end

%--------------------------------------------------------------------------
% Create output table
%--------------------------------------------------------------------------
if isempty(all_simMetrics)

    warning('No similarity metrics were computed.');
    
    thrSimMetrics = table();

else

    % Ensure matrix format
    all_simMetrics = reshape(all_simMetrics, [], 3);

    %----------------------------------------------------------------------
    % Create table
    %----------------------------------------------------------------------
    thrSimMetrics = table( ...
        all_row_titles, ...
        all_simMetrics(:,1), ...
        all_simMetrics(:,2), ...
        all_simMetrics(:,3), ...
        'VariableNames', { ...
            'Combination', ...
            'Jaccard', ...
            'Dice', ...
            'CohenKappa'});

    disp(thrSimMetrics);

end

%--------------------------------------------------------------------------
% Save Excel file
%--------------------------------------------------------------------------
if flag.saveExcel
    
    if nargin < 5 || isempty(excel_path)
        error('You must provide excel_path when saveExcel = 1');
    end
    
    writetable(thrSimMetrics, excel_path);
    
    fprintf('\nExcel file saved at:\n%s\n', excel_path);
end

end