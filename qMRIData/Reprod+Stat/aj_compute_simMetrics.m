function [J, D, K] = aj_compute_simMetrics(Y1, Y2)
%--------------------------------------------------------------------------
% Function to compute similarity metrics between two binary 3D matrices.
% This function computes three common similarity metrics used to assess the 
% overlap between two binary masks, typically in neuroimaging or other 
% biomedical contexts. The function calculates:
% 1. The Jaccard Index (J)
% 2. The Dice-Sørensen Coefficient (D)
% 3. Cohen's Kappa (K)
%
% These metrics compare the agreement between two binary 3D matrices, 
% which could represent regions of interest (ROIs), significant voxels, 
% or other binary data.
%
% The function assumes that the inputs `Y1` and `Y2` are binary 3D matrices 
% (e.g., derived from thresholded MPMs). 
% The output values provide measures of the overlap and agreement between 
% these matrices.

% INPUTS
% - Y1: A binary 3D matrix representing the first data set (e.g., a mask of 
%   significant voxels from the first map).
% - Y2: A binary 3D matrix representing the second data set (e.g., a mask 
%   of significant voxels from the second map).

% OUTPUTS
% - J: Jaccard Index, a measure of the similarity between the two binary masks 
%   defined as the ratio of the intersection to the union of the two sets.
% - D: Dice-Sørensen Coefficient, a measure of overlap between the two masks.
% - K: Cohen's Kappa, a measure of agreement between the two binary masks 
%   that accounts for chance agreement.

% PROCESS
% 1. The function first checks if the two input matrices `Y1` and `Y2` have the 
%    same size. If not, an error is raised.
% 2. It then thresholds the matrices to create binary masks (i.e., any non-zero 
%    values are considered significant).
% 3. The Jaccard Index is computed by dividing the number of overlapping significant 
%    voxels (intersection) by the total number of unique significant voxels (union) 
%    between the two masks.
% 4. The Dice-Sørensen Coefficient is calculated, which is a weighted measure 
%    of overlap that gives more importance to smaller values.
% 5. Cohen's Kappa is computed to assess the agreement between the two masks, 
%    taking into account the expected agreement by chance. The formula for Cohen's 
%    Kappa accounts for the probability of agreement by chance based on the 
%    proportions of significant voxels in each mask.
% 6. If any of the masks is empty (i.e., no significant voxels), warnings are 
%    generated, and the metrics are set to `NaN` (not a number), indicating that 
%    the metrics cannot be computed for empty masks.

% REFERENCES
% Jaccard Index: https://en.wikipedia.org/wiki/Jaccard_index
% Dice-Sørensen coefficient: https://en.wikipedia.org/wiki/Dice-Sørensen_coefficient
% Cohen's Kappa: https://en.wikipedia.org/wiki/Cohen%27s_kappa
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Ensure they are of the same dimensions
if ~isequal(size(Y1), size(Y2))
    error('aj_compute_similarity_metrics: NIfTI files must have the same dimensions.');
end

%% Threshold to create binary masks
binary_mask1 = Y1 ~= 0; % Significant voxels in file 1
binary_mask2 = Y2 ~= 0; % Significant voxels in file 2

%% Compute the Jaccard Index
intersection = nnz(binary_mask1 & binary_mask2); % Number of overlapping significant voxels
union = nnz(binary_mask1 | binary_mask2);        % Total number of unique significant voxels

if union == 0
    warning('aj_compute_similarity_metrics: Both masks are empty. Metrics are undefined, returning NaN.');
    J = NaN; % Jaccard index undefined for empty masks
    D = NaN; % Dice coefficient undefined for empty masks
    K = NaN; % Cohen's kappa undefined for empty masks
    return;
else
    J = intersection / union; % Jaccard index
end

%% Compute the Dice-Sørensen Coefficient
size_A = nnz(binary_mask1); % Total significant voxels in mask 1
size_B = nnz(binary_mask2); % Total significant voxels in mask 2
D = 2 * intersection / (size_A + size_B); % Dice coefficient

%% Compute Cohen's Kappa
total_voxels = numel(binary_mask1); % Total number of voxels
observed_agreement = nnz(binary_mask1 == binary_mask2) / total_voxels; % Fraction of matching voxels

% Expected agreement
pA = size_A / total_voxels; % Fraction of significant voxels in mask 1
pB = size_B / total_voxels; % Fraction of significant voxels in mask 2
expected_agreement = (pA * pB) + ((1 - pA) * (1 - pB));

if expected_agreement == 1
    warning('aj_compute_similarity_metrics: Cohen''s Kappa is undefined due to perfect expected agreement.');
    K = NaN; % Cohen's Kappa undefined if expected agreement is 1
else
    K = (observed_agreement - expected_agreement) / (1 - expected_agreement); % Cohen's Kappa
end

% fprintf('Similarity metrics between:\n');
% fprintf('  Jaccard index (J): %.4f\n', J);
% fprintf('  Dice-Sørensen coefficient (D): %.4f\n', D);
% fprintf('  Cohen''s Kappa (K): %.4f\n', K);
end
