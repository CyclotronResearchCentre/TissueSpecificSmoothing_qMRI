%% aj_run_statistical_tools
% =========================================================================
% aj_run_statistical_tools
% =========================================================================
% Main execution script for the statistical comparison framework used to
% evaluate spatial smoothing strategies in quantitative MRI (qMRI).
%
% This pipeline centralizes the execution of several complementary analyses
% designed to compare smoothing approaches such as:
%
%   - SUSAN
%   - TSPOON
%   - TWS
%
% across multiple qMRI contrasts and tissue classes.
%
% -------------------------------------------------------------------------
% IMPLEMENTED ANALYSES
% -------------------------------------------------------------------------
%
% 1. Cluster Extraction
%    Retrieves statistical cluster information from SPM/SnPM analyses.
%
% 2. Similarity Metrics
%    Computes overlap and similarity measures between thresholded
%    statistical maps obtained with different smoothing methods.
%
% 3. Bland–Altman Analysis
%    Evaluates voxelwise agreement and systematic bias between two
%    smoothing approaches.
%
% 4. Voxelwise Log-Likelihood (vwLL)
%    Computes voxelwise Gaussian log-likelihood maps from GLM residuals
%    and determines the optimal smoothing method at each voxel.
%
% -------------------------------------------------------------------------
% CONFIGURATION FILES
% -------------------------------------------------------------------------
% The pipeline relies on dedicated configuration files:
%
%   aj_config_stat_tools :
%       General statistical analysis settings
%
%   aj_config_ba :
%       Bland–Altman configuration
%
%   aj_config_vwll :
%       Voxelwise log-likelihood configuration
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% Depending on the enabled analyses, the pipeline generates:
%
%   - Excel summary tables
%   - Similarity metric reports
%   - Bland–Altman figures
%   - Voxelwise LL maps
%   - Best-smoothing-method maps
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Antoine Jacquemin
% GIGA-CRC In Vivo Imaging
% University of Liège, Belgium
% =========================================================================

%% ========================================================================
% INITIALIZATION
% =========================================================================
clear;
close all;
clc;

%% ========================================================================
% LOAD CONFIGURATION FILES
% =========================================================================
% Load all user-defined settings and paths required for:
%   - statistical tools
%   - Bland–Altman analyses
%   - voxelwise log-likelihood analyses
% =========================================================================

cfg_stat_tools = aj_config_stat_tools();

cfg_ba = aj_config_ba();

cfg_vwll = aj_config_vwll();

%% ========================================================================
% CLUSTER EXTRACTION
% =========================================================================
% Extract cluster-level statistical information from the selected
% SPM/SnPM analysis directory.
%
% Outputs typically include:
%   - cluster size
%   - peak coordinates
%   - corrected/unorrected p-values
% =========================================================================

clusterData = aj_get_clusterData( ...
    cfg_stat_tools.base_dir1, ...
    cfg_stat_tools.combination_names, ...
    cfg_stat_tools.flag, ...
    cfg_stat_tools.excel.clusterData);

%% ========================================================================
% SIMILARITY METRICS
% =========================================================================
% Compute similarity measures between two smoothing approaches using
% thresholded statistical maps.
%
% Typical metrics may include:
%   - Dice coefficient
%   - Jaccard index
%   - voxel overlap
% =========================================================================

simMetrics = aj_get_simMetrics( ...
    cfg_stat_tools.base_dir1, ...
    cfg_stat_tools.base_dir2, ...
    cfg_stat_tools.combination_names, ...
    cfg_stat_tools.flag, ...
    cfg_stat_tools.excel.simMetrics);

%% ========================================================================
% BLAND–ALTMAN ANALYSIS
% =========================================================================
% Perform voxelwise Bland–Altman analyses between two smoothing methods
% in order to assess:
%
%   - agreement
%   - systematic bias
%   - limits of agreement
%
% across qMRI contrasts and tissue classes.
% =========================================================================

results = aj_get_BlandAltman_prepare( ...
    cfg_ba.root_dir, ...
    cfg_ba.combination_names, ...
    cfg_ba.method1_name, ...
    cfg_ba.method2_name, ...
    cfg_ba.flag);

%% ========================================================================
% VOXELWISE LOG-LIKELIHOOD ANALYSIS
% =========================================================================
% Compute voxelwise Gaussian log-likelihood maps derived from GLM
% residuals and identify the smoothing method maximizing the likelihood
% at each voxel.
%
% Outputs:
%   - LL maps for each smoothing method
%   - voxelwise best-method maps
% =========================================================================

aj_get_voxelwise_log_likelihood(cfg_vwll);
