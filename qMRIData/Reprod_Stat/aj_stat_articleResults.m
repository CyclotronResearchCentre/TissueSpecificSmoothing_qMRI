% Main script to analyze the statistical parametric maps obtained during
% the reproduction of the article results (Callaghan et al (2014)). Some
% statistical approaches are implemented as Bland-Altman Plot, Theshold
% Scatter Plot, Similarity Metrics (Jaccard, Dice and Cohen's kappa),
% Cluster Level Comparison and Effective Smoothing.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Starter
close all; clear; clc;

% Choose between 'Perso PC', 'Desktop'
flags.users = 0;
paths = aj_default_GLM(flags); % Get user paths, set up spm12 env

% smoo_approachs = {'TWsmoot', 'TWS', 'TSPOON', 'SUSAN'};
smoo_approachs = {'TWS', 'TSPOON', 'SUSAN'};
nsmoo = length(smoo_approachs);

qmetrics = {'MTsat', 'PDmap', 'R1map', 'R2starmap'};
TCs = {'GM', 'WM'};
combinations = {};
for ii = 1:length(qmetrics)
    for j = 1:length(TCs)
        combinations{end+1} = [qmetrics{ii}, '_', TCs{j}];
    end
end

nCombi = length(combinations);

%% Getting information from xSPM
resutls_info = cell(nsmoo,1);
for i = 1:nsmoo
    params = cell(nCombi,1);
    for ii = 1:nCombi
        GLM_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM_rft1',smoo_approachs{i}),combinations{ii});
        
        % Get xSPM all the contrasts contained in this GLM_dir
        GLM_xSPM = aj_get_xSPM(GLM_dir);

        % Initialize new struct to store specific parameters
        params{ii} = struct();

        % Extract specific parameters from each xSPM
        for iii = 1:length(GLM_xSPM)
            params{ii}.swd{iii} = GLM_xSPM{iii}.swd;
            params{ii}.title{iii} = GLM_xSPM{iii}.title;
            params{ii}.u{iii} = GLM_xSPM{iii}.u;
            params{ii}.FWHM{iii} = GLM_xSPM{iii}.FWHM;
            params{ii}.nClusters{iii} = str2double(regexp(GLM_xSPM{iii}.STATstr, '\d+', 'match'));
        end
    end
    resutls_info{i} = params;
end

%% Get all thresholds T-values
% after statistical test using p<0.05 FWE corrected level
contrasts = {'T1', 'T2', 'F3'};

all_GLM_xSPM = cell(nsmoo,1);
for i = 1:nsmoo
    GLM_xSPM = cell(nCombi,1);
    for ii = 1:nCombi
        GLM_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM_rft0',smoo_approachs{i}),combinations{ii});
        GLM_xSPM{ii} = aj_get_xSPM(GLM_dir);
    end
    all_GLM_xSPM{i} = GLM_xSPM;
end

% Preallocate or initialize arrays to store the values
nIterations = nCombi*length(contrasts);
TWS_u = zeros(1, nIterations);
TSPOON_u = zeros(1, nIterations);
SUSAN_u = zeros(1, nIterations);

% Index to track where to insert values
index = 1;
for i = 1:nCombi
    % Retrieve the 3-cells array for the specified combination
    TWS_xSPM = all_GLM_xSPM{1}{i};
    TSPOON_xSPM = all_GLM_xSPM{2}{i};
    SUSAN_xSPM = all_GLM_xSPM{3}{i};
    
    for ii = 1:length(SUSAN_xSPM)
        % Assign the values directly to the preallocated arrays
        TWS_u(index) = TWS_xSPM{ii}.u;
        TSPOON_u(index) = TSPOON_xSPM{ii}.u;
        SUSAN_u(index) = SUSAN_xSPM{ii}.u;
        index = index + 1; % Increment the index
    end
end

%% Threshold T-values Scatter Plot
% Before run this section: have to get all T-values
% Can choose the 2 smoothing_method
% after statistical test using p<0.05 FWE corrected level

meth1 = 2;
meth2 = 3;

data_info = struct();
data_info.metric1 = smoo_approachs{meth1};
data_info.metric2 = smoo_approachs{meth2};
pth_out_dir = fullfile(paths.ds_dir,'derivatives','qMRI_results','rft1','thrScatterPlot');

% Divide indices into three groups
indices_T1 = 1:3:length(SUSAN_u); % Indices of T1 contrast
indices_T2 = 2:3:length(SUSAN_u); % Indices of T2 contrast
indices_F3 = 3:3:length(SUSAN_u); % Indices of F3 contrast

% Loop through each set of contrasts (T1, T2, F3)
for i = 1:length(contrasts)
    currentContrast = contrasts{i};
    
    % Get the corresponding data for the current index (T1, T2, F3)
    eval([data_info.metric1 '_u_' currentContrast ' = ' data_info.metric1 '_u(indices_' currentContrast ');']);
    eval([data_info.metric2 '_u_' currentContrast ' = ' data_info.metric2 '_u(indices_' currentContrast ');']);
    
    data_info.contrast = currentContrast;
    
    % Call the plotting function
    pth_out = fullfile(pth_out_dir,sprintf('%s-%s_%s.png',data_info.metric1, data_info.metric2, currentContrast));
    aj_scatterPlot(eval([data_info.metric1 '_u_' currentContrast]), eval([data_info.metric2 '_u_' currentContrast]),data_info,pth_out);
end

%% Cluster Level Comparison: Finding the numbers of significant voxels and clusters
pth_out = fullfile(paths.ds_dir,'derivatives','qMRI_results','rft1','ClusterData');

all_GLM_data = cell(nsmoo,1);
Smoot_GLM_clusterData = cell(nsmoo,1);
for i = 1:nsmoo
    GLM_data = cell(nCombi,1);
    GLM_clusterData = cell(nCombi,1);
    for ii = 1:nCombi
        GLM_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM_rft1',smoo_approachs{i}),combinations{ii});
        basename = sprintf('%s_%s',smoo_approachs{i},combinations{ii});
        [GLM_data{ii}, GLM_clusterData{ii}] = aj_get_clusterData(GLM_dir,pth_out,basename);
    end
    all_GLM_data{i} = GLM_data;
    Smoot_GLM_clusterData{i} = GLM_clusterData;
end

% Save the cell array into the .mat file
output_file = fullfile(pth_out,'all_GLM_data.mat');
save(output_file, 'all_GLM_data');

%% Optionnal: Merging all output files into ONE
% Specify the folder containing the .mat files
folderPath = pth_out;

% Get the list of .mat files in the folder
matFiles = dir(fullfile(folderPath, '*.mat'));

% Initialize a structure to store merged data
dataMatrix = [];
columnNames = {};

% Loop through each .mat file and load the data
for i = 1:length(matFiles)
    % Get the basename of the file (without extension)
    [~, baseName, ~] = fileparts(matFiles(i).name);
    
    % Load the .mat file
    data = load(fullfile(folderPath, matFiles(i).name));
    
    % Assuming each .mat file contains a single variable
    fieldNames = fieldnames(data);
    if length(fieldNames) == 1
        % Extract the data from the structure
        rowData = data.(fieldNames{1});
        
        % Ensure the data is numeric and convert if necessary
        if istable(rowData)
            % Convert table to array if it's a table
            rowData = table2array(rowData);
        elseif isstruct(rowData)
            % Extract numeric values from a struct if possible
            rowData = struct2array(rowData);
        end
        
        % Ensure the rowData is a row vector
        if iscolumn(rowData)
            rowData = rowData';
        end
        
        % Append the rowData as a column in the data matrix
        dataMatrix = [dataMatrix; rowData]; % Collect rows to build columns later
        
        % Add the basename to the columnNames array
        columnNames{end+1} = baseName;
        
        % Add the basename as a row name (line title)
        rowNames = fieldnames(data.(fieldNames{1}));
    else
        error('File %s contains more than one variable.', matFiles(i).name);
    end
end

% Transpose the data matrix to make it column-major
dataMatrix = dataMatrix';

% Convert to a table with column names
mergedTable = array2table(dataMatrix, 'VariableNames', columnNames, 'RowNames', rowNames(1:11));

% Save the merged table to a .mat file
save(fullfile(folderPath, 'merged_data_with_titles.mat'), 'mergedTable');


%% Thresholded Similarity Metrics
% after statistical test using p<0.05 FWE corrected level

meth1 = 1;
meth2 = 2;

flagBinary = 1; % make the matrix binary
all_binThrMatrix = cell(nsmoo,1);
for i = 1:nsmoo % orig_TWS is not relevant here
    binThrMatrix = cell(nCombi,1);
    for ii = 1:nCombi
        GLM_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM_rft1',smoo_approachs{i}),combinations{ii});
        binThrMatrix{ii} = aj_get_thrMatrix(GLM_dir, flagBinary);
    end
    all_binThrMatrix{i} = binThrMatrix;
end

contrast_labels = {'T1', 'T2', 'F3'};

% Initialize variables
all_simMetrics = [];
all_row_titles = [];

% Loop through combinations
for i = 1:nCombi
    % Get the binary threshold matrices for the current combination
    binMatrices1 = all_binThrMatrix{meth1}{i};
    binMatrices2 = all_binThrMatrix{meth2}{i};
    
    % Check that dimensions match
    if ~isequal(size(binMatrices1), size(binMatrices2))
        error('The dimensions of binary threshold matrices in combination %d do not match.', i);
    end
    
    for ii = 1:length(binMatrices1)
        row_title = sprintf('%s_%s',combinations{i},contrast_labels{ii});
        
        Y1 = binMatrices1{ii};
        Y2 = binMatrices2{ii};
        [J, D, K] = aj_compute_simMetrics(Y1, Y2);
        all_simMetrics = [all_simMetrics; J, D, K];
        all_row_titles = [all_row_titles; {row_title}];
    end
end

% Convert metrics to a table
thrSimMetrics = array2table(all_simMetrics, ...
    'VariableNames', {'Jaccard', 'Dice', 'CohenKappa'}, ...
    'RowNames', all_row_titles);

% Save the table as .mat file
pth_out = fullfile(paths.ds_dir,'derivatives','qMRI_results','rft1','SimMetrics');
output_file = fullfile(pth_out, sprintf('thrSimMetrics_%s_%s.mat',smoo_approachs{meth1},smoo_approachs{meth2}));
save(output_file, 'thrSimMetrics');

fprintf('Similarity metrics table saved to %s\n', output_file);

%% NOT RELEVANT ? - Similarity Metrics
% orig_TWS is not relevant here
SM1_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM',smoo_approachs{2}));
SM2_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM',smoo_approachs{3}));

% Define contrasts
contrasts = {'spmT_0001.nii', 'spmT_0002.nii', 'spmF_0003.nii'};
contrast_labels = {'T1', 'T2', 'F3'};

% Initialize variables
all_simMetrics = [];
all_row_titles = [];

for i = 1:length(contrasts)
    % Get file lists
    SM1_files = spm_select('FPListRec', SM1_dir, ['^', contrasts{i}, '$']);
    SM2_files = spm_select('FPListRec', SM2_dir, ['^', contrasts{i}, '$']);

    if size(SM1_files, 1) ~= size(SM2_files, 1)
        error('Not the same number of %s files in the folders', contrasts{i});
    end

    % Process each pair of files
    for ii = 1:size(SM1_files, 1)
        % Extract file parts for row title
        [~, Callaghan_parent, ~] = fileparts(fileparts(SM1_files(ii, :))); % same for both
        row_title = sprintf('%s_%s', Callaghan_parent, contrast_labels{i});

        % Compute similarity metrics
        Y1 = spm_read_vols(spm_vol(SM1_files(ii, :)));
        Y2 = spm_read_vols(spm_vol(SM2_files(ii, :)));
        [J, D, K] = aj_compute_simMetrics(Y1, Y2);

        % Append results
        all_simMetrics = [all_simMetrics; J, D, K]; %#ok<AGROW>
        all_row_titles = [all_row_titles; {row_title}]; %#ok<AGROW>
    end
end

% Convert metrics to a table
simMetrics = array2table(all_simMetrics, ...
    'VariableNames', {'Jaccard', 'Dice', 'CohenKappa'}, ...
    'RowNames', all_row_titles);

% Save the table as .mat file
pth_out = fullfile(paths.ds_dir,'derivatives','qMRI_results','Index');
output_file = fullfile(pth_out, 'simMetrics.mat');
save(output_file, 'simMetrics');

fprintf('Similarity metrics table saved to %s\n', output_file);

%% Bland-Altman plot
% TWsmoot is not relevant here -> subject-specific denominator
meth1 = 2;
meth2 = 3;

SM1_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM_rft1',smoo_approachs{meth1}));
SM2_dir = fullfile(paths.ds_dir,'derivatives',sprintf('AJ-%s_GLM_rft1',smoo_approachs{meth2}));
contrasts = {'spmT_0001.nii', 'spmT_0002.nii', 'spmF_0003.nii'};

% Create regular expressions for all contrasts
contrast_patterns = strcat('^', contrasts, '$');

% Use spm_select to get all files for all contrasts in one call
SM1_lists = cellfun(@(pattern) spm_select('FPListRec', SM1_dir, pattern), ...
                              contrast_patterns, 'UniformOutput', false);
SM2_lists = cellfun(@(pattern) spm_select('FPListRec', SM2_dir, pattern), ...
                           contrast_patterns, 'UniformOutput', false);

% Convert cell array of char arrays to a single char array
SM1_list = vertcat(SM1_lists{:});
SM2_list = vertcat(SM2_lists{:});

% Count the number of contrasts
SM1_nContrasts = size(SM1_list, 1);
SM2_nContrasts = size(SM2_list, 1);
if SM1_nContrasts ~= SM2_nContrasts
    error('The two lists of paths have not the same length.');
end

% Convert multi-row char arrays to cell arrays of strings
SM1_list_cells = cellstr(SM1_list);
SM2_list_cells = cellstr(SM2_list);

% Extract the relevant parts of the paths
SM1_relPaths = cellfun(@(x) extractAfter(x, 'GLM_rft1\'), SM1_list_cells, 'UniformOutput', false);
SM2_relPaths = cellfun(@(x) extractAfter(x, 'GLM_rft1\'), SM2_list_cells, 'UniformOutput', false);

% Compare the relative paths
are_relative_paths_same = isequal(sort(SM1_relPaths), sort(SM2_relPaths));

if ~are_relative_paths_same
    error('The relative paths are different.');
else
    basename_out = cellfun(@(x) erase(x, '.nii'), SM1_relPaths, 'UniformOutput', false);
    basename_out = cellfun(@(x) strrep(x, '\', '_'), basename_out, 'UniformOutput', false);
    basename_out = cellfun(@(x) sprintf('%s-%s_%s', smoo_approachs{meth1}, smoo_approachs{meth2}, x), basename_out, 'UniformOutput', false);
    plot_titles = cellfun(@(x) extractBefore(x, '\spm'), SM1_relPaths, 'UniformOutput', false);
    plot_titles = cellfun(@(x) strrep(x, '_', '-'), plot_titles, 'UniformOutput', false);
    plot_titles = cellfun(@(x) sprintf('%s-%s:%s', smoo_approachs{meth1}, smoo_approachs{meth2}, x), plot_titles, 'UniformOutput', false);
end

flag.drawPlot = 0;
flag.savePlot = 1;
all_mean_diff = cell(SM1_nContrasts,1);
all_std_diff = cell(SM1_nContrasts,1);
for i = 1:SM1_nContrasts
    matrix3D_1 = spm_read_vols(spm_vol(SM1_list(i,:)));
    matrix3D_2 = spm_read_vols(spm_vol(SM2_list(i,:)));
    pth_out = fullfile(paths.ds_dir,'derivatives','qMRI_results','rft1','BlandAltman',sprintf('%s.eps',basename_out{i}));
    plot_title = char(plot_titles{i});
    [all_mean_diff{i},all_std_diff{i}] = aj_BlandAltman(matrix3D_1, matrix3D_2, flag, pth_out, plot_title);
end
