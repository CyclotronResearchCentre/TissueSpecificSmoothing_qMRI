%% =========================================================================
%  Unified Bland–Altman Analysis for SPM / SnPM voxelwise comparisons
%
%  This function compares voxelwise statistical maps obtained using two
%  smoothing approaches and generates publication-ready Bland–Altman plots.
%
%  FEATURES
%  --------
%  - Supports:
%       * SnPM
%       * SPM RFT0
%       * SPM RFT1
%
%  - Automatic selection of:
%       * snpmT+
%       * snpmT-
%       * spmT_0001
%       * spmT_0002
%
%  - Tissue masking (GM / WM)
%  - Removal of:
%       * NaN
%       * Inf
%       * zero voxels
%
%  - Computes:
%       * mean difference
%       * SD difference
%       * limits of agreement
%       * linear regression
%       * number of voxels
%
%  - Saves:
%       * PNG (1000 dpi)
%       * EPS
%       * Excel summary
%
%  OUTPUT NAMING
%  -------------
%  ba_[spm_rft0/spm_rft1/snpm]_
%  [method1]_minus_[method2]_
%  [metric]_[GM/WM]
%
% =========================================================================

function results = aj_BlandAltman_unified( ...
    root_dir,...
    combination_names,...
    method1_name,...
    method2_name,...
    flag)

%% =========================================================================
% Defaults
%% =========================================================================

if nargin < 5
    flag.drawPlot  = 1;
    flag.savePlot  = 1;
    flag.saveExcel = 1;
    flag.paperReady = 1;
    flag.analysis = 'snpm';
    flag.thresholded = 0;
    flag.useLL = 0;   % 0 = Tmaps, 1 = LL_methodX
    flag.maskGM = '';
    flag.maskWM = '';
end

%% =========================================================================
% Define directories
%% =========================================================================

if flag.useLL == 0
    
    switch lower(flag.analysis)

        case 'snpm'
            base_dir1 = fullfile(root_dir,...
                sprintf('SnPM-%s_mreg-age',method1_name));

            base_dir2 = fullfile(root_dir,...
                sprintf('SnPM-%s_mreg-age',method2_name));

        case 'spm_rft0'
            base_dir1 = fullfile(root_dir,...
                sprintf('AJ-%s_GLM_rft0_age_linear',method1_name));

            base_dir2 = fullfile(root_dir,...
                sprintf('AJ-%s_GLM_rft0_age_linear',method2_name));

        case 'spm_rft1'
            base_dir1 = fullfile(root_dir,...
                sprintf('AJ-%s_GLM_rft1',method1_name));

            base_dir2 = fullfile(root_dir,...
                sprintf('AJ-%s_GLM_rft1',method2_name));

        otherwise
            error('Unknown analysis type.');
    end
    
elseif flag.useLL == 1
    
    method_map = struct( ...
        'SUSAN', 1, ...
        'TSPOON', 2, ...
        'TWS', 3);

    methodA = method_map.(upper(method1_name));
    methodB = method_map.(upper(method2_name));

    base_dir1 = fullfile(root_dir,flag.analysis);
    base_dir2 = fullfile(root_dir,flag.analysis);
else
    error('Unknown analysis type.');
end

%% =========================================================================
% Output dir
% =========================================================================

out_dir = fullfile(root_dir,...
    'results_BlandAltman',...
    sprintf('%s_minus_%s',method1_name,method2_name));

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% =========================================================================
% Init
% =========================================================================

results = struct();
excel_rows = {};

%% =========================================================================
% Loop
% =========================================================================

for i = 1:length(combination_names)

    combi = combination_names{i};

    fprintf('\nProcessing %s\n',combi);

    % ---------------------------------------------------------------------
    % Metric / tissue extraction
    % ---------------------------------------------------------------------

    parts = split(combi,'_');

    metric = strrep(parts{1},'map','');
    tissue = parts{2};

    %======================================================================
    % Select correct input type
    %======================================================================

    if flag.useLL == 1

        %==============================================================
        % LL METHOD FILES
        %==============================================================
        file1 = fullfile(base_dir1, combi, sprintf('LL_method%d.nii', methodA));
        file2 = fullfile(base_dir2, combi, sprintf('LL_method%d.nii', methodB));

    else

        %==============================================================
        % STAT MAPS (SnPM / SPM)
        %==============================================================

        if contains(combi, 'R2starmap') || contains(combi, 'PDmap')

            if strcmp(flag.analysis,'snpm')
                file1 = fullfile(base_dir1, combi, 'snpmT+.img');
                file2 = fullfile(base_dir2, combi, 'snpmT+.img');
            elseif strcmp(flag.analysis,'spm_rft0')
                file1 = fullfile(base_dir1, combi, 'spmT_0001.nii');
                file2 = fullfile(base_dir2, combi, 'spmT_0001.nii');
            elseif strcmp(flag.analysis,'spm_rft1')
                file1 = fullfile(base_dir1, combi, 'Tmap_thr_bin.nii');
                file2 = fullfile(base_dir2, combi, 'Tmap_thr_bin.nii');
            end

        elseif contains(combi, 'MTsat') || contains(combi, 'R1map')

            if strcmp(flag.analysis,'snpm')
                file1 = fullfile(base_dir1, combi, 'snpmT-.img');
                file2 = fullfile(base_dir2, combi, 'snpmT-.img');
            elseif strcmp(flag.analysis,'spm_rft0')
                file1 = fullfile(base_dir1, combi, 'spmT_0002.nii');
                file2 = fullfile(base_dir2, combi, 'spmT_0002.nii');
            elseif strcmp(flag.analysis,'spm_rft1')
                file1 = fullfile(base_dir1, combi, 'Tmap_thr_bin.nii');
                file2 = fullfile(base_dir2, combi, 'Tmap_thr_bin.nii');
            end

        else
            error('Unknown contrast type');
        end
    end
    
    %----------------------------------------------------------------------
    % Thresholded masks
    %----------------------------------------------------------------------
    mask_thr = [];

    if flag.thresholded == 1

        switch lower(flag.analysis)

            %==============================================================
            % SnPM
            %==============================================================
            case 'snpm'

                mask1_file = fullfile(base_dir1, combi, 'Tmap_thr_bin.nii');
                mask2_file = fullfile(base_dir2, combi, 'Tmap_thr_bin.nii');

            %==============================================================
            % SPM rft0
            %==============================================================
            case 'spm_rft0'

                d1 = dir(fullfile(base_dir1, combi, '*_FWE05_mask.nii'));
                d2 = dir(fullfile(base_dir2, combi, '*_FWE05_mask.nii'));

                if ~isempty(d1)
                    mask1_file = fullfile(d1(1).folder, d1(1).name);
                else
                    mask1_file = '';
                end

                if ~isempty(d2)
                    mask2_file = fullfile(d2(1).folder, d2(1).name);
                else
                    mask2_file = '';
                end

            %==============================================================
            % SPM rft1
            %==============================================================
            case 'spm_rft1'

                mask1_file = fullfile(base_dir1, combi, 'Tmap_thr_bin.nii');
                mask2_file = fullfile(base_dir2, combi, 'Tmap_thr_bin.nii');

            otherwise
                error('Unknown analysis type.');
        end

        %------------------------------------------------------------------
        % Check existence
        %------------------------------------------------------------------
        if ~exist(mask1_file,'file') || ~exist(mask2_file,'file')

            warning('Thresholded mask missing for %s → skipping', combi);
            continue;

        end

        %------------------------------------------------------------------
        % Load threshold masks
        %------------------------------------------------------------------
        M1 = spm_read_vols(spm_vol(mask1_file)) > 0;
        M2 = spm_read_vols(spm_vol(mask2_file)) > 0;

        % Common threshold mask
        mask_thr = M1 & M2;

    end

    %% ---------------------------------------------------------------------
    % Load maps
    % ---------------------------------------------------------------------

    Y1 = spm_read_vols(spm_vol(file1));
    Y2 = spm_read_vols(spm_vol(file2));

    %% ---------------------------------------------------------------------
    % Flatten
    % ---------------------------------------------------------------------

    vec1 = Y1(:);
    vec2 = Y2(:);
    
    %======================================================================
    % Apply anatomical mask (GM / WM)
    %======================================================================

    if contains(combi,'_GM')
        mask_fn = flag.maskGM;
    elseif contains(combi,'_WM')
        mask_fn = flag.maskWM;
    else
        error('Cannot determine tissue type (GM/WM) for %s', combi);
    end

    if ~isempty(mask_fn)

        Vmask = spm_vol(mask_fn);
        mask = spm_read_vols(Vmask) > 0;
        mask = mask(:);

        vec1 = vec1(mask);
        vec2 = vec2(mask);

    else
        warning('No mask provided, skipping masking step');
    end

    %% ---------------------------------------------------------------------
    % Remove invalid voxels
    % ---------------------------------------------------------------------

    valid = ...
        ~isnan(vec1) & ~isnan(vec2) & ...
        ~isinf(vec1) & ~isinf(vec2) & ...
        vec1 ~= 0 & vec2 ~= 0;

    %----------------------------------------------------------------------
    % Add threshold mask if requested
    %----------------------------------------------------------------------
    if flag.thresholded == 1

        mask_thr = mask_thr(:);

        valid = valid & mask_thr;

    end

    vec1 = vec1(valid);
    vec2 = vec2(valid);

    nPoints = numel(vec1);

    %% ---------------------------------------------------------------------
    % Bland–Altman
    % ---------------------------------------------------------------------

    mean_vals = (vec1 + vec2)/2;
    diff_vals = vec1 - vec2;

    mean_diff = mean(diff_vals);
    std_diff  = std(diff_vals);

    loa_upper = mean_diff + 1.96*std_diff;
    loa_lower = mean_diff - 1.96*std_diff;

    %% ---------------------------------------------------------------------
    % Linear fit
    % ---------------------------------------------------------------------

    p = polyfit(mean_vals,diff_vals,1);

    xfit = linspace(min(mean_vals),max(mean_vals),100);
    yfit = polyval(p,xfit);

    %% ---------------------------------------------------------------------
    % Save results
    % ---------------------------------------------------------------------

    results.(combi).mean_diff = mean_diff;
    results.(combi).std_diff  = std_diff;
    results.(combi).loa_upper = loa_upper;
    results.(combi).loa_lower = loa_lower;
    results.(combi).fit       = p;
    results.(combi).nPoints   = nPoints;

    %% ---------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------

    figureHandle = figure;

    if flag.paperReady
        scatter(mean_vals,diff_vals,4,...
            [0.8 0.8 0.8],'filled');
    else
        scatter(mean_vals,diff_vals,6,'filled');
    end

    hold on;

    h1 = yline(mean_diff,'k-','LineWidth',1.5);
    h2 = yline(loa_upper,'k--','LineWidth',1.2);
    h3 = yline(loa_lower,'k--','LineWidth',1.2);

    h4 = plot(xfit,yfit,'k:','LineWidth',2);

    xlabel('Mean of paired values');

    ylabel(sprintf('Difference of paired values (%s - %s)', ...
        method1_name,...
        method2_name));

    if ~flag.paperReady

        title(sprintf('%s: %s - %s (%s %s)',...
            upper(flag.analysis),...
            method1_name,...
            method2_name,...
            metric,...
            tissue));

    end

    legend( ...
        [h1 h2 h3 h4], ...
        { ...
        sprintf('\\mu = %.3f',mean_diff),...
        sprintf('+1.96SD = %.3f',loa_upper),...
        sprintf('-1.96SD = %.3f',loa_lower),...
        sprintf('Fit: y = %.3fx + %.3f | n = %d',...
        p(1),p(2),nPoints)...
        },...
        'Location','best');

    grid on;

    %% ---------------------------------------------------------------------
    % Output naming
    % ---------------------------------------------------------------------
    
    suffix_thr = '';

    if flag.thresholded == 1
        suffix_thr = '_thr';
    end
    
    suffix_mask = '';
    
    if flag.useLL
        data_type = 'LL';
        suffix_mask = '_mask';
    else
        data_type = flag.analysis;
    end

    basename = sprintf( ...
        'ba_%s_%s_minus_%s_%s_%s%s%s', ...
        data_type, ...
        method1_name, ...
        method2_name, ...
        metric, ...
        tissue, ...
        suffix_thr,...
        suffix_mask);

    %% ---------------------------------------------------------------------
    % Save figures
    % ---------------------------------------------------------------------

    if flag.savePlot

        png_out = fullfile(out_dir,[basename '.png']);
        eps_out = fullfile(out_dir,[basename '.eps']);

        exportgraphics(figureHandle,png_out,'Resolution',1000);

        exportgraphics(figureHandle,eps_out);

    end

    %% ---------------------------------------------------------------------
    % Excel rows
    % ---------------------------------------------------------------------

    excel_rows = [excel_rows;
        {basename,...
        mean_diff,...
        std_diff,...
        loa_upper,...
        loa_lower,...
        p(1),...
        p(2),...
        nPoints}];

    if ~flag.drawPlot
        close(gcf);
    end

end

%% =========================================================================
% Excel export
% =========================================================================

if flag.saveExcel

    T = cell2table(excel_rows,...
        'VariableNames',...
        {'Name',...
         'MeanDiff',...
         'StdDiff',...
         'UpperLoA',...
         'LowerLoA',...
         'Slope',...
         'Intercept',...
         'Npoints'});

    excel_out = fullfile(out_dir,...
        sprintf('BlandAltman_%s_minus_%s%s%s.xlsx',...
        method1_name,...
        method2_name,...
        suffix_thr,...
        suffix_mask));

    writetable(T,excel_out);

end

fprintf('\nBland–Altman analysis completed.\n');

end