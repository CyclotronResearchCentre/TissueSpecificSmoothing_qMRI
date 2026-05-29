function [results, excel_rows] = aj_get_BlandAltman_compute( ...
    file1, file2, combi, ...
    method1_name, method2_name, basename, title_name, ...
    flag, out_dir, mask_thr, excel_rows)

% Load maps
Y1 = spm_read_vols(spm_vol(file1));
Y2 = spm_read_vols(spm_vol(file2));

% Flatten
vec1 = Y1(:);
vec2 = Y2(:);

% Apply anatomical mask (GM / WM)
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
    fprintf('\nNo mask provided, skipping masking step\n');
end

% Remove invalid voxels
valid = ...
    ~isnan(vec1) & ~isnan(vec2) & ...
    ~isinf(vec1) & ~isinf(vec2) & ...
    vec1 ~= 0 & vec2 ~= 0;

% Add threshold mask if requested
if flag.thresholded
    
    if isempty(mask_thr)
        error('mask_thr is required when thresholded = 1');
    end

    mask_thr = mask_thr(:);
    valid = valid & mask_thr;
end

vec1 = vec1(valid);
vec2 = vec2(valid);

nPoints = numel(vec1);

% Bland–Altman

mean_vals = (vec1 + vec2)/2;
diff_vals = vec1 - vec2;

mean_diff = mean(diff_vals);
std_diff  = std(diff_vals);

loa_upper = mean_diff + 1.96*std_diff;
loa_lower = mean_diff - 1.96*std_diff;

% Linear fit
p = polyfit(mean_vals,diff_vals,1);

xfit = linspace(min(mean_vals),max(mean_vals),100);
yfit = polyval(p,xfit);

% Save results
results.(combi).mean_diff = mean_diff;
results.(combi).std_diff  = std_diff;
results.(combi).loa_upper = loa_upper;
results.(combi).loa_lower = loa_lower;
results.(combi).fit       = p;
results.(combi).nPoints   = nPoints;

% Plot
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
    title(title_name);
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

% Save figures
if flag.savePlot

    png_out = fullfile(out_dir,[basename '.png']);
    eps_out = fullfile(out_dir,[basename '.eps']);

    exportgraphics(figureHandle,png_out,'Resolution',1000);

    exportgraphics(figureHandle,eps_out);

end

% Excel rows
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