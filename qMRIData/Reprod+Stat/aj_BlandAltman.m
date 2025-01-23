function [mean_diff,std_diff] = aj_BlandAltman(matrix3D_1, matrix3D_2, flag, pth_out, plot_title)
%--------------------------------------------------------------------------
% Function to compute the Bland–Altman statistical methods.
% The Bland–Altman statistical method is used to assess agreement between 2
% methods of clinical measurement. In qMRI GLM constrasts context, the 
% measurement has to be discussed -> which data are relevant?
% Here relevant data are defined as the space region where a voxel has a at
% least one value (not zero, NaN or Inf) in one of the two contrast
% matrices.
%
% INPUTS
% matrix3D_1:   A 3D matrix containing the first set of measurements.
% matrix3D_2:   A 3D matrix containing the second set of measurements.
% flag:         (optional) a structure with flags that control whether a 
%               plot is drawn or saved.
% pth_out:      (optional) the output path where the Bland-Altman plot will
%               be saved, if flag.savePlot is set to 1. If not specified, 
%               the current working directory is used.
% plot_title:   (optional) just a char chain for plot title.
%
% OUTPUTS
% mean_diff:    The mean of the differences between the two sets of measurements.
% std_diff:     The standard deviation of the differences between the two sets.
%
% PROCESS
% The script flattens the 3D matrices into 1D vectors and removes the values 
% that are NaN, Inf, or zero in either of the two input matrices. It then 
% calculates the mean and difference for each pair of corresponding values.
% Bland-Altman statistics are computed:
%     - Mean difference
%     - Standard deviation of the differences
%     - Upper and lower limits of agreement (mean ± 1.96*std)
%
% If the flag.drawPlot is set to 1, a Bland-Altman plot will be generated 
% showing the differences between the measurements and the limits of agreement.
% If the flag.savePlot is set to 1, the plot will be saved to the specified 
% path in high resolution (PNG format).
%
% REFERENCE
% Bland JM, Altman DG (1986), https://doi.org/10.1016/S0140-6736(86)90837-8
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Dealing with inputs
if nargin<2
    error('Not enough inputs.');
end
if nargin<3 || isempty(flag)
    flag.drawPlot = 0;
    flag.savePlot = 0;
end
if nargin<4
    pth_out=pwd;
%     fprintf('Default output path: %s', pth_out);
end
if nargin<5
    plot_title = 'Bland-Altman Plot';
end
if size(matrix3D_1) ~= size(matrix3D_2)
    error('The two input matrices have not the same size.');
end

% Flatten volumes into 1D arrays
full_vector1 = matrix3D_1(:);
full_vector2 = matrix3D_2(:);

%% Do the job
% Identify indices where at least one of the two vectors are NaN or Inf
nan_inf_indices1 = find(isnan(full_vector1) | isinf(full_vector1));
nan_inf_indices2 = find(isnan(full_vector2) | isinf(full_vector2));
nan_inf_indices = unique([nan_inf_indices1; nan_inf_indices2]);

% Identify indices where both vectors are zero
zero_indices1 = find(full_vector1 == 0);
zero_indices2 = find(full_vector2 == 0);
% common_zero_indices = intersect(zero_indices1, zero_indices2);
all_zero_indices = unique([zero_indices1; zero_indices2]);

% Points of the forms (0, x) or (y, 0)
% zeroOnlyInVector1 = setdiff(zero_indices1, zero_indices2);
% disp(numel(zeroOnlyInVector1));
% zeroOnlyInVector2 = setdiff(zero_indices2, zero_indices1);
% disp(numel(zeroOnlyInVector2));

% Keep indices that are NOT in the remove_indices
% remove_indices = unique([nan_inf_indices; common_zero_indices]);
remove_indices = unique([nan_inf_indices; all_zero_indices]);
% keep_indices = setdiff(all_zero_indices, remove_indices);
keep_indices = setdiff(1:length(full_vector1), remove_indices);
vector1 = full_vector1(keep_indices);
vector2 = full_vector2(keep_indices);

% Compute the mean and the difference for each pair of values
mean_vector = (vector1+vector2)/2;
diff_vector = vector1-vector2;

% Compute Bland–Altman statistics
mean_diff = mean(diff_vector);
std_diff = std(diff_vector);
upper_limit = mean_diff + 1.96 * std_diff;
lower_limit = mean_diff - 1.96 * std_diff;

if flag.drawPlot || flag.savePlot
    % Plot the Bland-Altman data points
    figureHandle = figure;
    scatter(mean_vector, diff_vector, 2, [0.5, 0.5, 0.5], 'filled');
    hold on;

    % Plot the mean difference (blue dashed line)
    yline(mean_diff, 'b--', 'LineWidth', 1.5);
    % Text above the mean difference line (name)
    text(max(mean_vector), mean_diff + 0.1, 'Mean difference', 'Color', [0, 0, 1], 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    % Text below the mean difference line (value)
    text(max(mean_vector), mean_diff - 0.1, sprintf('%f', mean_diff), 'Color', [0, 0, 1], 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    % Plot the upper limit of agreement (red dashed line)
    yline(upper_limit, 'r--', 'LineWidth', 1.5);
    % Text above the upper limit line (name)
    text(max(mean_vector), upper_limit + 0.1, '+1.96 SD', 'Color', [1, 0, 0], 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    % Text below the upper limit line (value)
    text(max(mean_vector), upper_limit - 0.1, sprintf('%f', upper_limit), 'Color', [1, 0, 0], 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    % Plot the lower limit of agreement (red dashed line)
    yline(lower_limit, 'r--', 'LineWidth', 1.5);
    % Text above the lower limit line (name)
    text(max(mean_vector), lower_limit + 0.1, '-1.96 SD', 'Color', [1, 0, 0], 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    % Text below the lower limit line (value)
    text(max(mean_vector), lower_limit - 0.1, sprintf('%f', lower_limit), 'Color', [1, 0, 0], 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    % Customize the plot
    xlabel('Average of paired measurements');
    ylabel('Difference of paired measurements');
    title(plot_title);
    legend({sprintf('Data Points (n = %d)', length(mean_vector)), ...
            'Mean Difference','Upper Limit','Lower Limit'},...
           'Location', 'best');

    grid on;
    xlim([min(mean_vector) - 1, max(mean_vector) + 1]);
    ylim([min(diff_vector) - 1, max(diff_vector) + 1]);
    hold off;
    
    if flag.savePlot
        exportgraphics(figureHandle, pth_out, 'Resolution', 450); % file format defined in pth_out
    end
    if ~flag.drawPlot
        close(gcf);
    end
end

%% Fit a linear regression to the difference vs mean data
% p = polyfit(mean_vector, diff_vector, 1); % 1st degree polynomial (linear fit)
% 
% hold on;
% plot(mean_vector, polyval(p, mean_vector), 'k--', 'LineWidth', 1.5); % Black dashed line for linear fit
% hold off;
% 
% disp(['Slope of the regression line: ', num2str(p(1))]);

%% Find the index of the closest point to the given (mean_value, diff_value)
% mean_value = -2.78492;
% diff_value = 5.56984;
% 
% [~, index_closest] = min(sqrt((mean_vector - mean_value).^2 + (diff_vector - diff_value).^2));
% original_index = keep_indices(index_closest); 
% 
% disp(['The original index corresponding to the point (', num2str(mean_value), ', ', num2str(diff_value), ') is:']);
% disp(original_index);
% fprintf('Original values in the matrices: (%d ; %d)\n',full_vector1(original_index),full_vector2(original_index));

end