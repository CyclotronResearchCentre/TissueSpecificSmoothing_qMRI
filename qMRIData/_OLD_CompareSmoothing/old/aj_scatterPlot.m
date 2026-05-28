function aj_scatterPlot(vector1, vector2, data_info, pth_out)
% This script computes a scatter plot built on two vectors. This plot
% is then saved.
%
% INPUT
% vector1:  Vector containing the numbers to plot according to the x-axis.
% vector2:  Vector containing the numbers to plot according to the y-axis.
% data_info:    List of characters containing the two names/types of
%               provided data in vector 1 and 2.
% pth_out:  Character chain containg the output path to save the plot.
%
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Dealing with inputs
if nargin < 3
    error('Not enough inputs.');
end
if nargin < 4
    pth_out = [];
    disp('No figure will be saved (no input pth_out).');
end
if ~isvector(vector1) || ~isvector(vector2)
    error('At least one input is not a vector.');
end

% Names for each data point
pointNames = {'MTsat GM', 'MTsat WM', 'PD GM', 'PD WM', ...
              'R1 GM', 'R1 WM', 'R2star GM', 'R2star WM'};

if length(vector1) ~= length(pointNames) || length(vector2) ~= length(pointNames)
    error('The number of points in vector1 and vector2 must match the number of point names.');
end

%% Plot the data
minVal = min([vector1, vector2]);
maxVal = max([vector1, vector2]);
diff = (maxVal - minVal) / 3;
minVal_plot = minVal - diff;
maxVal_plot = maxVal + diff;

figureHandle = figure;
s = scatter(vector1, vector2, 30, 'r', 'filled');
hold on;

% Plot diagonal line -> identity line
id_line = plot([minVal_plot, maxVal_plot], [minVal_plot, maxVal_plot], 'k-', 'LineWidth', 1);

% Plot vertical dashed line to min value
plot([min(vector1), min(vector1)], [minVal_plot, min(vector2)], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off'); 
% Plot horizontal dashed line to min value
plot([minVal_plot, min(vector1)], [min(vector2), min(vector2)], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Plot vertical dashed line to max value
plot([max(vector1), max(vector1)], [minVal_plot, max(vector2)], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off'); 
% Plot horizontal dashed line to max value
plot([minVal_plot, max(vector1)], [max(vector2), max(vector2)], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Adjust axis limits
xlim([minVal_plot, maxVal_plot]);
ylim([minVal_plot, maxVal_plot]);

grid on;

legend([s, id_line], {'Data Points', 'Identity Line'}, 'Location', 'best');

% Add annotations for each data point
for i = 1:length(vector1)
    text(vector1(i), vector2(i), pointNames{i}, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'blue', ...
        'Margin', 1); % Position en bas à droite
end

% Set labels and title
xlabel([data_info.metric1, ' T-values']);
ylabel([data_info.metric2, ' T-values']);
title(['Comparison plot for ', data_info.metric1, ' vs ', data_info.metric2, ' T-values']);

% Save figure if a path is provided
if ~isempty(pth_out)
    exportgraphics(figureHandle, pth_out, 'Resolution', 300);  % High-resolution PNG
    fprintf('Figure saved in %s\n', pth_out);
end

end

%% With Interactive Legend
% function aj_scatterPlot(vector1,vector2,data_info,pth_out)
% %% Dealing with inputs
% if nargin<3
%     error('Not enough inputs.');
% end
% if nargin<4
%     pth_out=[];
%     disp('No figure will be saved (no input pth_out).');
% end
% if ~isvector(vector1) || ~isvector(vector2)
%     error('At least one input is not a vector.');
% end
% 
% %% Do the job
% minVal_x = min(vector1);
% minVal_y = min(vector2);
% maxVal_x = max(vector1);
% maxVal_y = max(vector2);
% 
% minVal = min([vector1, vector2]);
% maxVal = max([vector1, vector2]);
% diff = (maxVal - minVal)/3;
% minVal_plot = minVal-diff;
% maxVal_plot = maxVal+diff;
% 
% figureHandle = figure;
% s = scatter(vector1, vector2, 30, 'r', 'filled');
% hold on;
% 
% % Plot diagonal line -> identity line
% id_line = plot([minVal_plot, maxVal_plot], [minVal_plot, maxVal_plot], 'k-', 'LineWidth', 1);
% 
% % Plot vertical dashed line to min value
% plot([minVal_x, minVal_x], [minVal_plot, minVal_y], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off'); 
% % Plot horizontal dashed line to min value
% plot([minVal_plot, minVal_x], [minVal_y, minVal_y], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
% 
% % Plot vertical dashed line to max value
% plot([maxVal_x, maxVal_x], [minVal_plot, maxVal_y], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off'); 
% % Plot horizontal dashed line to max value
% plot([minVal_plot, maxVal_x], [maxVal_y, maxVal_y], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
% 
% % Adjust axis limits to fit the data+diff
% xlim([minVal_plot, maxVal_plot]);
% ylim([minVal_plot, maxVal_plot]);
% 
% grid on;
% 
% legend([s, id_line], {'Data Points', 'Identity Line'}, 'Location', 'best');
% 
% % Add interactive datatips showing the qMRI metrics of each point
% dcm = datacursormode; % Enable data cursor mode 
% set(dcm, 'Enable', 'on');
% 
% % Customize the data tip 
% % Create a custom update function for the data tip 
% function txt = customDataTip(~, event_obj) 
%     % Get the position of the data cursor 
%     pos = get(event_obj, 'Position');
%     % Find the index of the closest data point 
%     [~, index] = min(abs(vector1 - pos(1)));
%     combi = char('MTsat GM', 'MTsat WM','PDmap GM','PDmap WM','R1map GM','R1map WM','R2starmap GM','R2starmap WM');
%     % Create the custom data tip text 
%     txt = {['Name: ', combi(index,:)], ...
%         [metric1 ': ', num2str(pos(1))], ... 
%         [metric2 ': ', num2str(pos(2), '%.2f')]
%         };
% end
% 
% set(dcm, 'UpdateFcn', @customDataTip); % Set the custom data tip function
% 
% metric1 = data_info.metric1;
% xlabel([metric1 ' T-values']);
% metric2 = data_info.metric2;
% ylabel([metric2 ' T-values']);
% contrast = data_info.contrast;
% title(['Comparison plot for ' metric1 ' vs ' metric2 ' T-values']);
% % title(['Comparison plot for ' metric1 ' vs ' metric2 ' T-values of ' contrast ' contrast']);
% 
% if ~isempty(pth_out)
%     exportgraphics(figureHandle, pth_out, 'Resolution', 300);  % High-resolution PNG
% %     exportgraphics(figureHandle, pth_out);  % EPS format
%     fprintf('Figure saved in %s\n', pth_out);
% end
% 
% end
