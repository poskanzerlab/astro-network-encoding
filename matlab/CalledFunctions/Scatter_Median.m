function Scatter_Median( TS_mean_preVpost, TS_idx, TITLE, XLABEL, YLABEL)
% Michelle Cahill 20230814
% Scatter_Median creates a scatter plot from values in TS_mean_preVpost and
% plots the median +- iqr for the population

%   TS_mean_preVpost: nx2 matrix with n= # of tseries/pairs, column 1:
%   metric pre-stim column 2: metric post-stim
%   TS_idx: the indices of t-series to plot (identified pairs)
%   TITLE: title for the plot (as a string)
%   YLABEL: ylabel for the plot (as a string)
%   YLIMS: ylims for the plot [yMin, yMax]
MAX = max(TS_mean_preVpost, [], 'all');
pre_median = median(TS_mean_preVpost(TS_idx,1), 'omitnan');
post_median = median(TS_mean_preVpost(TS_idx,2), 'omitnan');

pre_qt = quantile(TS_mean_preVpost(TS_idx,1), [.25 .75]);
post_qt = quantile(TS_mean_preVpost(TS_idx,2), [.25 .75]);

% if strcmp(stdOrsem, 'std')
%     pre_EB = nanstd(TS_mean_preVpost(TS_idx,1));
%     post_EB = nanstd(TS_mean_preVpost(TS_idx,2));
% elseif strcmp(stdOrsem, 'sem')
%     pre_EB = nanstd(TS_mean_preVpost(TS_idx,1))/(sqrt(length(TS_idx)-1));
%     post_EB = nanstd(TS_mean_preVpost(TS_idx,2))/(sqrt(length(TS_idx)-1));
% end

figure()
hold on
scatter(TS_mean_preVpost(TS_idx,1),TS_mean_preVpost(TS_idx,2), 'k', 'filled', 'MarkerFaceAlpha', 0.5)
plot([0 MAX], [0, MAX], 'k--')

scatter(pre_median, post_median, 50, 'r', 'filled')
plot(pre_qt, [post_median, post_median], '-r', 'LineWidth', 1) %IQR for x-axis data
plot([pre_median, pre_median], post_qt, '-r', 'LineWidth', 1) %IQR for y-axis data

% 
% X_median = errorbar(pre_median, post_median, pre_EB, 'horizontal');
% X_median.Marker = 'o';
% X_median.Color = 'red';
% X_median.MarkerFaceColor = 'red';
% X_median.MarkerSize = 10;
% 
% Y_median = errorbar(pre_median, post_median, post_EB, 'vertical');
% Y_median.Marker = 'o';
% Y_median.Color = 'red';
% Y_median.MarkerFaceColor = 'red';
% Y_median.MarkerSize = 10;

title([TITLE, 'median, 25th & 75th percentile'])


xlabel(XLABEL)
ylabel(YLABEL)
axis square

clear pre_median post_median p_RuBiGABA str post_x post_y pre_x pre_y
end