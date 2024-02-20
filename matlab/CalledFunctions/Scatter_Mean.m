function Scatter_Mean( TS_mean_preVpost, TS_idx, TITLE, XLABEL, YLABEL, stdOrsem )
% Michelle Cahill 20230814
% Scatter_Mean creates a scatter plot from values in TS_mean_preVpost and
% plots the mean +- errorbars for the population

%   TS_mean_preVpost: nx2 matrix with n= # of tseries/pairs, column 1:
%   metric pre-stim column 2: metric post-stim
%   TS_idx: the indices of t-series to plot (identified pairs)
%   TITLE: title for the plot (as a string)
%   YLABEL: ylabel for the plot (as a string)
%   YLIMS: ylims for the plot [yMin, yMax]
%   stdOrsem: string 'std' or 'sem' to determine the type of error bars
MAX = max(TS_mean_preVpost, [], 'all');
pre_mean = mean(TS_mean_preVpost(TS_idx,1), 'omitnan');
post_mean = mean(TS_mean_preVpost(TS_idx,2), 'omitnan');
if strcmp(stdOrsem, 'std')
    pre_EB = nanstd(TS_mean_preVpost(TS_idx,1));
    post_EB = nanstd(TS_mean_preVpost(TS_idx,2));
elseif strcmp(stdOrsem, 'sem')
    pre_EB = nanstd(TS_mean_preVpost(TS_idx,1))/(sqrt(length(TS_idx)-1));
    post_EB = nanstd(TS_mean_preVpost(TS_idx,2))/(sqrt(length(TS_idx)-1));
end

figure()
hold on
scatter(TS_mean_preVpost(TS_idx,1),TS_mean_preVpost(TS_idx,2), 'k', 'filled', 'MarkerFaceAlpha', 0.5)
plot([0 MAX], [0, MAX], 'k--')

X_mean = errorbar(pre_mean, post_mean, pre_EB, 'horizontal');
X_mean.Marker = 'o';
X_mean.Color = 'red';
X_mean.MarkerFaceColor = 'red';
X_mean.MarkerSize = 10;

Y_mean = errorbar(pre_mean, post_mean, post_EB, 'vertical');
Y_mean.Marker = 'o';
Y_mean.Color = 'red';
Y_mean.MarkerFaceColor = 'red';
Y_mean.MarkerSize = 10;

if strcmp(stdOrsem, 'std')
    title([TITLE, 'mean+-std'])
elseif strcmp(stdOrsem, 'sem')
    title([TITLE, 'mean+-sem'])
end

xlabel(XLABEL)
ylabel(YLABEL)

clear pre_mean post_mean p_RuBiGABA str post_x post_y pre_x pre_y
end