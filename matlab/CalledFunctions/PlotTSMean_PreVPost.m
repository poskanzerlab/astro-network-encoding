function PlotTSMean_PreVPost( TS_mean_preVpost, TS_idx, TITLE, YLABEL, YLIMS, stdOrsem )
% Michelle Cahill 20200601
%PlotTSMean_PreVPost plots the mean value of each t-series pre-stim and
%   post-stim & the mean+-std or +- sem of all t-series. A paired t-test is carried
%   out for identified pre/post pairs 
%
%   TS_mean_preVpost: nx2 matrix with n= # of tseries/pairs, column 1:
%   metric pre-stim column 2: metric post-stim
%   TS_idx: the indices of t-series to plot (identified pairs)
%   TITLE: title for the plot (as a string)
%   YLABEL: ylabel for the plot (as a string)
%   YLIMS: ylims for the plot [yMin, yMax]
%   stdOrsem: string 'std' or 'sem' to determine the type of error bars

figure()
[~, p_val] = ttest(TS_mean_preVpost(TS_idx,1), TS_mean_preVpost(TS_idx,2)); %paired t-test
hold on
plot(TS_mean_preVpost(TS_idx,:)', 'k-o', 'MarkerFaceColor', 'k')
if strcmp(stdOrsem, 'std')
    pre_mean = errorbar(.8, nanmean(TS_mean_preVpost(TS_idx,1)), nanstd(TS_mean_preVpost(TS_idx,1)), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 2);
    post_mean = errorbar(2.2, nanmean(TS_mean_preVpost(TS_idx,2)), nanstd(TS_mean_preVpost(TS_idx,2)), '-o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 2);
    pre_EB = nanstd(TS_mean_preVpost(TS_idx,1));
    post_EB = nanstd(TS_mean_preVpost(TS_idx,2));
elseif strcmp(stdOrsem, 'sem')
    pre_mean = errorbar(.8, nanmean(TS_mean_preVpost(TS_idx,1)), nanstd(TS_mean_preVpost(TS_idx,1))/(sqrt(length(TS_idx)-1)), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 2);
    post_mean = errorbar(2.2, nanmean(TS_mean_preVpost(TS_idx,2)), nanstd(TS_mean_preVpost(TS_idx,2))/(sqrt(length(TS_idx)-1)), '-o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 2);
    pre_EB = nanstd(TS_mean_preVpost(TS_idx,1))/(sqrt(length(TS_idx)-1));
    post_EB = nanstd(TS_mean_preVpost(TS_idx,2))/(sqrt(length(TS_idx)-1));
end

% text(.45, nanmean(TS_mean_preVpost(TS_idx,1)), sprintf('%.3g',nanmean(TS_mean_preVpost(TS_idx,1))))
% text(2.35, nanmean(TS_mean_preVpost(TS_idx,2)), sprintf('%.3g',nanmean(TS_mean_preVpost(TS_idx,2))))
text(.45, nanmean(TS_mean_preVpost(TS_idx,1)), {sprintf('%.3g',nanmean(TS_mean_preVpost(TS_idx,1))), sprintf('%.3g',pre_EB)})
text(2.35, nanmean(TS_mean_preVpost(TS_idx,2)), {sprintf('%.3g',nanmean(TS_mean_preVpost(TS_idx,2))), sprintf('%.3g',post_EB)})

str = sprintf('p = %d', p_val);
annotation('textbox',[.4 .5 .3 .3],'String',str,'FitBoxToText','on', 'EdgeColor', [1 1 1]);

if strcmp(stdOrsem, 'std')
    title([TITLE, 'mean+-std'])
elseif strcmp(stdOrsem, 'sem')
    title([TITLE, 'mean+-sem'])
end
legend([pre_mean, post_mean], {'pre-stim', 'post-stim'})
ylabel(YLABEL)
xticks([1 2])
xticklabels({'pre-stim', 'post-stim'})
if ~isempty(YLIMS)
    axis([0 3 YLIMS(1) YLIMS(2)])
else
    xlim([0 3])
end
clear pre_mean post_mean p_RuBiGABA str post_x post_y pre_x pre_y

end

% cd(save_dir)
% saveas(gcf, sprintf('Cyto_RuBiGABA2P_%dpulse_TempCoOccur_TTest_%dS_%d.tif',pulse_num, time_window, date))
% saveas(gcf, sprintf('Cyto_RuBiGABA2P_%dpulse_TempCoOccur_TTest_%dS_%d.svg',pulse_num, time_window, date))
% close
% cd(start_dir)

