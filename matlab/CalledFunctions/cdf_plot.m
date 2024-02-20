function cdf_plot(val, linecolors, XLABEL, TITLE, LEGEND_LABELS)
% cdf_plot takes an array of values and calculates the cumulative
% probability of that parameter. MC 20201217
% val: 1xn cell array of where each column of the cell array holds the
% values to calculate the CDF for a single condition
% linecolors: nx3 matrix for the colors to plot each condition in

min_edge = floor(min(cellfun((@min), val)));
max_edge = ceil(max(cellfun((@max), val)));

figure()
hold on
for cond = 1:size(val,2)
    bin_edges = min_edge:.1:max_edge;
    [val_count, val_edges] = histcounts(val{1,cond}, bin_edges);
    val_pdf = val_count ./ sum(val_count);
    val_edges = val_edges(1:end-1);
    val_cdf = cumsum(val_count);
    val_cdf = val_cdf ./ val_cdf(end);
    plot(val_edges, val_cdf, 'Color', linecolors(cond,:), 'LineWidth', 1.5);
end
xlabel(XLABEL)
ylabel('cumulative probability')
title(TITLE)
legend(LEGEND_LABELS)
% xlabel('event area (um^2)')
% ylabel('percent of total count')
% title(sprintf('Area: 100uM Baclofen tseries %d-%d', first_tseries, last_tseries))
% legend({'prestim', 'poststim'})
% set(gca, 'xscale', 'log')
end