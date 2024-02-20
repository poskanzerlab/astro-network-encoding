function LinePlot_SummaryStat_IndData(mat2plot, TS, TITLE, condition, YLIMS, YLABEL, SUMMARY, stdOrsem, COLOR)
%LinePlot_SummaryStat_IndData: 
% 1. Takes a matrix of values
% 2. Plots each row of the matrix as a grey line plot (repeated measures as columns, ex: multiple imaging rounds) 
% 3. Plots a summary statistic in red for each column: either mean+- errorbars or
%    median & IQR (ex. average value in each imaging round)
% 4. Run repeated measures ANOVA or Friedman test to compare values across
%    columns
% mat2plot: (r, c) matrix, in which rows are individual subject (ex. FOV) and columns are
%            repeated measures for that subject (ex. rounds of imaging)
% TS:        array indicating which rows from mat2plot to include 
% condition: cell array containing a label for each column in mat2plot
% SUMMARY:   'mean' or 'median'
% stdOrsem:  'std' or 'sem'; only relevant when SUMMARY == 'mean'
% COLOR:     'r' or any plot color to be the summary stat color
%
%   Michelle Cahill 20230927

mat2plot_selectTS = mat2plot(TS,:);

if strcmp(SUMMARY, 'mean') %if plotting mean as the summary statistic carry out a parametric test (repeated measures ANOVA)
    %Repeated measures ANOVA
    temp_table = array2table(mat2plot_selectTS);
    temp_table.Properties.VariableNames = condition;
    COLUMNS = [1:size(mat2plot_selectTS, 2)];
    MODEL_SPECS = sprintf('%s-%s ~1', condition{1}, condition{end});
    RM = fitrm(temp_table, MODEL_SPECS, 'WithinDesign', COLUMNS);
    Rounds_ranova = ranova(RM);
%     Rounds_multcomp = multcompare(RM, 'Time');
elseif strcmp(SUMMARY, 'median') %if plotting median as the summary statistic carry out a non-parametric test (equivalent to repeated measures ANOVA)
    %Friedman Test (non-parametric)
    [p_val, ~, stats] = friedman(mat2plot_selectTS, 1);
end

figure()
hold on
plot([1:size(mat2plot_selectTS, 2)], mat2plot_selectTS, 'Color', [0.75 0.75 0.75]) %individual data will be plotted in the background as grey

if strcmp(SUMMARY, 'mean')
    summary_data = mean(mat2plot_selectTS, 1, 'omitnan');
    SD_data = std(mat2plot_selectTS, 0, 1, 'omitnan');
    sem_data = SD_data ./ (sqrt(size(mat2plot_selectTS, 1) - 1));
    if strcmp(stdOrsem, 'std')
        errorbar(summary_data,  SD_data, COLOR, 'Marker', 'o', 'MarkerFaceColor', COLOR, 'MarkerEdgeColor', 'none')
        tempTITLE = [TITLE, 'mean+-std'];
    elseif strcmp(stdOrsem, 'sem')
        errorbar(summary_data,  sem_data, COLOR, 'Marker', 'o', 'MarkerFaceColor', COLOR, 'MarkerEdgeColor', 'none')
        tempTITLE = [TITLE, 'mean+-sem'];
    end
    tempTITLE = [tempTITLE, sprintf('F(%d,%d)=%.2f', Rounds_ranova.DF(1), Rounds_ranova.DF(2), Rounds_ranova.F(1)), sprintf('RM ANOVA p: %.15f', Rounds_ranova.pValue(1))];
elseif strcmp(SUMMARY, 'median')
    summary_data = median(mat2plot_selectTS, 1, 'omitnan');
    IQR = quantile(mat2plot_selectTS, [.25 .75]); %25th (row 1) and 75th (row 2) percentile for each column of mat2plot 
    plot(summary_data, COLOR, 'Marker', 'o', 'MarkerFaceColor', COLOR, 'MarkerEdgeColor', 'none')
    plot([1:size(mat2plot_selectTS, 2);1:size(mat2plot_selectTS, 2)], IQR, 'r')
    tempTITLE = [TITLE, 'median, 25th & 75th percentile'];
    tempTITLE = [tempTITLE, sprintf('Friedman p: %.15f', p_val)];
end

XTICKLABELS = cellfun(@(x) strrep(x, '_', ' '), condition, 'UniformOutput', 0);

title(tempTITLE)

xlim([0 size(mat2plot_selectTS, 2)+1])
xticks([1:size(mat2plot_selectTS, 2)])
xticklabels(XTICKLABELS)
if ~isempty(YLIMS)
    ylim(YLIMS)
end
ylabel(YLABEL)
end
