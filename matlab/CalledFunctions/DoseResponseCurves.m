function [CORRECTED_STATS] = DoseResponseCurves(change_struct,tseries_idx_struct,TreatmentLevels, condition,...
    TITLE, YLABEL, YLIMS, XLABEL, XLIMS,...
    DATACOLORS, stdOrsem, LW, STATS, pairedOrunpaired)
% DoseResponseCurves plots line plots for each condition with mean +- error
% for each treatment level. Individual data points for each t-series are
% plotted behind mean values
% Michelle Cahill 20210903

% change_stuct: A structured array with fields for each condition. Within
% each field are the values per t-series or region (intended to be change
% values post-pre for a time bin defined outside this function). T-series
% here can be from multiple treatment levels
% tseries_idx_struct: A structured array with field names condition +
% treatment value ('Baclofen5'), each field has an mx1 matrix with the
% t-series indices for that condition & concentration
% TreatmentLevels: nx1 matrix with the doses/concentrations to plot
% condition: cx1 cell array with the conditions to include on the plot
% {'Baclofen'; 'tACPD'; 'Baclofen_CBX'; 'tACPD_CBX'}
%   TITLE: title for the plot (as a string)
%   YLABEL: ylabel for the plot (as a string)
%   YLIMS: ylims for the plot [yMin, yMax]
%   XLABEL: xlabel for the plot (as a string)
%   XLIMS: ylims for the plot [xMin, xMax] [0, length(pulse_num) + 1]
%   DATACOLORS: specification for how each condition will be plotted
%   {'-mo', '-co'}
%   stdOrsem: 'std' or 'sem' for the errorbars
%   LW: Width of the lines
%   STATS: 0 or 1 to perform paired t-test between conditions at each
%   treatment level and adjust for multiple comparisons with
%   Holm-Bonferroni correction
%   pairedOrunpaired: 'paired' to perform ttest (paired), 'unpaired' to
%   perform ttest2 (unpaired)

figure()
hold on
temp_change_perTV = struct();

for cond = 1:length(condition)
    temp_changeMeanError = nan(2, length(TreatmentLevels)); %row 1: treatment level means row 2:
    
    temp = DATACOLORS{cond};
    %         temp_color = temp(2);
    %         temp_linestyle = temp(1);
    %         temp_marker = temp(3);
    temp_color = extractBefore(extractAfter(temp, length(temp)-2), 2);
    temp_linestyle = extractBefore(temp, length(temp)-1);
    temp_marker = extractAfter(temp, length(temp)-1);
    
    for TV = TreatmentLevels
        TV_idx = find(TreatmentLevels == TV);
        concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
        temp_TS_idx = tseries_idx_struct.(concentration_TS_idx);
        temp_change_perTV.(concentration_TS_idx).AllData = change_struct.(condition{cond})(temp_TS_idx);
        temp_changeMeanError(1, TV_idx) = nanmean(temp_change_perTV.(concentration_TS_idx).AllData);
        if strcmp(stdOrsem, 'std')
            temp_changeMeanError(2, TV_idx) = nanstd(temp_change_perTV.(concentration_TS_idx).AllData);
        elseif strcmp(stdOrsem, 'sem')
            temp_changeMeanError(2, TV_idx) = nanstd(temp_change_perTV.(concentration_TS_idx).AllData)/...
                (sqrt(length(temp_TS_idx) - 1));
        end
        XVALS =  ones(length(temp_TS_idx), 1) * TV;
        scatter(XVALS, temp_change_perTV.(concentration_TS_idx).AllData, [],...
            temp_color, 'fill', 'MarkerFaceAlpha', .4)
        clear XVALS TV_idx concentration_TS_idx temp_TS_idx
    end
    e(cond) = errorbar(TreatmentLevels, temp_changeMeanError(1, :), temp_changeMeanError(2, :), '-o');
    e(cond).Color = temp_color;
    e(cond).MarkerFaceColor = temp_color;
    e(cond).LineStyle = temp_linestyle;
    e(cond).LineWidth = LW;
    e(cond).Marker = temp_marker;
    %     errorbar_struct(cond) = e;
    clear temp temp_color temp_linestyle temp_marker
end
if strcmp(stdOrsem, 'std')
    %     title({TITLE, 'mean+-std'})
    title([TITLE, 'mean+-std'])
elseif strcmp(stdOrsem, 'sem')
    %     title({TITLE, 'mean+-sem'})
    title([TITLE, 'mean+-sem'])
end
DATALABELS = cellfun(@(x) strrep(x, '_', ' '), condition, 'UniformOutput', 0);
XTICKLABELS = strsplit(num2str(TreatmentLevels));
xlim(XLIMS)
ylim(YLIMS)
hline(0, '--k')
xticks(TreatmentLevels)
xticklabels(XTICKLABELS)
xlabel(XLABEL)
ylabel(YLABEL)
legend(e, DATALABELS, 'Location', 'best')

% If STATS == 1 carry out paired t-test comparing between conditions at each treatment
% value, adjusting for multiple comparisons using Holm-Bonferroni correction
% with a FWER of 0.05
CORRECTED_STATS = struct();
if STATS == 1
    PVALS = nan(1, length(TreatmentLevels));
    for TV = TreatmentLevels
        TV_idx = find(TreatmentLevels == TV);
        COND1 = sprintf('%s%d',condition{1}, TV);
        COND2 = sprintf('%s%d',condition{2}, TV);
        if strcmp(pairedOrunpaired, 'paired')
            [~, PVALS(1,TV_idx)] = ttest(temp_change_perTV.(COND1).AllData, temp_change_perTV.(COND2).AllData);
        elseif strcmp(pairedOrunpaired, 'unpaired')
            [~, PVALS(1,TV_idx)] = ttest2(temp_change_perTV.(COND1).AllData, temp_change_perTV.(COND2).AllData);
        end
    end
    [CORRECTED_STATS.pvals, CORRECTED_STATS.alpha, CORRECTED_STATS.h] = fwer_holmbonf(PVALS, 0.05);
    XTICKLABELS_PVAL = cell(1, length(TreatmentLevels));
    for TV = TreatmentLevels
        TV_idx = find(TreatmentLevels == TV);
        XTICKLABELS_PVAL{1, TV_idx} = sprintf('%d p=%.4f : %d ', CORRECTED_STATS.h(TV_idx), CORRECTED_STATS.pvals(TV_idx), TV);
        clear TV_idx
    end
    xticklabels(XTICKLABELS_PVAL)
    xtickangle(45)
end
end

