function BarPlotGrouped_PreVPostTSMeanERROR_IndTrials(STRUCT, CONDITION, COLORS, TITLE, YLABEL, XLABEL, ERRORBARS)
%BarPlotGrouped_PreVPostTSMeanERROR_IndTrials plots a grouped bar chart with a bar for pre-stim
% mean and post-stim mean for each condition
% Individual trials will be plotted over the bars +- errorbars connecting corresponding
% trials between conditions for each group
% A 2-sample t-test will be run between pre and post values for each condition. P-values
% will be corrected for multiple comparisons using  a Bonferroni-Holm correction
% This function calls fwer_holmbonf (MATLAB>Add-ons>Multiple Testing Toolbox> fwer_holmbonf)
% Michelle Cahill 20211012
%   STRUCT: a structured array with each field containing a matrix of values
%           rows are individual trials, columns are various measurments for
%           that trial. The mean will be taken across all rows for the the
%           bar plots. ex. FIELD1 = [a b]; a = number of t-series, b =
%           principal components. The average will be taken across t-series
%           and the average for each principal component for each
%           designated condition will be plotted.
%           The values from conditions of corresponding PCs will be
%           compared with a two-sample t-test and corrected for multiple
%           comparisons using Bonferroni-Holm correction
%   CONDITION: the fieldnames in STRUCT to compare (designed for 2
%              conditions) {'CONDITION1', 'CONDITION2'}
%   COLORS: [2 3] matrix with rgb color specifications
%   TITLE, YLABEL, XLABEL: strings for the title, ylabel and xlabel


MEAN_PER_PC = nan(length(CONDITION), 2); %rows will be conditions, columns will be pre-stim and post-stim values
ERRORBARS_PER_PC = nan(length(CONDITION), 2); %rows will be conditions, columns will be pre-stim and post-stim values
for cond = 1:length(CONDITION)
    MEAN_PER_PC(cond,:) = nanmean(STRUCT.(CONDITION{cond}), 1);
    if strcmp(ERRORBARS, 'std')
        ERRORBARS_PER_PC(cond,:) = std(STRUCT.(CONDITION{cond}), 0, 1, 'omitnan');
    elseif strcmp(ERRORBARS, 'sem')
        temp_STD = std(STRUCT.(CONDITION{cond}), 0, 1, 'omitnan');
        N = size(STRUCT.(CONDITION{cond}), 1);
        ERRORBARS_PER_PC(cond,:) = temp_STD / sqrt(N - 1);
        clear temp_STD N
    end
end

%To properly place errorbars on bar chart
%Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(MEAN_PER_PC);

%Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

PVal_PER_GROUP = nan(length(CONDITION),1);
for cond = 1:length(CONDITION)
    [~, PVal_PER_GROUP(cond,1)] = ttest(STRUCT.(CONDITION{cond})(:,1), STRUCT.(CONDITION{cond})(:,2));
end

% Correcting p-values with Bonferroni-Holm correction
% (MATLAB>Add-ons>Multiple Testing Toolbox> fwer_holmbonf)
[CORRECTED_STATS.pvals, ~, CORRECTED_STATS.h] = fwer_holmbonf(PVal_PER_GROUP, 0.05);

XTICKLABELS_PVAL = cell(1, length(CONDITION));
for cond = 1:length(CONDITION)
    XTICKLABELS_PVAL{1, cond} = sprintf('%d p=%.4f : %s ', CORRECTED_STATS.h(cond), CORRECTED_STATS.pvals(cond), CONDITION{cond});
end

%         DATALABELS = cellfun(@(x) strrep(x, '_', ' '), CONDITION, 'UniformOutput', 0);
DATALABELS = {'pre-stim', 'post-stim'};
MaxVals = [];
for cond = 1:length(CONDITION)
    MaxVals = [MaxVals, max(STRUCT.(CONDITION{cond}), [], 'all')];
end
MaxY = max(MaxVals);

figure()
b = bar(MEAN_PER_PC, 'hist');
b(1).FaceColor = COLORS(1, :);
b(2).FaceColor = COLORS(2, :);
hold on

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, MEAN_PER_PC(:,i), ERRORBARS_PER_PC(:,i), 'k', 'linestyle', 'none');
end

COND = 1;
for cond = 1:length(CONDITION)
    X_COND1 = COND - 0.15;
    X_COND2 = COND + 0.15;
    Y_VAL_PAIRS = [STRUCT.(CONDITION{cond})(:,1), STRUCT.(CONDITION{cond})(:,2)];
    keepTS = sum(~isnan(Y_VAL_PAIRS),2) == 2;
    plot([X_COND1, X_COND2], Y_VAL_PAIRS, 'Color', [0.25 0.25 0.25])
    %             plot([X_COND1, X_COND2], Y_VAL_PAIRS(keepTS, :), 'Color', [0.6 0.6 0.6])
    COND = COND + 1;
    clear X_COND1 X_COND2 Y_VAL_PAIRS keepTS
end
ylabel(YLABEL)
xlabel(XLABEL)
xticklabels(XTICKLABELS_PVAL)
xtickangle(45)
legend(DATALABELS, 'Location', 'Best')
title(TITLE)
ylim([0 MaxY + (MaxY*.2)])
end