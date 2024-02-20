function PlotMean_IndDataPts_MultCond( DIST_STRUCT, TITLE, YLABEL, cmap, cmap_num,...
    stdOrsem, STATS )
% Michelle Cahill 20211201: modified from violinplot_STATS.m
%   PlotMean_IndDataPts_MultCond plots the mean value of each t-series for each
%   condition & the std or sem of all t-series. A statistical test of the 
%   user's choosing is carried out to test the difference between groups
%
%   DIST_STRUCT: a structured array where each field is a group to be plotted
%   TITLE: title for the plot (as a string)
%   YLABEL: ylabel for the plot (as a string)
%   cmap: mx3 array designating different colors
%   cmap_num: for each group, which color should be chosen from cmap
%   stdOrsem: string 'std' or 'sem' to determine the type of error bars
%   STATS: 'ANOVA' or 'RankSum'

CONDITION = fieldnames(DIST_STRUCT);

dist = [];
C = [];
for cond = 1:length(CONDITION)
    dist = [dist; DIST_STRUCT.(CONDITION{cond})];
    C = [C; repelem([{CONDITION{cond}}], [length(DIST_STRUCT.(CONDITION{cond}))])'];
end

if strcmp(STATS, 'ANOVA')
    [ANOVA_pval, ~, ANOVA_stats] = anova1(dist, C);
    MULTCOMP = multcompare(ANOVA_stats);
elseif strcmp(STATS, 'RankSum')
    ranksum_pval = ranksum(DIST_STRUCT.(CONDITION{1}), DIST_STRUCT.(CONDITION{2}));
end

DATALABELS = cellfun(@(x) strrep(x, '_', ' '), CONDITION, 'UniformOutput', 0);


figure()
set(gcf,'position', [100, 100, 800, 800])
if strcmp(STATS, 'ANOVA') 
    subplot(4,1,1:2)
elseif strcmp(STATS, 'RankSum')   
    subplot(8,1,1:7)
end
hold on
for cond = 1:length(CONDITION)
    XVALS = (ones(length(DIST_STRUCT.(CONDITION{cond})),1) * cond) + 0.05;
    scatter(XVALS, DIST_STRUCT.(CONDITION{cond}), 'Marker', 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmap(cmap_num(cond),:),...
        'MarkerFaceAlpha', 0.4)

    if strcmp(stdOrsem, 'std')
        errorbar(cond-0.05, nanmean(DIST_STRUCT.(CONDITION{cond})), nanstd(DIST_STRUCT.(CONDITION{cond})), '-o',...
            'Color',cmap(cmap_num(cond),:), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmap(cmap_num(cond),:),...
            'MarkerSize', 8, 'LineWidth', 2, 'CapSize', 0);
    elseif strcmp(stdOrsem, 'sem')
        errorbar(cond-0.05, nanmean(DIST_STRUCT.(CONDITION{cond})), nanstd(DIST_STRUCT.(CONDITION{cond}))/(sqrt(length(DIST_STRUCT.(CONDITION{cond}))-1)), '-o',...
            'Color',cmap(cmap_num(cond),:), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmap(cmap_num(cond),:),...
            'MarkerSize', 8, 'LineWidth', 2, 'CapSize', 0);
    end
end

% if strcmp(STATS, 'ANOVA')
% str = sprintf('p = %d', p_val);
% annotation('textbox',[.4 .5 .3 .3],'String',str,'FitBoxToText','on', 'EdgeColor', [1 1 1]);

title([TITLE, sprintf('mean+-%s %s', stdOrsem, STATS)])

ylabel(YLABEL)
xticks(1:length(CONDITION))
xticklabels(DATALABELS)
xtickangle(45)
xlim([0 length(CONDITION)+1])
% breakyaxis([12,46])


if strcmp(STATS, 'ANOVA')
    subplot(4,1,3) %write the results from the 1-way ANOVA
    text(0.01, 0.05, sprintf('1-way ANOVA p-val: %d', ANOVA_pval))
    set(gca, 'Visible', 'off')
    
    h = subplot(4,1,4); %write the results from the multiple comparisons analysis from the ANOVA statistics
    hPos = get(h,'Position');
    T = array2table(MULTCOMP, 'VariableNames', {'G1', 'G2', 'lower_CI',...
        'estimate', 'upper_CI', 'p_val'});
    
    % Get the table in string form.
    TString = evalc('disp(T)');
    % Use TeX Markup for bold formatting and underscores.
    TString = strrep(TString,'<strong>','\bf');
    TString = strrep(TString,'</strong>','\rm');
    TString = strrep(TString,'_','\_');
    % Get a fixed-width font.
    FixedWidth = get(0,'FixedWidthFontName');
    % Output the table using the annotation command.
    annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
        'FontName',FixedWidth,'Units','Normalized','Position',hPos);
    title('tukey-kramer multiple comparison')
elseif strcmp(STATS, 'RankSum')
    subplot(8,1,8) %write the results from the ranksum test
    text(0.01, 0.05, sprintf('Rank sum p-val: %d', ranksum_pval))
    set(gca, 'Visible', 'off')
end

end