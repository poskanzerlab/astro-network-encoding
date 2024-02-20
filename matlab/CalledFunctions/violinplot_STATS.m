function violinplot_STATS(DIST_STRUCT, NORMBY, TITLE, YLABEL,...
    cmap, cmap_num, STATS)
% violinplot_STATS plots a violin plot for each field of DIST_STRUCT and
% carries out a one-way ANOVA and multiple comparisons analysis to test
% differences between groups (fields in DIST_STRUCT) or a rank sum test, writing the results
% from the tests in subplots below the violin plot
%
% DIST_STRUCT: a structured array where each field is a group to be plotted
% NORMBY:
% TITLE: string for the  plot title
% YLABEL: string for y-label
% cmap: mx3 array designating different colors
% cmap_num: for each group, which color should be chosen from cmap
% STATS: 'ANOVA' or 'RankSum'
% Michelle Cahill 20210816

% carry out a 1-way ANOVA for unbalanced design and multiple comparisons to
% figure out if any groups are significantly different
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
% if ~isempty(NORMBY)
%     YLABEL = {'percent field active', sprintf('%s', strrep(NORMBY, '_', ''))};
% else
%     YLABEL = '# active px';
% end

figure();
set(gcf,'position', [100, 100, 800, 800])
if strcmp(STATS, 'ANOVA') 
    subplot(4,1,1:2)
elseif strcmp(STATS, 'RankSum')   
    subplot(8,1,1:7)
end
v = violinplot(DIST_STRUCT);
for cond = 1:length(CONDITION)
    if ~isempty(cmap)
        v(cond).ViolinColor = cmap(cmap_num(cond),:);
        v(cond).ScatterPlot.MarkerEdgeColor = cmap(cmap_num(cond),:);
    end
    v(cond).ViolinAlpha = 0.25;
end
title(TITLE)
xticklabels(DATALABELS)
xtickangle(45)
ylabel(YLABEL)
% ylim([0 2.5])

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

