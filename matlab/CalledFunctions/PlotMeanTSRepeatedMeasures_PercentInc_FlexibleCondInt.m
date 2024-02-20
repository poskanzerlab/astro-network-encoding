function PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt( PerIn_ByTS_preVpost, pulse_idx, pulse_num, TITLE, YLABEL, YLIMS, XLABEL, XLIMS, XTICKS,XTICKLABELS,...
    DATALABELS, DATACOLORS, stdOrsem, cmap, cmap_num)
% Michelle Cahill 20201006, updated 20201019
%PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt creates a line plot showing mean +-
%std or sem for different conditions across different experimental
%intensities (pulse number, concentration etc)
%
%   PerIn_ByTS_preVpost: 1xn cell array, where each cell contains and mxt
%   numerical array for the percent increase of a certain parameter from pre
%   to post stim. n = the number of conditions to plot (ie GABA vs Glu)
%   m = the number of tseries for that condition (ie all GABA TS)
%   t = x-second time bins post uncaging
%   pulse_idx: 1xn cell array with a cell for each condition containing a
%   1xp cell array, where each cell contains the indices of the
%   tseries included in a particular experimental intensity (ie pulse
%   number)p = the number of different experimental intensities.
%   pulse_num: 1xp matrix listing intensity values (ie [1 5 10] for pulses)
%   TITLE: title for the plot (as a string)
%   YLABEL: ylabel for the plot (as a string)
%   YLIMS: ylims for the plot [yMin, yMax]
%   XLABEL: xlabel for the plot (as a string)
%   XLIMS: xlims for the plot [xMin, xMax] [0, t + 1]
%   XTICKS: [1:length(pulse_num)]
%   XTICKLABELS: strings to label each tickmark {'1', '5', '10'}
%   DATALABELS: strings to label each condition {'RuBiGABA', RuBiGlu'}
%   DATACOLORS: specification for how each condition will be plotted
%   {'-mo', '-co'}
%   stdOrsem: 'std' or 'sem' for the errorbars
%   cmap: if DATACOLORS is emtpy, specify the color for each condition
%   using a matrix of rgb values
%   cmap_num: if DATACOLORS is empty, specify which row of cmap to pull
%   each condition color from

figure()
hold on

for cond = 1:length(PerIn_ByTS_preVpost)
    % For each pulse intensity record mean, std and SEM for tseries for
    % each time bin
    % pulse_idx = {pulse_1_idx, pulse_5_idx, pulse_10_idx};
    % pulse_num = [1 5 10];
    pulse_idx_for_cond = pulse_idx{cond};
    for pulse = 1:length(pulse_num)
        PerIn_AllTimeBins = nan(3,size(PerIn_ByTS_preVpost{1,cond},2)); %columns are time bins row1:mean row2:std row3:SEM
        PerIn = PerIn_ByTS_preVpost{1,cond};
        
        PerIn_Finite = PerIn(pulse_idx_for_cond{pulse},:);
        PerIn_Finite = PerIn_Finite(isfinite(PerIn_Finite(:,1)), :); %discard t-series with infinite values in the first time bin
        PerIn_AllTimeBins(1, :) = nanmean(PerIn_Finite,1);
        PerIn_AllTimeBins(2, :) = nanstd(PerIn_Finite,1);
        PerIn_AllTimeBins(3, :) = nanstd(PerIn_Finite,1)/(sqrt(length(pulse_idx_for_cond{pulse}) - 1));
        %         PerIn_AllPulses(1, pulse) = nanmean(PerIn(pulse_idx{pulse}));
        %         PerIn_AllPulses(2, pulse) = nanstd(PerIn(pulse_idx{pulse}));
        %         PerIn_AllPulses(3, pulse) = nanstd(PerIn(pulse_idx{pulse}))/(sqrt(length(pulse_idx{pulse}) - 1));
        
        
        if strcmp(stdOrsem, 'std')
            %             e = errorbar(PerIn_AllTimeBins(1, :),  PerIn_AllTimeBins(2, :), DATACOLORS{cond*pulse});
            e = errorbar(XTICKS, PerIn_AllTimeBins(1, :),  PerIn_AllTimeBins(2, :), '-o');
        elseif strcmp(stdOrsem, 'sem')
            %             e = errorbar(PerIn_AllTimeBins(1, :),  PerIn_AllTimeBins(3, :), DATACOLORS{cond*pulse});
            e = errorbar(XTICKS, PerIn_AllTimeBins(1, :),  PerIn_AllTimeBins(3, :), '-o');
        end
        
        if ~isempty(DATACOLORS)
            temp = DATACOLORS{cond*pulse};
            %         temp_color = temp(2);
            %         temp_linestyle = temp(1);
            %         temp_marker = temp(3);
            temp_color = extractBefore(extractAfter(temp, length(temp)-2), 2);
            temp_linestyle = extractBefore(temp, length(temp)-1);
            temp_marker = extractAfter(temp, length(temp)-1);
            
            e.Color = temp_color;
            e.MarkerFaceColor = temp_color;
            e.MarkerEdgeColor = 'none';
            e.LineWidth = 1.5;
            e.LineStyle = temp_linestyle;
            e.Marker = temp_marker;
        else
            e.Color = cmap(cmap_num(cond*pulse),:);
            e.MarkerFaceColor = cmap(cmap_num(cond*pulse),:);
            e.MarkerEdgeColor = 'none';
            e.LineWidth = 1.5;
        end
    end
end

if strcmp(stdOrsem, 'std')
    %     title({TITLE, 'mean+-std'})
    title([TITLE, 'mean+-std'])
elseif strcmp(stdOrsem, 'sem')
    %     title({TITLE, 'mean+-sem'})
    title([TITLE, 'mean+-sem'])
end
xlim(XLIMS)
ylim(YLIMS)
hline(0, '--k')
xticks(XTICKS)
if ~isempty(XTICKLABELS)
    xticklabels(XTICKLABELS)
end
xlabel(XLABEL)
ylabel(YLABEL)
legend(DATALABELS, 'Location', 'best')

