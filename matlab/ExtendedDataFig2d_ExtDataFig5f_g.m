% Measuring baseline event frequencies in chosen cells/regions:
% Calculate baseline event frequency per cell/region using the res file for each tseries
% Modified from CompBLEventFrequency_XSecStatBins_MultCond_MC20230830.m
% Only written for param_name 'Number'
%
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%     Section E: time_window = 30;
%  3. Run the following sections of code:
%
%  For Extended Data Fig. 2d: baseline event frequency (all events) prior to GABA and
%                             glutamate uncaging in directly stimulated astrocytes
%      A. Combine average baseline event frequency across t-series/recordings 
%         of the same condition:
%         time_window = 30; param_name = 'Number'; EVENT_TYPE = {'all'}; EVENT_IDX = {'event_idx_pre'}; WhichCells = 1;
%         UseAreaThresh = 0;ExcludeRampingCells = 1; RampingCellThresh = 'pval_1';
%         ExcludeBaseline = 60; pulse_num = 10;
%      B. Scatter plots of baseline event frequency for two chosen
%         conditions with median and IQR, comparing the two conditions with a
%         paired test: 
%         condition = {'RuBiGABA_WT', 'RuBiGlu_WT'}; STAT = 'nonparametric';
%
%  Extended Data Fig. 5g: baseline event frequency (separated into static
%             and propagative events) prior to GABA and glutamate uncaging in
%             the same neighboring cells
%      A. Combine average baseline event frequency across t-series/recordings 
%         of the same condition:
%         time_window = 30; param_name = 'Number'; 
%         EVENT_TYPE = {'propagative', 'static'}; EVENT_IDX = {'event_idx_pre_prop', 'event_idx_pre_stat'};
%         WhichCells = 2; UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1';
%         ExcludeBaseline = 60; pulse_num = 10;
%      B. Scatter plots of baseline event frequency for two chosen
%         conditions with median and IQR, comparing the two conditions with a
%         paired test: 
%         condition = {'RuBiGABA_WT', 'RuBiGlu_WT'}; STAT = 'nonparametric';
%
%  Extended Data Fig. 5f: baseline event frequency (separated into static
%             and propagative events) prior to GABA and glutamate uncaging.
%             Comparing WT and Cx43 slices (not paired recordings)
%      A. Combine average baseline event frequency across t-series/recordings 
%         of the same condition:
%         time_window = 30; param_name = 'Number'; 
%         EVENT_TYPE = {'propagative', 'static'}; EVENT_IDX = {'event_idx_pre_prop', 'event_idx_pre_stat'};
%         WhichCells = 2; UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1';
%         ExcludeBaseline = 60; pulse_num = 10;
%      C. Violin plots comparing baseline event frequency between WT and
%         CX43 networks:
%         condition2plot = {'RuBiGABA', 'RuBiGlu'}; EXCLUDE = {'CBX', 'CBX'};
%         COLORS = {'RdPu', 'Blues'}; COLOR_Scheme = {'seq', 'seq'}; STATS = {'RankSum', 'RankSum'};
%
% Michelle Cahill 20240115
%% A. Combine average baseline values across t-series/recordings of the same condition
% Choose to combine uncaging cell values, average neighboring cell values
% or individual neighboring cell values into 
% 'TS_combined_AvgBL': with a field for all events, static events and
%                      propagative events
%
% For Ext. Data Fig. 2d: time_window = 30; param_name = 'Number';
%        EVENT_TYPE = {'all'}; EVENT_IDX = {'event_idx_pre'}; WhichCells = 1;
%        UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1';
%        ExcludeBaseline = 60; pulse_num = 10;
%
% For Ext. Data Fig. 5f-g: time_window = 30; param_name = 'Number'; 
%       EVENT_TYPE = {'propagative', 'static'}; EVENT_IDX = {'event_idx_pre_prop', 'event_idx_pre_stat'};
%       WhichCells = 2; UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1';
%       ExcludeBaseline = 60; pulse_num = 10;

time_window = 30; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

param_name = 'Number'; % What parameter to plot? 'Number' (Not tested for the following:'Area', 'Duration' or 'Amplitude')
EVENT_TYPE = {'propagative', 'static'}; %{'all', 'propagative', 'static'};
EVENT_IDX = {'event_idx_pre_prop', 'event_idx_pre_stat'}; %{'event_idx_pre','event_idx_pre_prop', 'event_idx_pre_stat'};

WhichCells = 2; %1=uncaging cell, 2=all neighbors combined, 3=each neighboring cell (matched between conditions)

UseAreaThresh = 0; %0= take all tseries from that pulse number, 1= take only tseries exceeding a specific threshold
ResBasedOn = 'UncagingEvents'; %'UncagingEvents', 'NeighboringEvents'
PostStimTime = []; %[]; %time period post-stim for determining responders (should match time_window used in dFF_Responders_MC20220908.m)
ResNonRes = 1; %0 to plot only non-responders, 1 to plot responders (if UseAreaThresh == 1)

ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 60; %number of seconds at the beginning of baseline to exclude (must be a multiple of time_window)

pulse_num = 10;

condition = fieldnames(mydata);
TSERIES = struct();
for cond = 1:length(condition)
    if exist('pulse_num', 'var')
        if contains(data_sets{1,2}, 'Cx43')
            pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num)); %For WT, CX43 and CBX dataset
        else
            pulse_cond = strcat(condition{cond}, sprintf('_%d',pulse_num)); %For laser uncaging control and receptor antagonist dataset
        end
    else
        pulse_cond = extractAfter(condition{cond}, '_');
    end
    if UseAreaThresh
%         if isempty(PostStimTime)
%             RESPONDER_TS = dFFResponders.poststim;
%         else
            temptime = sprintf('poststim%ds', PostStimTime);
            RESPONDER_TS = dFFResponders.(temptime);
            clear temptime
%         end
        
        if ResNonRes == 1
            TSERIES.(condition{cond}) = intersect(RESPONDER_TS.(ResBasedOn).(condition{cond}), pulse_x_idx.(pulse_cond)); %return the t-series of responders within that pulse condition
        elseif ResNonRes == 0
            TSERIES.(condition{cond}) = setdiff(pulse_x_idx.(pulse_cond), RESPONDER_TS.(ResBasedOn).(condition{cond})); %return the t-series of non-responders within that pulse condition
        end
    else
        TSERIES.(condition{cond}) = intersect(pulse_x_idx.(pulse_cond), 1:length(mydata));
    end
    clear pulse_cond RESPONDER_TS
end
clear cond temp_idx

for ET = 1:length(EVENT_TYPE)
    for cond = 1:length(condition)
        TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}) = []; %rows will be tseries, column 1 will be avg BL
    end
end

for ET = 1:length(EVENT_TYPE)
    for cond = 1:length(condition)
        for tseries = 1:length(mydata)
            if ~isempty(mydata(tseries).(condition{cond}))
                if ExcludeRampingCells
                    if isfield(mydata(tseries).(condition{cond}), 'RampingCellIdx')
                        if isempty(RampingCellThresh)
                            ramping_cell_idx = mydata(tseries).(condition{cond}).RampingCellIdx;
                        else
                            ramping_cell_idx = mydata(tseries).(condition{cond}).RampingCellIdx.(RampingCellThresh);
                        end
                    else
                        ramping_cell_idx = [];
                    end
                else
                    ramping_cell_idx = [];
                end
                
                BLBinsExclude = ExcludeBaseline / time_window; %how many baseline bins should be excluded?

                temp_EvtIdx = mydata(tseries).(condition{cond}).(TB).(EVENT_IDX{ET});
                temp_EvtNum = cell2mat(cellfun(@(x) numel(x), temp_EvtIdx, 'UniformOutput', 0));
                
                if WhichCells == 1 %To include only the uncaging cells
                    Cell_idx = find(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk == 1);
                    if ~ismember(Cell_idx, ramping_cell_idx) && ismember(tseries, TSERIES.(condition{cond})) %if the uncaging cell is not a ramping cell and it comes from a tseries included in TSERIES
                        TS_val = nanmean(temp_EvtNum(Cell_idx, BLBinsExclude+1:end));
                        TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}) = [TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}); TS_val];
                    else %if the uncaging cell is a ramping cell or from a tseries that shouldn't be included, place nans in the row for that t-series
                        TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}) = [TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}); NaN];
                    end
                    clear Cell_idx TS_val ramping_cell_idx BLBinsExclude
                elseif WhichCells == 2 %To include neighboring cells with activity (average per t-series)
                    if ismember(tseries, TSERIES.(condition{cond})) %if this a tseries that should be included
                        Cell_idx = find(isnan(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk)); %Identifying the neighboring cells
                        Cell_idx = Cell_idx(~ismember(Cell_idx, ramping_cell_idx)); %exclude any ramping cells from the group (if ExcludeRampingCells = 1)
                        sum_evts = sum(mydata(tseries).(condition{cond}).numevents_PreVPost,2); %Sum the events pre & post in each cell
                        Nonzero_cells = Cell_idx(sum_evts(Cell_idx) > 0); %Modified 3/21/2022 Find the indices of neighboring cells with at least one event
                        TS_val = nanmean(temp_EvtNum(Nonzero_cells, BLBinsExclude+1:end),'all');
                        TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}) = [TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}); TS_val];
                    else %if this a tseries that should be excluded
                        TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}) = [TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}); NaN];
                    end
                    clear Cell_idx TS_val sum_evts Nonzero_cells ramping_cell_idx TS_val BLBinsExclude
                elseif WhichCells == 3 %To include neighboring cells (each as a separate row), including NaNs for individual neighboring cells that ramp/with no activity
                    AllCell_AvgBL = nanmean(temp_EvtNum(:, BLBinsExclude+1:end), 2); %Record the average baseline of all cells in FOV
                    if ismember(tseries, TSERIES.(condition{cond})) %if this a tseries that should be included
                        UncagingCell_idx = find(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk == 1);
                        sum_evts = sum(mydata(tseries).(condition{cond}).numevents_PreVPost,2); %Sum the events pre & post in each cell
                        ZeroCells = find(sum_evts == 0);
                        AllCell_AvgBL(UncagingCell_idx) = NaN; %put a nan in for the uncaging cell
                        AllCell_AvgBL(ramping_cell_idx) = NaN;  %put a nan in for any ramping cells
                        AllCell_AvgBL(ZeroCells) = NaN;  %put a nan in for any cells with no detected AQuA-events
                        TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}) = [TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}); AllCell_AvgBL];
                    else %if this a tseries that should be excluded
                        AllCell_AvgBL(:) = NaN;  %put a nan in for any cells with no detected AQuA-events
                        TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}) = [TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}); AllCell_AvgBL];
                    end
                    clear AllCell_AvgBL UncagingCell_idx sum_evts ZeroCells ramping_cell_idx BLBinsExclude
                end
                clear temp_EvtNum temp_EvtIdx
            else
                TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}) = [TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond}); NaN];
            end
        end
    end
end

clear ET cond
%% B. Compare baseline parameters across conditions with paired measurements(currently designed for 2 conditions)
%  Combine pre-stim data and plot baseline values for two chosen conditions as 
%  scatter plot with median and IQR, comparing the two conditions using 
%  t-test (STAT = 'parametric') or wilcoxon sign rank (STAT = 'nonparametric')
%  Uses WhichCells and EVENT_TYPE specified in Section A. 
% 
% Calls the function: 'Scatter_Median.m'
%
% For Ext. Data Fig. 2d and Ext. Data Fig. 5g: 
% condition = {'RuBiGABA_WT', 'RuBiGlu_WT'}; STAT = 'nonparametric';

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig5g';
condition = {'RuBiGABA_WT', 'RuBiGlu_WT'};
STAT = 'nonparametric'; %'parametric' or 'nonparametric'

if WhichCells == 1
    CELLS = 'uncaging_cells';
elseif WhichCells == 2
    CELLS = 'neighboring_cells_FOV_Avg';
elseif WhichCells == 3
    CELLS = 'neighboring_cells_ind';
end

for ET = 1:length(EVENT_TYPE)
    temp_comb = [];
    for cond = 1:length(condition)
        temp_comb = [temp_comb, TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond})];
    end
    
    struct2plot = temp_comb;
    TS = 1:size(temp_comb, 1);
    
    if strcmp(STAT, 'parametric')
        [~, p_val] = ttest(struct2plot(TS,1), struct2plot(TS,2)); %paired t-test
    elseif strcmp(STAT, 'nonparametric')
        p_val = signrank(struct2plot(TS,1), struct2plot(TS,2)); %Wilcoxon signed rank test
    end
    
    TITLE = {sprintf('%s: %s', indicator, experiment_type),...
        sprintf('%s event %s per %dS bins: %s', EVENT_TYPE{ET}, lower(param_name), time_window, strrep(CELLS, '_', ' ')),...
        sprintf('avg baseline: %s p-val: %.15f', STAT, p_val)};
    
    %plot pre and post-stim values as scatter plots with median and
    %interquartile range
    DataLabels = cellfun(@(x) sprintf('%s Pre', x), condition, 'UniformOutput', 0);
    XLABEL = (sprintf('%s %s', DataLabels{1}, param_name));
    YLABEL = (sprintf('%s %s', DataLabels{2}, param_name));
    Scatter_Median(struct2plot, TS, TITLE, XLABEL, YLABEL)
%     xlim([0 20])
%     ylim([0 20])
    
    cd(save_dir)
    saveas(gcf, sprintf('%s_%sEvt%s%dSBins_ScatMedianIQRCompAvgPre_%s.tif', indicator, EVENT_TYPE{ET}, lower(param_name), time_window, CELLS))
    saveas(gcf, sprintf('%s_%sEvt%s%dSBins_ScatMedianIQRCompAvgPre_%s.svg', indicator, EVENT_TYPE{ET}, lower(param_name), time_window, CELLS))
    close
    cd(start_dir)
    
    clear TITLE YLABEL  XLABEL struct2plot TS temp_comb cond DataLabels p_val
end

clear ET CELLS ERRORBARS condition start_dir save_dir STAT
%% C. Extended Data Fig. 5f: Compare baseline parameters across multiple unpaired conditions
%    Show baseline frequencies as a violin plot for each chosen condition
%  
% Calls the function violinplot_STATS.m and cbrewer.m
%
% For Ext. Data Fig. 5f: condition2plot = {'RuBiGABA', 'RuBiGlu'}; EXCLUDE = {'CBX', 'CBX'};
%   COLORS = {'RdPu', 'Blues'}; COLOR_Scheme = {'seq', 'seq'}; STATS = {'RankSum', 'RankSum'};

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig5f';

condition2plot = {'RuBiGABA', 'RuBiGlu'}; % 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' Which conditions should be plotted
EXCLUDE = {'CBX', 'CBX'}; %Within that dataset, should any condition be left out?
COLORS = {'RdPu', 'Blues'}; %Colors from cbrewer: 'RdPu', 'Blues', 'YlGnBu' %'RdPu', 'Blues', 'MgCy', 'MgCy', 'MgCy'
COLOR_Scheme = {'seq', 'seq'}; %Color schemes from cbrewer
STATS = {'RankSum', 'RankSum'}; %'ANOVA' or 'RankSum' for each comparison

if WhichCells == 1
    CELLS = 'uncaging_cells';
elseif WhichCells == 2
    CELLS = 'neighboring_cells_FOV_Avg';
elseif WhichCells == 3
    CELLS = 'neighboring_cells_ind';
end

clear condition
for cp = 1:length(condition2plot)
    if strcmp(condition2plot{cp}, 'AllConditions')
        temp_condition = fieldnames(mydata);
        if ~isempty(EXCLUDE{cp})
            idx = 1:length(temp_condition);
            exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
            idx = setdiff(idx, exclude_idx);
            condition = temp_condition(idx, 1);
        else
            condition = temp_condition;
        end  
        clear temp_condition idx exclude_idx
%         condition = fieldnames(mydata);
    else
        temp_condition = fieldnames(mydata);
        if ~contains(condition2plot{cp}, '_')
            idx = find(cell2mat(cellfun(@(x) contains(x, condition2plot{cp}), temp_condition, 'UniformOutput', 0)));
        else %to pull out conditions with different strings (for example: All RuBiGABA conditions and NoRuBi_LaserUncagingControl would be used with 'RuBiGABA_NoRuBi'
            idx_1 = cell2mat(cellfun(@(x) contains(x, extractBefore(condition2plot{cp}, '_')), temp_condition, 'UniformOutput', 0));
            idx_2 = cell2mat(cellfun(@(x) contains(x, extractAfter(condition2plot{cp}, '_')), temp_condition, 'UniformOutput', 0));
            idx = find(idx_1 | idx_2);
        end
        if ~isempty(EXCLUDE{cp})
            exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
            idx = setdiff(idx, exclude_idx);
        end
        condition = temp_condition(idx, 1);
        clear temp_condition idx exclude_idx idx_1 idx_2
    end
    
%     pulse_num = pulses{cp};
    
    if ~isempty(COLOR_Scheme{cp})
        cmap = cbrewer(COLOR_Scheme{cp}, COLORS{cp}, 8,'PCHIP');
        if contains(condition2plot{cp}, '_')
            cmap(1,:) = [0.5 0.5 0.5];
        end
        if contains(condition2plot{cp}, '_') || strcmp(COLOR_Scheme{cp}, 'div')
            %             cmap_num = flip(floor(linspace(1, 8, length(condition)*length(pulse_num))));
            cmap_num = flip(floor(linspace(1, 8, length(condition))));
        else
            %             cmap_num = flip(floor(linspace(3, 8, length(condition)*length(pulse_num))));
            cmap_num = flip(floor(linspace(3, 8, length(condition))));
        end
    else
        cmap = [1 0 1; 0 1 1]; %magneta and cyan
        cmap_num = [1, 2];
    end
    
    for ET = 1:length(EVENT_TYPE)
        struct2plot = struct();
        for cond = 1:length(condition)
            struct2plot.(condition{cond}) = TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond});
        end
        clear cond
        YLABEL = sprintf('event %s /%dS pre-stim', lower(param_name), time_window);
        TITLE = {sprintf('%s: %s', indicator, experiment_type),...
            sprintf('%s event %s per %dS bins: %s', EVENT_TYPE{ET}, lower(param_name), time_window, strrep(CELLS, '_', ' ')),...
            'avg baseline'};
        
        violinplot_STATS(struct2plot, [], TITLE, YLABEL,...
            cmap, cmap_num, STATS{cp})
        
        if contains(EVENT_TYPE{ET}, 'propagative')
            if strcmp(STATS, 'ANOVA')
                subplot(4,1,1:2)
                ylim([0 2.5])
            elseif strcmp(STATS, 'RankSum')
                subplot(8,1,1:7)
                ylim([0 2.5])
            end
        elseif contains(EVENT_TYPE{ET}, 'static')
            if strcmp(STATS, 'ANOVA')
                subplot(4,1,1:2)
                ylim([0 20])
            elseif strcmp(STATS, 'RankSum')
                subplot(8,1,1:7)
                ylim([0 20])
            end
        end
        
        cd(save_dir)
        saveas(gcf, sprintf('%s_Ex%s_%sEvt%s%dSBins_VioAvgPre_%s_%s.tif', condition2plot{cp}, EXCLUDE{cp}, EVENT_TYPE{ET}, lower(param_name), time_window, CELLS, STATS{cp}))
        saveas(gcf, sprintf('%s_Ex%s_%sEvt%s%dSBins_VioAvgPre_%s_%s.svg', condition2plot{cp}, EXCLUDE{cp}, EVENT_TYPE{ET}, lower(param_name), time_window, CELLS, STATS{cp}))
        close
        cd(start_dir)

        clear struct2plot cond TITLE YLABEL 
    end
    clear ET cmap cmap_num condition 
end

clear condition2plot EXCLUDE COLORS COLOR_Scheme STATS cp CELLS start_dir save_dir
%%
clear TSERIES EVENT_IDX EVENT_TYPE ExcludeBaseline ExcludeRampingCells param_name...
    param_post param_pre param_unit PlotIndTS PostStimTime RampingCellThresh ResBasedOn...
    ResNonRes TB time_window tseries UseAreaThresh WhichCells YLIMS pulse_num TS_combined_AvgBL