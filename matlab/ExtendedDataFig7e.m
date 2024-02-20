%% Extended Data Fig. 7e: compare average baseline event frequency for neighboring cells across multiple recording rounds in the same FOV 
% Modified from CompBLEventFrequency_XSecStatBins_MultCond_MC20230830.m
% sections 1 & 6
%
% For Extended Data Fig. 7e:
%  1. Load 'AQuA_CytoGCaMP_2PMultiRoundGluUncaging_70AU_ExtDataFig7.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B, C & D: as is
%  3. Run ExtendedDataFig7_PreppingDataStruct.m: time_window = 30;
%  4. Run the following sections of code:
%     A. Combine average baseline event frequencies across recordings for
%        neighboring cells in each FOV:
%        time_window = 30; param_name = 'Number'; EVENT_TYPE = {'all'}; EVENT_IDX = {'event_idx_pre'};
%        WhichCells = 2; UseAreaThresh = 0; ExcludeRampingCells = 1;
%        RampingCellThresh = 'pval_1'; pulse_num = 10;
%        ExcludeBaseline = 30; (in this dataset ExcludeBaseline = 30 corresponds to excluding the first 60s of baseline, since
%                the first 30s time bin has been cropped in ExtendedDataFig7_PreppingDataStruct.m)
%    B. Plot individual recordings with lines connecting paired measurements 
%       and a summary statistic (mean+- error bars or median & IQR).
%       Test for significant differences between measurements with repeated
%       measures ANOVA (parametric) or with Friedman Test (non-parametric):
%              ERRORBARS = 'sem'; STAT = 'parametric';

% Michelle Cahill 20240130
%% Section A. Combine average baseline frequencies across t-series/recordings
% Choose which cells to use:
%     -uncaging cells (directly stimulated by uncagin)
%     -neighboring cells (average frequency per FOV)
%     -neighboring cells (individual cells)
% Combine event frequencies of chosen cells into 'TS_combined_AvgBL': with a field for all events, static events and
% propagative events
%
% For Extended Data Fig. 7e:
%   time_window = 30; param_name = 'Number'; EVENT_TYPE = {'all'}; EVENT_IDX = {'event_idx_pre'};
%   WhichCells = 2; PlotIndTS = 0; UseAreaThresh = 0; ExcludeRampingCells = 1;
%   RampingCellThresh = 'pval_1'; pulse_num = 10;
%   ExcludeBaseline = 30; (in this dataset ExcludeBaseline = 30 corresponds to excluding the first 60s of baseline, since
%       the first 30s time bin has been cropped in ExtendedDataFig7_PreppingDataStruct.m)

time_window = 30; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

param_name = 'Number'; % What parameter to plot? 'Number' 
EVENT_TYPE = {'all'}; %Which type of events should be included? Multiple event types can be run in parallel: {'all', 'propagative', 'static'}; 
EVENT_IDX = {'event_idx_pre'}; %These should correspond to EVENT_TYPE: EVENT_TYPE = {'all'}, EVENT_IDX = {'event_idx_pre'}; EVENT_TYPE = {'propagative'}, EVENT_IDX = {'event_idx_pre_prop'}; EVENT_TYPE = {'static'}, EVENT_IDX = {'event_idx_pre_stat'};   {'event_idx_pre','event_idx_pre_prop', 'event_idx_pre_stat'};

WhichCells = 2; %1=uncaging cell, 2=all neighbors averaged for each FOV, 3=each neighboring cell (matched between conditions)

UseAreaThresh = 0; %1 to plot only responding or non-responding cells/t-series, 0 to plot all tseries/recordings in a condition
ResBasedOn = 'UncagingEvents'; %If UseAreaThresh == 1, should responders be based on 'UncagingEvents' or 'NeighboringEvents'
PostStimTime = []; %If UseAreaThresh == 1, time period post-stim for determining responders (should match time_window used in dFF_Responders_MC20220908.m)
ResNonRes = 1; %If UseAreaThresh == 1, 0 to plot only non-responders, 1 to plot responders 

ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 30; %number of seconds at the beginning of baseline to exclude (must be a multiple of time_window)

pulse_num = 10;

condition = fieldnames(mydata);
TSERIES = struct();
for cond = 1:length(condition)
    if exist('pulse_num', 'var')
        if exist('data_sets', 'var') && contains(data_sets{1,2}, 'Cx43')
            pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num)); %For WT, CX43 and CBX dataset
        else
            pulse_cond = strcat(condition{cond}, sprintf('_%d',pulse_num)); %For laser uncaging control, receptor antagonist and multi-round glutamate uncaging dataset
        end
    else
        pulse_cond = extractAfter(condition{cond}, '_');
    end
    if UseAreaThresh
        temptime = sprintf('poststim%ds', PostStimTime);
        RESPONDER_TS = dFFResponders.(temptime);
        clear temptime
        
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
%% Section B. Compare baseline parameters across multiple paired measurements (>2 measurements/columns)
% Plotting individual recordings with lines connecting paired measurements (for
% example, multiple rounds at same FOV) and a summary statistic (mean+- error bars or
% median & IQR)
% Test for significant differences between measurements with repeated
% measures ANOVA (parametric) or with Friedman Test (non-parametric)
%
% Calls the function LinePlot_SummaryStat_IndData.mat
%
% For Extended Data Fig. 7e:
%      ERRORBARS = 'sem'; STAT = 'parametric';

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig7e';

condition = fieldnames(mydata);
ERRORBARS = 'sem';
STAT = 'parametric'; %'parametric': plot mean+-errorbars and run repeated measures ANOVA or 'nonparametric': plot median & IQR and run Friedman Test

if WhichCells == 1
    CELLS = 'uncaging_cells';
elseif WhichCells == 2
    CELLS = 'neighboring_cells_FOV_Avg';
elseif WhichCells == 3
    CELLS = 'neighboring_cells_ind';
end

for ET = 1:length(EVENT_TYPE)
    mat2plot = [];
    for cond = 1:length(condition)
        mat2plot = [mat2plot, TS_combined_AvgBL.(EVENT_TYPE{ET}).(condition{cond})];
    end
    clear cond
    mat2plot_NaNremoved = mat2plot(all(~isnan(mat2plot), 2), :); %remove FOV where any rounds were NaN
    TS = [1:size(mat2plot_NaNremoved, 1)];
    YLABEL = sprintf('event %s /%dS pre-stim', lower(param_name), time_window);
    TITLE = {sprintf('%s: %s', indicator, experiment_type),...
        sprintf('%s event %s per %dS bins: %s', EVENT_TYPE{ET}, lower(param_name), time_window, strrep(CELLS, '_', ' '))};
    
    if strcmp(STAT, 'parametric')
        SUMMARY = 'mean';
    elseif strcmp(STAT, 'nonparametric')
        SUMMARY = 'median';
    end
    
    LinePlot_SummaryStat_IndData(mat2plot_NaNremoved, TS, TITLE, condition, [], YLABEL, SUMMARY, ERRORBARS, 'r')
    
    cd(save_dir)
    saveas(gcf, sprintf('%s_%sEvt%s%dSBins_CompAvgPre_%s_%s_%s.tif', indicator, EVENT_TYPE{ET}, lower(param_name), time_window, CELLS, SUMMARY, STAT))
    saveas(gcf, sprintf('%s_%sEvt%s%dSBins_CompAvgPre_%s_%s_%s.svg', indicator, EVENT_TYPE{ET}, lower(param_name), time_window, CELLS, SUMMARY, STAT))
    close
    cd(start_dir)
    clear mat2plot mat2plot_NaNremoved TS TITLE YLABEL SUMMARY cond 
end
clear ET CELLS ERRORBARS condition start_dir save_dir STAT
%%
clear TSERIES EVENT_IDX EVENT_TYPE ExcludeBaseline ExcludeRampingCells param_name...
    param_post param_pre param_unit PlotIndTS PostStimTime RampingCellThresh ResBasedOn...
    ResNonRes TB time_window tseries UseAreaThresh WhichCells YLIMS pulse_num TS_combined_AvgBL