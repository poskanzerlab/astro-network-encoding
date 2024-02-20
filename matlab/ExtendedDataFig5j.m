%% Extended Data Fig. 5j: Avgerage pre- and post-stim propagative event frequency 
%  in neighboring astrocytes split by "low" and "high" baseline propagation probability
%  Modified from PropEvtFreq_SplitLowHighBLProp_MC20231030.m
%
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%     Section E: time_window = 30;
%  3. Run the following code:
%     A. Aggregate together eligible cells for each condition and sort on
%        baseline propagation probability for each cell (out of all baseline events in a cell, what proportion of them are propagative): 
%        time_window = 30; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; 
%        ExcludeBaseline = 90; TotalBaseline = 150; PostStimPeriod = 120;
%        pulse_num = 10; UseAreaThresh = 0;
%     B. Plot the average pre- and post-stim propagative event frequency for
%        eligible cells of each condition, split into lower and upper half of
%        baseline propagation probabilities:
%        STAT = 'nonparametric';
%
% Michelle Cahill 20240115
%% A. Aggregate eligible cells (neighboring cells, non-ramping, with some propagative events in the baseline period (60-0s pre-stim))
% One struct with each eligible cell as its own row
%
% For Ext. Data Fig. 5j: time_window = 30; ExcludeRampingCells = 1;
%    RampingCellThresh = 'pval_1'; ExcludeBaseline = 90; TotalBaseline = 150;
%    PostStimPeriod = 120; pulse_num = 10; UseAreaThresh = 0;
%
time_window = 30; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 90; %number of seconds at the beginning of baseline to exclude (must be a multiple of time_window)
TotalBaseline = 150; %number of seconds recorded for baseline
PostStimPeriod = 120; % number of seconds to include post-stim (must be a multiple of time_window)
UseAreaThresh = 0; %0= take all tseries from that pulse number, 1= take only tseries exceeding a specific threshold

pulse_num = 10;

condition = fieldnames(mydata);
condition = condition(1:2); %only carrying out these calculations on WT data

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

for cond = 1:length(condition)
    PropProb_PropEvtRatePrePost.(condition{cond}) = []; %rows will be individual eligible cells, column 1: baseline propagation probability, 2: avg pre-stim (60-0s pre) propagative event rate, 3: post-stim (0-120s) propagative event rate
end

%Aggregate eligible cells (neighboring cells, non-ramping, with some propagative events in the baseline period (60-0s pre-stim))
% One struct with each eligible cell as its own row
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
            PostStimBinsInclude = PostStimPeriod / time_window; % how many post-stim bins to include?
            
            temp_BLPropFreq = mydata(tseries).(condition{cond}).(TB).FreqGrowingProp_Pre; %rows are individual cells, columns are x second time bins pre-stim
            temp_SumBLPropEvt = sum(temp_BLPropFreq(:, BLBinsExclude+1:end), 2); %For each cell, sum the number of baseline propagative events in included baseline bins
            temp_PostStimPropFreq = mydata(tseries).(condition{cond}).(TB).FreqGrowingProp_Post; %rows are individual cells, columns are x second time bins post-stim
            
            temp_BLAllFreq = mydata(tseries).(condition{cond}).(TB).numevents_Pre;
            temp_SumBLAllEvt = sum(temp_BLAllFreq(:, BLBinsExclude+1:end), 2); %For each cell, sum the number of baseline events (all events) in included baseline bins
            
            Cell_idx = find(isnan(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk)); %Identifying the neighboring cells
            Cell_idx = Cell_idx(~ismember(Cell_idx, ramping_cell_idx)); %exclude any ramping cells from the group (if ExcludeRampingCells = 1)
            %             sum_evts = sum(mydata(tseries).(condition{cond}).numevents_PreVPost,2); %Sum the events pre & post in each cell
            Active_cells = Cell_idx(temp_SumBLPropEvt(Cell_idx) > 0); %Find the indices of neighboring cells with at least one baseline propagative event
            %             NoBLProp_cells = find(temp_SumBLPropEvt == 0);
            %             IncludeCellIdx = intersect(Active_cells,NoBLProp_cells);
            
            temp_PropagationProbability = (temp_SumBLPropEvt(Active_cells) ./ temp_SumBLAllEvt(Active_cells)); %For each cell calculate the proportion of baseline events that are propagative
            temp_PropEvtFreqPre = mean(temp_BLPropFreq(Active_cells, BLBinsExclude+1:end), 2, 'omitnan');
            temp_PropEvtFreqPost = mean(temp_PostStimPropFreq(Active_cells, 1:PostStimBinsInclude), 2, 'omitnan');
            
            PropProb_PropEvtRatePrePost.(condition{cond}) = [PropProb_PropEvtRatePrePost.(condition{cond});...
                temp_PropagationProbability, temp_PropEvtFreqPre, temp_PropEvtFreqPost]; %Include each eligible cell as it's own row
            
            clear Cell_idx sum_evts Active_cells ...
                ramping_cell_idx BLBinsExclude temp_BLPropFreq temp_SumBLPropEvt temp_PostStimPropFreq...
                PostStimBinsInclude temp_BLAllFreq temp_SumBLAllEvt temp_PropagationProbability...
                temp_PropEvtFreqPre temp_PropEvtFreqPost
        end
    end
end
clear cond

% Sort eligible cells in each condition by propagation probability of
% baseline events
SortedPropProb_PropEvtRatePrePost = struct();
for cond = 1:length(condition)
    temp_struct = PropProb_PropEvtRatePrePost.(condition{cond});
    [~, sortingIDX] = sort(temp_struct(:,1));
    SortedPropProb_PropEvtRatePrePost.(condition{cond}) = temp_struct(sortingIDX, :);
    clear temp_struct sortingIDX
end

clear cond UseAreaThresh
%% B. Plot the average pre- and post-stim propagative event frequency for
%     eligible cells of each condition, split into lower and upper half of
%     baseline propagation probabilities (SortedPropProb divided in half)
%
% For Ext. Data Fig. 5j: STAT = 'nonparametric';

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig5j';
STAT = 'nonparametric';

for cond = 1:length(condition)
    temp_AllCellEvtRate = SortedPropProb_PropEvtRatePrePost.(condition{cond});
    numTotalCells = size(temp_AllCellEvtRate, 1);
    Cutoff_IDX = ceil(numTotalCells/2);
    LOWER = 1:Cutoff_IDX;
    UPPER = Cutoff_IDX+1:numTotalCells;
    LOWER_UPPER = {LOWER, UPPER};
    for SPLIT = 1:length(LOWER_UPPER) %Plot event rate separately for cells in the lower and upper 50th percentile
        temp_EvtRate = temp_AllCellEvtRate(LOWER_UPPER{SPLIT}, 2:3); %Take the average pre- and post-stim event frequencies for the batch of cells
        CELL_nums = 1:size(temp_EvtRate, 1);
        YLIMS = [0 7];
        YLABEL = sprintf('Propagative event number/%dS', time_window);
        if SPLIT == 1
            CELL_group = 'Lower';
        elseif SPLIT == 2
            CELL_group = 'Upper';
        end
        
        TITLE = {sprintf('%s propagative event number in neighbors with BL prop activity', strrep(condition{cond}, '_', ' ')),...
            sprintf('%s 50th percentile of cells baseline event propagative probability', CELL_group),...
            sprintf('avg/%dS time bins %dS pre %dS post', time_window, TotalBaseline-ExcludeBaseline, PostStimPeriod),...
            sprintf('n = %d cells; exclude ramping cells: %d',size(temp_EvtRate, 1), ExcludeRampingCells)};
        
        PlotTSMean_PreVPost_ChooseStat(temp_EvtRate, CELL_nums, TITLE,...
            YLABEL, YLIMS, 'sem', STAT)
        
        cd(save_dir)
        saveas(gcf, sprintf('%sAvgPropNum%dSNeighWBLProp_%sHalfCellsByBLPropProb_%dSPre_%dSPost_%s.tif', condition{cond}, time_window, CELL_group, TotalBaseline-ExcludeBaseline, PostStimPeriod,  STAT))
        saveas(gcf, sprintf('%sAvgPropNum%dSNeighWBLProp_%sHalfCellsByBLPropProb_%dSPre_%dSPost_%s.svg', condition{cond}, time_window, CELL_group, TotalBaseline-ExcludeBaseline, PostStimPeriod,  STAT))
        close
        cd(start_dir)
        
        clear temp_EvtRate CELL_nums YLIMS YLABEL CELL_group TITLE
    end
    clear SPLIT temp_AllCellEvtRate numTotalCells Cutoff_IDX LOWER UPPER LOWER_UPPER
end

clear cond STAT start_dir save_dir
%%
clear ExcludeBaseline ExcludeRampingCells PostStimPeriod PropProb_PropEvtRatePrePost...
    RampingCellThresh SortedPropProb_PropEvtRatePrePost TB time_window TotalBaseline...
    tseries ans pulse_num TSERIES