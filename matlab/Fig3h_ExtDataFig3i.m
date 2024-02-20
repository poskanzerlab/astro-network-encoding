%% Fig. 3h & Extended Data Fig. 3i: Change in event number / 30s in the 
%  local astrocyte network binned into 50um bands radiating out from the
%  uncaging site (Fig. 3h), including permutation testing to see which time
%  windows have event frequency changes significantly different from
%  chance.
%  Ext. Data Fig. 3i: total event count per band
% Bin events into time bins of a chosen duration, in concentric circles
% radiating out from the uncaging site (with band width of choice),
% including events from chosen regions.
% - Designed for multiple conditions, can include or exclude events from the
%   uncaging cell.
% - Distance of event from landmark is the minimum distance of the event in
%   the event's first frame
% Modified from
% ShollLikeAnalysis_LocAtEvtStart_MultCond_XSecBin_MC20220615.m (Section
% 1,2&3)
%
%  For Fig. 3h & Ext. Data Fig. 3i
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%     Section E: time_window = 30;
%  3. Run the following sections of code:
%
%   A. For each t-series, record the event indices within chosen time windows
%      and chosen area bands excluding events from ramping cells if specified: 
%      area_band = 50; Area_Bands = 25:area_band:275; time_window = 30; WE = 'NeighboringEvents';
%      UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh =
%      'pval_1'; lmk_idx = 1; pulse_num = 10; 
%  B. For chosen event features/parameters, record the events within each
%     time bin and area band per tseries:
%     param_name = {'Number'};
%  C. For each condition, for each area band, aggregate across t-series/recordings
%     and perform a permutation test and calculate p-values to determine which 
%     time windows have post-stim changes signficantly different from chance:
%     param_name = {'Number'}; plotchange = [1]; ExcludeBaseline = 60; PerformPerm = 1; perm_num = 10000;
%  D. Fig. 3h: Semicircle bullseye plots showing change in event frequency
%     in 30s bins post-stim for all area bands of each condition. 
%     condition2plot = {'RuBiGABA', 'RuBiGlu'}; EXCLUDE = {'C', 'C'};
%  E. Ext. Data Fig. 3i: Total number of events per area band for each
%     condition
%     condition2plot = {'RuBiGABA', 'RuBiGlu'}; EXCLUDE = {'C', 'C'};
%
% Michelle Cahill 20240115
%% A. For each t-series, record the event indices within chosen time windows
%     and chosen area bands excluding events from ramping cells if specified 
%     This will be stored in: mydata(tseries).(condition{cond}).(TB).(AB).event_idx
%
%  For Fig. 3h & Ext. Data Fig. 3i:
%   area_band = 50; Area_Bands = 25:area_band:275; time_window = 30; WE = 'NeighboringEvents';
%   UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh =
%   'pval_1'; lmk_idx = 1; pulse_num = 10; 
%    ResNonRes, ResBasedOn and PostStimTime aren't relevant when
%    UseAreaThresh = 0
%
area_band = 50; %distance in um
AB = sprintf('AreaBands%dum', area_band);
Area_Bands = 25:area_band:275; %0:area_band:100; %25:area_band:275 %0:area_band:250; %50; %250: 

time_window = 30; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

WE = 'NeighboringEvents'; %'NeighboringEvents', 'UncagingEvents', 'AllEvents' which events should be included?

UseAreaThresh = 0; %0= take all tseries from that pulse number, 1= take only tseries exceeding a specific threshold
ResNonRes = 1; %0 to plot only non-responders, 1 to plot responders (if UseAreaThresh == 1)
ResBasedOn = 'UncagingEvents' ; % If UseAreaThresh == 1, should responder TS be based on response of uncaging cell or neighboring network? 'NeighboringEvents' or 'UncagingEvents'
PostStimTime = []; %% If UseAreaThresh == 1, time period post-stim for determining responders (should match time_window used in 'Fig2_dFFResponders.m')

ExcludeRampingCells = 1;
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run

condition = fieldnames(mydata);
lmk_idx = 1; %Distance of events will be measured from this landmark drawn in AQuA
pulse_num = 10;

for tseries = 1:length(mydata)
    for cond = 1:length(condition)
        if ~isempty(mydata(tseries).(condition{cond}))
            uncaging_frame = mydata(tseries).(condition{cond}).uncagingframe;
            num_frames = round(time_window / mydata(tseries).(condition{cond}).SecPerFrame); % convert the time window from seconds to frames
            num_bins_pre = ceil((mydata(tseries).(condition{cond}).uncagingframe-1)/num_frames); %How many time bins are needed to cover the entire pre-period
            num_bins_post = ceil((mydata(tseries).(condition{cond}).totalframes - mydata(tseries).(condition{cond}).uncagingframe)/num_frames); %How many time bins are needed to cover the entire post-period
            
            dist2lmk = cellfun(@(x) x(1, lmk_idx), mydata(tseries).(condition{cond}).res.ftsFilter.region.landmarkDist.distPerFrame); %distance of the edge closest to the uncaging pt at the start of each evt
            
            %Collapsing the event idx across regions into a single array
            event_idx_pre = cell(1, num_bins_pre);
            event_idx_post = cell(1, num_bins_post);
            
            if strcmp(WE, 'AllEvents')
                regions = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk, 1);
            elseif strcmp(WE, 'NeighboringEvents')
                regions = find(isnan(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk)); %Identifying the neighboring cells
            elseif strcmp(WE, 'UncagingEvents')
                regions = find(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk == 1); %Identifying the uncaging cell
            end
            
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
            regions = regions(~ismember(regions, ramping_cell_idx));
            
            for t = 1:num_bins_pre
                temp_idx_pre = [];
                temp_idx_post = [];
                for r = regions
                    temp_idx_pre = [temp_idx_pre, mydata(tseries).(condition{cond}).(TB).event_idx_pre{r, t}];
                    temp_idx_post = [temp_idx_post, mydata(tseries).(condition{cond}).(TB).event_idx_post{r, t}];
                end
                clear r
                event_idx_pre{1, t} = temp_idx_pre;
                event_idx_post{1, t} = temp_idx_post;
                clear temp_idx_pre temp_idx_post
            end
            clear t regions
            
            %Sorting events by band
            evt_idx_by_band_pre = cell(length(Area_Bands)-1,num_bins_pre); %rows: 10um bins columns: x-sec bins pre-stim
            evt_idx_by_band_post = cell(length(Area_Bands)-1,num_bins_post); %rows: 10um bins columns: x-sec bins post-stim
            
            for bin = 1:length(Area_Bands)-1
                bin_idx = find(dist2lmk >= Area_Bands(bin) & dist2lmk < Area_Bands(bin + 1)); %evt idx for all events within this distance bin from the uncaging site
                evt_idx_by_band_pre(bin,:) = cellfun(@(x) bin_idx(ismember(bin_idx, x)), event_idx_pre, 'UniformOutput', 0); %evt idx for evts pre-stim with in the bin
                evt_idx_by_band_post(bin,:) = cellfun(@(x) bin_idx(ismember(bin_idx, x)), event_idx_post, 'UniformOutput', 0); %evt idx for evts post-stim with in the bin
                clear bin_idx
            end
            mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).AreaBands = Area_Bands;
            mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).evt_idx_by_band_pre = evt_idx_by_band_pre;
            mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).evt_idx_by_band_post = evt_idx_by_band_post;
            clear bin evt_idx_by_band_pre evt_idx_by_band_post event_idx_pre event_idx_post...
                dist2lmk uncaging_frame num_frames num_bins_pre num_bins_post ramping_cell_idx
        end
    end
end
clear cond tseries
%% B. For chosen event features/parameters, record the events within each time bin and area band per tseries
% For Fig. 3h & Ext. Data Fig. 3i: param_name = {'Number'};

param_name = {'Number'};% Event features to include {'Area', 'Duration', 'Amplitude'};

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            event_idx_pre = mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).evt_idx_by_band_pre;
            event_idx_post = mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).evt_idx_by_band_post;
            
            for param = 1:length(param_name)
                if strcmp(param_name{param}, 'Number')
                    param_pre = 'numevents_Pre'; %Which parameter would you like to aggregate?
                    param_post = 'numevents_Post';
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_pre) = cell2mat(cellfun(@(x) numel(x), event_idx_pre, 'UniformOutput', 0));
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_post) = cell2mat(cellfun(@(x) numel(x), event_idx_post, 'UniformOutput', 0));
                    
                elseif strcmp(param_name{param}, 'Area')
                    param_pre = 'avgarea_Pre'; %Which parameter would you like to aggregate?
                    param_post = 'avgarea_Post';
                    all_pre = extractAfter(param_pre, 'avg');
                    all_post = extractAfter(param_post, 'avg');
                    Area_PerBinPre = cellfun(@(x) mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(x), event_idx_pre, 'UniformOutput', 0);
                    Area_PerBinPost = cellfun(@(x) mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(x), event_idx_post, 'UniformOutput', 0);
                    
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_pre) = cell2mat(cellfun(@(x) mean(x), Area_PerBinPre, 'UniformOutput', 0));
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_post) = cell2mat(cellfun(@(x) mean(x), Area_PerBinPost, 'UniformOutput', 0));
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(all_pre) = Area_PerBinPre;
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(all_post) = Area_PerBinPost;
                    clear Area_PerBinPre Area_PerBinPost all_pre all_post
                    
                elseif strcmp(param_name{param}, 'Duration')
                    param_pre = 'avgduration_Pre'; %Which parameter would you like to aggregate?
                    param_post = 'avgduration_Post';
                    all_pre = extractAfter(param_pre, 'avg');
                    all_post = extractAfter(param_post, 'avg');
                    duration = (mydata(tseries).(condition{cond}).res.ftsFilter.loc.t1 - mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0) * mydata(tseries).(condition{cond}).SecPerFrame;
                    Duration_PerBinPre = cellfun(@(x) duration(x), event_idx_pre, 'UniformOutput', 0);
                    Duration_PerBinPost = cellfun(@(x) duration(x), event_idx_post, 'UniformOutput', 0);
                    
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_pre) = cell2mat(cellfun(@(x) mean(x), Duration_PerBinPre, 'UniformOutput', 0));
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_post) = cell2mat(cellfun(@(x) mean(x), Duration_PerBinPost, 'UniformOutput', 0));
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(all_pre) = Duration_PerBinPre;
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(all_post) = Duration_PerBinPost;
                    clear Duration_PerBinPre Duration_PerBinPost all_pre all_post duration
                    
                elseif strcmp(param_name{param}, 'Amplitude')
                    param_pre = 'avgdffMax_Pre'; %Which parameter would you like to aggregate?
                    param_post = 'avgdffMax_Post';
                    all_pre = extractAfter(param_pre, 'avg');
                    all_post = extractAfter(param_post, 'avg');
                    dffMax_PerBinPre = cellfun(@(x) mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(x), event_idx_pre, 'UniformOutput', 0);
                    dffMax_PerBinPost = cellfun(@(x) mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(x), event_idx_post, 'UniformOutput', 0);
                    
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_pre) = cell2mat(cellfun(@(x) mean(x), dffMax_PerBinPre, 'UniformOutput', 0));
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_post) = cell2mat(cellfun(@(x) mean(x), dffMax_PerBinPost, 'UniformOutput', 0));
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(all_pre) = dffMax_PerBinPre;
                    mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(all_post) = dffMax_PerBinPost;
                    clear dffMax_PerBinPre dffMax_PerBinPost all_pre all_post
                end
                clear param_pre param_post
            end
            clear param event_idx_pre event_idx_post
        end
    end
end

clear tseries cond
%% C. For each condition, for each area band, aggregate across t-series/recordings
% and perform a permutation test and calculate p-values to determine which 
% time windows have post-stim changes signficantly different from chance 
% For the aggregated data: tseries are rows, time bins are columns)
% 
% For Fig. 3h: param_name = {'Number'}; plotchange = [1]; 
%       ExcludeBaseline = 60; PerformPerm = 1; perm_num = 10000;
%
% Calls the function 'MeanFt_XSecStatBins_PermEventOnsetByAreaBand.m'

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\Fig3H';

param_name = {'Number'}; %{'Number', 'Area', 'Duration', 'Amplitude'};
plotchange = [1]; %0 = plot pre- and post-stim bins, 1 = plot post-stim bins change from average baseline (this will be used to determine a 1- or 2-tailed test for p-values and to determine what to plot in the following section [only compatible with plotchange = 1])

ExcludeBaseline = 60; %time ins seconds to exclude from the beginning of the baseline period, for activity to stabilize
BLBinsExclude = ExcludeBaseline / time_window; %how many baseline bins should be excluded?

TS_combined_PrePost = struct();
TS_combined_change = struct();
Responder = struct();

PerformPerm = 1; %0: don't perform any permutation test 1: perform a permutation test
perm_num = 10000; %How many permutations should be carried out for each condition?
if PerformPerm == 1
    XSecBinsPrePost_PermMean = struct(); %Keep track of the permuted means for each X-sec bin pre- and post-stim
    XSecBinsChangeAvgBaseline_PermMean = struct(); %Keep track of the permuted means for the avgerage change from baseline each X-sec post-stim bin
end

% 1. Aggregate together t-series/recordings for each condition and perform a
%    permutation test
for param = 1:length(param_name)
    if strcmp(param_name{param}, 'Number')
        param_pre = 'numevents_Pre'; %Which parameter would you like to aggregate?
        param_post = 'numevents_Post';
        per_change = 'EvtCt_PerChange'; %What are the units of this parameter? put in (), e.g. '(s)'
        param_unit = '';
        if strcmp(WE, 'UncagingEvents')
            YLIMS = [-1 3]; %[-10 10];
        else
            YLIMS = [-1 30]; %[-10 10];
        end
    elseif strcmp(param_name{param}, 'Area')
        param_pre = 'avgarea_Pre'; %Which parameter would you like to aggregate?
        param_post = 'avgarea_Post';
        per_change = 'area_PerChange';
        param_unit = '(um^2)'; %What are the units of this parameter? put in (), e.g. '(s)'
        YLIMS = [-50 200]; %[-10 10];
    elseif strcmp(param_name{param}, 'Duration')
        param_pre = 'avgduration_Pre'; %Which parameter would you like to aggregate?
        param_post = 'avgduration_Post';
        per_change = 'duration_PerChange';
        param_unit = '(s)'; %What are the units of this parameter? put in (), e.g. '(s)'
        YLIMS = [-10 15]; %[-20 40] %[-10 10];
    elseif strcmp(param_name{param}, 'Amplitude')
        param_pre = 'avgdffMax_Pre'; %Which parameter would you like to aggregate?
        param_post = 'avgdffMax_Post';
        per_change = 'dffMax_PerChange';
        param_unit = '(dF/F max)'; %What are the units of this parameter? put in (), e.g. '(s)'
        YLIMS = [-5 5]; %[-20 20];
    end
    
    condition = fieldnames(mydata);
    NumPreStimBins = size(mydata(1).(condition{1}).(TB).(AB)(lmk_idx).(WE).evt_idx_by_band_pre, 2) - BLBinsExclude;
    NumPostStimBins = size(mydata(1).(condition{1}).(TB).(AB)(lmk_idx).(WE).evt_idx_by_band_post, 2);
    
    for cond = 1:length(condition)
        for bin = 1:length(Area_Bands)-1
            ConditionPlusAreaBand = sprintf('%s_%d_%d', condition{cond}, Area_Bands(bin), Area_Bands(bin+1));
            TS_combined_change.(param_name{param}).(ConditionPlusAreaBand) = NaN(length(mydata), NumPostStimBins); %each condition+area band will have a matrix where rows will be tseries, columns will be 30s time bins post- avg(pre)
            TS_combined_PrePost.(param_name{param}).(ConditionPlusAreaBand) = NaN(length(mydata), NumPreStimBins+NumPostStimBins); %each condition+area band will have a matrix where rows will be tseries, columns will be 30s time bins of all included pre-stim time bins followed by all post-stim time bins
        clear ConditionPlusAreaBand
        end
        Responder.(condition{cond}) = []; %keep track of which tseries had at least one post-event above the baseline area threshold
    end
    
    for cond = 1:length(condition)
        for tseries = 1:length(mydata)
            if ~isempty(mydata(tseries).(condition{cond}))
                if UseAreaThresh
                    %             if mydata(tseries).(condition{cond}).Above3stdBLArea.PostEvtsYN && GABA_latencyVSdist(tseries,1) <= TimeThreshold
                    %                 Responder.(condition{cond}) = [Responder.(condition{cond}), tseries];
                    %             end
                    if AboveThresh.(condition{cond})(tseries) == 1
                        Responder.(condition{cond}) = [Responder.(condition{cond}), tseries];
                    end
                end
                
                
                for bin = 1:length(Area_Bands)-1
                    ConditionPlusAreaBand = sprintf('%s_%d_%d', condition{cond}, Area_Bands(bin), Area_Bands(bin+1));
                    temp_meanBaseline = mean(mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_pre)(bin, BLBinsExclude+1:end), 'omitnan');
                    TS_combined_PrePost.(param_name{param}).(ConditionPlusAreaBand)(tseries,:) = [mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_pre)(bin, BLBinsExclude+1:end), mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_post)(bin, :)];
                    TS_combined_change.(param_name{param}).(ConditionPlusAreaBand)(tseries,:) = mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_post)(bin, :) - temp_meanBaseline; 
                    clear ConditionPlusAreaBand temp_meanBaseline
                end
            end
        end
    end
    
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
    
    if PerformPerm == 1 %only include t-series with the chosen pulse number
        condition = fieldnames(mydata);
%         condition = condition(1:2); %Only perform the permutation tests on the first 2 conditions (WT)
        for cond = 1:length(condition)
            for bin = 1:length(Area_Bands)-1
                sprintf('START: %s', datetime)
                ConditionPlusAreaBand = sprintf('%s_%d_%d', condition{cond}, Area_Bands(bin), Area_Bands(bin+1));
                [XSecBinsPrePost_PermMean.(param_name{param}).(ConditionPlusAreaBand),XSecBinsChangeAvgBaseline_PermMean.(param_name{param}).(ConditionPlusAreaBand)] = MeanFt_XSecStatBins_PermEventOnsetByAreaBand(mydata,...
                    condition{cond}, TSERIES.(condition{cond}), time_window, ExcludeBaseline, perm_num, param_name{param}, WE, NumPreStimBins, NumPostStimBins,...
                    AB, lmk_idx, bin);
                sprintf('END: %s', datetime)
                clear ConditionPlusAreaBand
            end
        end
    end
    clear param_pre param_post param_unit YLIMS per_change cp ab
end

% 2. Calculate p-values comparing experimental mean from permuted mean for
%    each condition and area band
if PerformPerm == 1
    for param = 1:length(param_name)
        if strcmp(param_name{param}, 'Number')
            param_unit = '';
        elseif strcmp(param_name{param}, 'Area')
            param_unit = '(um^2)'; %What are the units of this parameter? put in (), e.g. '(s)'
        elseif strcmp(param_name{param}, 'Duration')
            param_unit = '(s)'; %What are the units of this parameter? put in (), e.g. '(s)'
        elseif strcmp(param_name{param}, 'Amplitude')
            param_unit = '(dF/F max)'; %What are the units of this parameter? put in (), e.g. '(s)'
        end
        condition = fieldnames(TS_combined_PrePost.(param_name{param}));
        for cond = 1:length(condition)
            temp_idx = strfind(condition{cond}, '_');
            NT_COND = extractBefore(condition{cond}, temp_idx(2)); %extractBefore(condition{cond}, temp_idx(2)); temp_idx(2) used in data for paper temp_idx(3) used to compare 1 and 10 pulse uncaging responses
            AreaBand = sprintf('%s-%s um from uncaging', extractBefore(extractAfter(condition{cond}, temp_idx(2)), '_'), extractAfter(condition{cond}, temp_idx(3))); % sprintf('%s-%s um from uncaging', extractBefore(extractAfter(condition{cond}, temp_idx(3)), '_'), extractAfter(condition{cond}, temp_idx(4))); %
            for PT = plotchange
                if PT == 0
                    exp_data = TS_combined_PrePost.(param_name{param}).(condition{cond});
                    exp_mean = mean(exp_data(TSERIES.(NT_COND),:), 1, 'omitnan');
                    perm_mean = XSecBinsPrePost_PermMean.(param_name{param}).(condition{cond});
                    
%                     XTICKLABELS = {'-150-120', '-120-90', '-90-60', '-60-30', '-30-0', '0-30', '30-60','60-90','90-120','120-150'};
%                     XTICKLABELS = XTICKLABELS(BLBinsExclude+1:end);
%                     YLABEL = sprintf('mean %s of events %s', lower(param_name{param}), param_unit);
%                     YMIN = -1;
%                     pchange = '';
%                     if exist('pulse_num', 'var')
%                         TITLE = {sprintf('%s %s %s %d-100ms pulses', indicator, strrep(NT_COND, '_', ' '), WE, pulse_num),...
%                             sprintf('event %s per %dS bins: %s', lower(param_name{param}), time_window, AreaBand)};
%                     else
%                         TITLE = {sprintf('%s %s-100ms pulses %s', indicator, strrep(NT_COND, '_', ' '), WE),...
%                             sprintf('event %s per %dS bins: %s', lower(param_name{param}), time_window, AreaBand)};
%                     end
                elseif PT == 1
                    exp_data = TS_combined_change.(param_name{param}).(condition{cond});
                    exp_mean = mean(exp_data(TSERIES.(NT_COND),:), 1, 'omitnan');
                    perm_mean = XSecBinsChangeAvgBaseline_PermMean.(param_name{param}).(condition{cond});
                    
%                     XTICKLABELS = {'0-30', '30-60','60-90','90-120','120-150'};
%                     YLABEL = sprintf('change in mean %s of events %s', lower(param_name{param}), param_unit);
%                     YMIN = -2;
%                     pchange = 'Change';
%                     if exist('pulse_num', 'var')
%                         TITLE = {sprintf('%s %s %s %d-100ms pulses', indicator, strrep(NT_COND, '_', ' '), WE, pulse_num),...
%                             sprintf('event %s per %dS bins: %s', lower(param_name{param}), time_window, AreaBand)};
%                     else
%                         TITLE = {sprintf('%s %s-100ms pulses %s', indicator, strrep(NT_COND, '_', ' '), WE),...
%                             sprintf('event %s per %dS bins: %s', lower(param_name{param}), time_window, AreaBand)};
%                     end
                end
                
                if PT == 0
                    PermGreaterExp = perm_mean >= exp_mean; %if plotting actual value, test specifically for increase from permuted mean (1-tailed)
                elseif PT == 1
                    PermGreaterExp = abs(perm_mean) >= abs(exp_mean); %if plotting change from baseline, test if changes in experimental mean are greater than expected by chance (tests for increase or decrease 2-tailed)
                end
                pval_PerFrame.(condition{cond}) = (sum(PermGreaterExp, 1)+1) ./ (perm_num+1);
                
                % correct for multiple comparisons using Benjamini Yekutieli FDR
                % method with FDR of 5%
                [h.(param_name{param}).(condition{cond}), ~, ~, adj_pval_PerFrame.(condition{cond})] = fdr_bh(pval_PerFrame.(condition{cond}), 0.05, 'dep');
                sig_frames.(condition{cond}) = find(h.(param_name{param}).(condition{cond}));
                
%                 %                 XTICKLABELS_PVAL = cell(1, length(XTICKLABELS));
%                 %                 for BINS = 1:length(XTICKLABELS)
%                 %                     XTICKLABELS_PVAL{1, BINS} = sprintf('%d p=%.4f : %s ', h.(param_name{param}).(condition{cond})(BINS),adj_pval_PerFrame.(condition{cond})(BINS), XTICKLABELS{BINS});
%                 %                 end
%                 %                 clear BINS
%                 XTICKLABELS_PVAL = cell(3, length(XTICKLABELS));
%                 for BINS = 1:length(XTICKLABELS)
%                     XTICKLABELS_PVAL{1, BINS} = sprintf('%s ', XTICKLABELS{BINS});
%                     XTICKLABELS_PVAL{2, BINS} = sprintf('p=%.4f ', pval_PerFrame.(condition{cond})(BINS));
%                     XTICKLABELS_PVAL{3, BINS} = sprintf('adj p=%.4f ', adj_pval_PerFrame.(condition{cond})(BINS));
%                     clear idx
%                 end
%                 clear BINS
%                 tickLabels = MultiLineLabels(XTICKLABELS_PVAL); %allows for multiline xlabels
% 
%                 Structs2Plot = {perm_mean, exp_data};
%                 TS2Plot = [{{1:perm_num}}, {{TSERIES.(NT_COND)}}];
%                 
%                 MAX_VAL = max(max(nanmean(perm_mean, 1)), max(nanmean(exp_mean, 1)));
%                 BUFFER_OrderMag = floor(log(abs(MAX_VAL))./log(10));
%                 BUFFER = 10^BUFFER_OrderMag;
%                 Pval_placement = MAX_VAL + BUFFER;
%                 clear MAX_VAL BUFFER_OrderMag BUFFER
%                 
%                 YLIMS = [YMIN Pval_placement+2];
%                 
%                 DataLabels = {'permutated data', strrep(condition{cond}, '_', '+')};
%                 if ~isempty(sig_frames.(condition{cond}))
%                     DataLabels = [DataLabels, 'p<0.05'];
%                 end
%                 
%                 if strcmp(extractBefore(NT_COND, '_'), 'RuBiGABA') %if it's a GABA condition, the color will be magenta
%                     DataColors = {'-ko', '-mo'};
%                 elseif strcmp(extractBefore(NT_COND, '_'), 'RuBiGlu') %if it's a Glu condition, the color will be cyan
%                     DataColors = {'-ko', '-co'};
%                 end
%                 
%                 %Plot permutated mean in black and experimental mean in magenta or
%                 %cyan
%                 PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt( Structs2Plot, TS2Plot, 1,...
%                     TITLE, YLABEL, YLIMS, sprintf('%dS time bins', time_window), [0 size(Structs2Plot{1,1},2)+1], [1:1:size(Structs2Plot{1,1},2)],...
%                     tickLabels, '', DataColors, ERRORBARS)
% %                 xtickangle(45)
%                 
%                 %plot the bins that are significantly above shuffled data
%                 SF = sig_frames.(condition{cond});
%                 scatter(SF, ones(1, length(sig_frames.(condition{cond})))* Pval_placement,  'k', 'fill')
%                 
%                 legend(DataLabels, 'Location', 'best')
%                 set(gcf,'position', [100, 200, 800, 400])
%                 
%                 cd(save_dir)
%                 if exist('pulse_num', 'var')
%                     saveas(gcf, sprintf('%s%s%s%dPul_%sEvt%s%dSB_Mean%s_Perm%d.tif', indicator, condition{cond}, WE, pulse_num, pchange, param_name{param}, time_window, upper(ERRORBARS), perm_num))
%                     saveas(gcf, sprintf('%s%s%s%dPul_%sEvt%s%dSB_Mean%s_Perm%d.svg', indicator, condition{cond}, WE, pulse_num, pchange, param_name{param}, time_window, upper(ERRORBARS), perm_num))
%                 else
%                     saveas(gcf, sprintf('%s%s%s_%sEvt%s%dSB_Mean%s_Perm%d.tif', indicator, condition{cond}, WE, pchange, param_name{param}, time_window, upper(ERRORBARS), perm_num))
%                     saveas(gcf, sprintf('%s%s%s_%sEvt%s%dSB_Mean%s_Perm%d.svg', indicator, condition{cond}, WE, pchange, param_name{param}, time_window, upper(ERRORBARS), perm_num))
%                 end
%                 close
%                 cd(start_dir)
                
                clear  DataColors DataLabels pchange YLABEL XTICKLABELS TS2Plots Structs2Plot...
                    SF Pval_placement exp_mean perm_mean PermGreaterExp pval_PerFrame adj_pval_PerFrame...
                    sig_frames YLIMS YMIN TITLE TS2Plot XTICKLABELS_PVAL exp_data
            end
            clear temp_idx NT_COND AreaBand
        end
        clear param_unit
    end
end
clear ERRORBARS cond plotchange PT param condition bin tseries
%% D. Fig. 3h: Semicircle bullseye plots for all area bands of each condition on the same plot
%     uses significance values calculated during permutation testing in
%     Section C
%     For Fig. 3h: condition2plot = {'RuBiGABA', 'RuBiGlu'}; 
%                  EXCLUDE = {'C', 'C'}; This will ensure that all area
%                  bands from WT slices will be plotted, all Cx43 and CBX
%                  bands will be excluded
% 
% Calls the function 'SemiCircleBullsEye.m': 
%        which calls createBullseye, cbrewer and natsortfiles
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig3h';
for param = 1:length(param_name)
    condition2plot = {'RuBiGABA', 'RuBiGlu'};%{'RuBiGABA_WT_1_', 'RuBiGlu_WT_1_', 'RuBiGABA_WT_10_', 'RuBiGlu_WT_10_'}; % 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' Which conditions should be plotted
    EXCLUDE = {'C', 'C'}; %Within that dataset, should any condition be left out?
    
    clear condition
    for cp = 1:length(condition2plot)
        if strcmp(condition2plot{cp}, 'AllConditions')
            condition = fieldnames(TS_combined_change.(param_name{param}));
        else
            temp_condition = fieldnames(TS_combined_change.(param_name{param}));
            idx = find(cell2mat(cellfun(@(x) contains(x, condition2plot{cp}), temp_condition, 'UniformOutput', 0)));
            if ~isempty(EXCLUDE{cp})
                exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
                idx = setdiff(idx, exclude_idx);
            end
            condition = temp_condition(idx, 1);
            clear temp_condition idx exclude_idx
        end

        Struct2Plot = struct();
        for cond = 1:length(condition)
            temp_idx = strfind(condition{cond}, '_');
            NT_COND = extractBefore(condition{cond}, temp_idx(2)); %extractBefore(condition{cond}, temp_idx(2)); temp_idx(2) used in data for paper temp_idx(3) used to compare 1 and 10 pulse uncaging responses
            
            temp_TS_idx = TSERIES.(NT_COND);
            Struct2Plot.(condition{cond}) = TS_combined_change.(param_name{param}).(condition{cond})(temp_TS_idx, :);
            %                 if contains(condition{cond}, 'GABA')
            %                     cmap = [cmap; 1 0 1]; %magenta
            %                 elseif contains(condition{cond}, 'Glu')
            %                     cmap = [cmap; 0 1 1]; %cyan
            %                 end
            %                 cmap_num = [cmap_num, cond];
            clear temp_TS_idx temp_idx NT_COND
        end

%         YLABEL = sprintf('change in peak %s (post-pre)', lower(param_name{param}));

%         YLABEL = sprintf('peak post-stim %dS bin', time_window);
        YTICKS = 1:size(Struct2Plot.(condition{1}),2);
        YTICKLABELS = {};
        for YT = YTICKS
            temp_YT_start = (YT-1) * time_window;
            temp_YT_end = (YT) * time_window;
            temp_YT = sprintf('%d-%dS', temp_YT_start, temp_YT_end);
            YTICKLABELS = [YTICKLABELS, temp_YT];
            clear temp_YT_end temp_YT_start temp_YT
        end
        clear YT
        XLABEL = 'distance from uncaging site (um)';
        if isempty(EXCLUDE{cp})
            XTICKLABELS = Area_Bands;
        else
            OB = str2double(extractAfter(EXCLUDE{cp}, '_'));
            XTICKLABELS = setdiff(Area_Bands, OB);
            clear OB
        end
        
        if exist('pulse_num', 'var')
            TITLE = {sprintf('%s %s: %s %d-100ms pulses', indicator, WE, strrep(condition2plot{cp},'_', ' '), pulse_num),...
                sprintf('change in event %s per %dS bins', lower(param_name{param}), time_window)};
        else
            TITLE = {sprintf('%s %s: %s-100ms pulses', indicator, WE, strrep(condition2plot{cp},'_', ' ')),...
                sprintf('change in event %s per %dS bins', lower(param_name{param}), time_window)};
        end
        
        RANGE = 4; %8; 
        SemiCircleBullsEye(Struct2Plot, h.(param_name{param}), RANGE, {'div','RdGy'},...
            XTICKLABELS, YTICKLABELS, param_name{param})
        title(TITLE)
        xlabel(XLABEL)

        cd(save_dir)
        if exist('pulse_num', 'var')
            saveas(gcf, sprintf('%s_%sEx%s_%dPul_ChEvt%s%dSBins_BullsEye_Perm%d.tif', indicator, condition2plot{cp}, EXCLUDE{cp}, pulse_num, lower(param_name{param}), time_window, perm_num))
            saveas(gcf, sprintf('%s_%sEx%s_%dPul_ChEvt%s%dSBins_BullsEye_Perm%d.svg', indicator, condition2plot{cp}, EXCLUDE{cp}, pulse_num, lower(param_name{param}), time_window, perm_num))
        else
            saveas(gcf, sprintf('%s_%sEx%s_ChEvt%s%dSBins_BullsEye_Perm%d.tif', indicator, condition2plot{cp}, EXCLUDE{cp}, lower(param_name{param}), time_window, perm_num))
            saveas(gcf, sprintf('%s_%sEx%s_ChEvt%s%dSBins_BullsEye_Perm%d.svg', indicator, condition2plot{cp}, EXCLUDE{cp}, lower(param_name{param}), time_window, perm_num))
        end
        close
        cd(start_dir)
        clear cond TITLE YLABEL Struct2Plot YTICKS YTICKLABELS...
            RANGE XLABEL XTICKLABELS ARC_LABEL
    end
    clear cp condition2plot PAIRED EXCLUDE PLOT...
        ERRORBARS condition cond
end
clear param ERRORBARS
%% E. Extended Data Fig. 3i: Total number of events in each area band for each condition
%     For Ext. Data Fig. 3i: condition2plot = {'RuBiGABA', 'RuBiGlu'}; 
%                  EXCLUDE = {'C', 'C'}; This will ensure that all area
%                  bands from WT slices will be plotted, all Cx43 and CBX
%                  bands will be excluded
%
% Calls the function violinplot

TOTAL_NUM_PerAreaBand = struct();
F = fieldnames(TS_combined_PrePost.Number);
for f = 1:length(F)
    temp = TS_combined_PrePost.Number.(F{f});
    TOTAL_NUM_PerAreaBand.(F{f}) = sum(temp, 2); %total number of events in that area band per t-series
    clear temp
end
clear f F

save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig3i';
for param = 1:length(param_name)
    condition2plot = {'RuBiGABA', 'RuBiGlu'};%{'RuBiGABA_WT_1_', 'RuBiGlu_WT_1_', 'RuBiGABA_WT_10_', 'RuBiGlu_WT_10_'}; % 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' Which conditions should be plotted
    EXCLUDE = {'C', 'C'}; %Within that dataset, should any condition be left out?
    
    clear condition
    for cp = 1:length(condition2plot)
        if strcmp(condition2plot{cp}, 'AllConditions')
            condition = fieldnames(TS_combined_change.(param_name{param}));
        else
            temp_condition = fieldnames(TS_combined_change.(param_name{param}));
            idx = find(cell2mat(cellfun(@(x) contains(x, condition2plot{cp}), temp_condition, 'UniformOutput', 0)));
            if ~isempty(EXCLUDE{cp})
                exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
                idx = setdiff(idx, exclude_idx);
            end
            condition = temp_condition(idx, 1);
            clear temp_condition idx exclude_idx
        end
        
        Struct2Plot = struct();
        AREABANDS = {};
        for cond = 1:length(condition)
            temp_idx = strfind(condition{cond}, '_');
            NT_COND = extractBefore(condition{cond}, temp_idx(2)); %extractBefore(condition{cond}, temp_idx(3));
            
            temp_TS_idx = TSERIES.(NT_COND);
            AREABANDS = [AREABANDS, extractAfter(condition{cond}, temp_idx(2))]; %[AREABANDS, extractAfter(condition{cond}, temp_idx(3))];
            Struct2Plot.(condition{cond}) = TOTAL_NUM_PerAreaBand.(condition{cond})(temp_TS_idx);
            clear temp_TS_idx temp_idx NT_COND 
        end
        
        DataLabels = cellfun(@(x) strrep(x, '_', '-'), AREABANDS, 'UniformOutput', 0); 
        
        if exist('pulse_num', 'var')
            TITLE = {sprintf('%s %s: %s %d-100ms pulses', indicator, WE, strrep(condition2plot{cp},'_', ' '), pulse_num),...
                sprintf('total number of events %d um area band', area_band)};
        else
            TITLE = {sprintf('%s %s: %s-100ms pulses', indicator, WE, strrep(condition2plot{cp},'_', ' ')),...
                sprintf('total number of events %d um area band', area_band)};
        end

        YLABEL = 'total number of events';
        
        figure();
        v = violinplot(Struct2Plot);
        for cond = 1:length(condition)
            v(cond).ViolinAlpha = 0.25;
        end
        title(TITLE)
        xticklabels(DataLabels)
        xlabel('Distance from uncaging site (um)')
        xtickangle(45)
        ylabel(YLABEL)
        
        cd(save_dir)
        if exist('pulse_num', 'var')
            saveas(gcf, sprintf('%s_%sEx%s_%dPul_TotalEvtNum%dAB_Violin_%s.tif', indicator, condition2plot{cp}, EXCLUDE{cp}, pulse_num, area_band, WE))
            saveas(gcf, sprintf('%s_%sEx%s_%dPul_TotalEvtNum%dAB_Violin_%s.svg', indicator, condition2plot{cp}, EXCLUDE{cp}, pulse_num, area_band, WE))
        else
            saveas(gcf, sprintf('%s_%sEx%s_TotalEvtNum%dAB_Violin_%s.tif', indicator, condition2plot{cp}, EXCLUDE{cp}, area_band, WE))
            saveas(gcf, sprintf('%s_%sEx%s_TotalEvtNum%dAB_Violin_%s.svg', indicator, condition2plot{cp}, EXCLUDE{cp}, area_band, WE))
        end
        close
        cd(start_dir)
        clear cond TITLE YLABEL Struct2Plot YTICKS YTICKLABELS...
            XLABEL XTICKLABELS v AREABANDS DataLabels
    end
    clear cp condition2plot PAIRED EXCLUDE PLOT...
        ERRORBARS condition cond
end
clear param ERRORBARS TOTAL_NUM_PerAreaBand
%%
clear AB area_band Area_Bands bin condition event_idx_post event_idx_pre lmk_idx param param_name...
    param_post param_pre param_unit per_change pulse_num ResNonRes Responder save_dir start_dir TB ThreshDir time_window...
    TS_combined TS_combined_PrePost TS_combined_change tseries UseAreaThresh WE YLIMS condition2plot condition2exclude plotchange...
    STAT stat WTS BLBinsExclude ExcludeBaseline ExcludeRampingCells h NumPostStimBins NumPreStimBins PerformPerm perm_num...
    PostStimTime pulse_num RampingCellThresh ResBasedOn ResNonRes Responder save_dir start_dir TSERIES XSecBinsChangeAvgBaseline_PermMean XSecBinsPrePost_PermMean...
    ans