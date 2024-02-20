%% Fig. 2l-o & Extended Data Fig. 2f: Event number/30s near and far from the uncaging site
%  pre- and post-uncaging.
% Bin events into time bins of a chosen duration, in concentric circles
% radiating out from the uncaging site (with band width of choice),
% including events from chosen regions.
% - Designed for multiple conditions, can include or exclude events from the
%   uncaging cell.
% - Distance of event from landmark is the minimum distance of the event in
%   the event's first frame
% Modified from
% ShollLikeAnalysis_LocAtEvtStart_MultCond_XSecBin_MC20220615.m (Section
% 1,2&4)
%
% For Fig. 2l-o & Ext. Data Fig. 2f:
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%     Section E: time_window = 30;
%  3. Run 'Fig2_dFFResponders.m': WE = 'UncagingEvents'; ExcludeBaseline = 60; time_window = [];
%  4. Run the following sections of code:
%   A. For each t-series, record the event indices within chosen time windows
%      and chosen area bands excluding events from ramping cells if specified: 
%      area_band = 10; Area_Bands = 0:area_band:100; time_window = 30;
%      WE = 'UncagingEvents'; UseAreaThresh = 1; ResNonRes = 1; ResBasedOn = 'UncagingEvents'; PostStimTime = [];
%      ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; lmk_idx = 1;
%  B. For chosen event features/parameters, record the events within each
%     time bin and area band per tseries:
%     param_name = {'Number'};
%  C. For each condition, combine events across designated area bands for
%     measurements "near" vs "far" from uncaging, combine across t-series & plot
%    1. Aggregate across t-series and combine event features across area bands
%       for averages near vs far from the uncaging site
%    2. Plot value of average event feature near vs far from uncaging
%       in time bins of chosen duration (change from baseline) of a single condition on the same
%       graph (Fig. 2l & n)
%    3. Plot value of average event feature either near or far from uncaging
%       average pre-stim vs specific post-stim x sec time bin per t-series
%       (mean +- sem) for each condition (Currently for 4th post-stim time
%        bin) (Fig. 2m & o, Ext. Data Fig. 2f):
%       param_name = {'Number'}; pulse_num = 10; plotchange = 1; GroupBy= {'tseries'};
%       STAT = 'nonparametric'; AreaCutOff = 10; ExcludeBaseline = 60; Post_TimeBin_Idx = 4;
%
% Michelle Cahill 20240115.m
%% A. For each t-series, record the event indices within chosen time windows
%     and chosen area bands excluding events from ramping cells if specified 
%     This will be stored in: mydata(tseries).(condition{cond}).(TB).(AB).event_idx
%
%  For Fig. 2l-o & Ext. Data Fig. 2f:
%  area_band = 10; Area_Bands = 0:area_band:100; time_window = 30;
%  WE = 'UncagingEvents'; UseAreaThresh = 1; ResNonRes = 1; ResBasedOn = 'UncagingEvents'; PostStimTime = [];
%  ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; lmk_idx = 1;

area_band = 10; %distance in um
AB = sprintf('AreaBands%dum', area_band);
Area_Bands = 0:area_band:100; %distance of the first area band from the uncaging site and distance of the outer most area band from the uncaging site

time_window = 30; %time window (in seconds) for static, nonoverlapping time bins taken before and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

WE = 'UncagingEvents'; %'NeighboringEvents', 'UncagingEvents', 'AllEvents' which events should be included?

UseAreaThresh = 1; %0= take all tseries from that pulse number, 1= take only tseries exceeding a specific threshold
ResNonRes = 1; %0 to plot only non-responders, 1 to plot responders (if UseAreaThresh == 1)
ResBasedOn = 'UncagingEvents' ; % If UseAreaThresh == 1, should responder TS be based on response of uncaging cell or neighboring network? 'NeighboringEvents' or 'UncagingEvents'
PostStimTime = []; %% If UseAreaThresh == 1, time period post-stim for determining responders (should match time_window used in 'Fig2_dFFResponders.m')

ExcludeRampingCells = 1;
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run

condition = fieldnames(mydata);
lmk_idx = 1; %Distance of events will be measured from this landmark drawn in AQuA

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
% For Fig. 2l-o & Ext. Data Fig. 2f: param_name = {'Number'};

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

%% C. For each condition, combine events across designated area bands for
%  measurements "near" vs "far" from uncaging, combine across t-series & plot
% 1. Aggregate across t-series and combine event features across area bands
%    for averages near vs far from the uncaging site
% 2. Plot value of average event feature near vs far from uncaging
%    in time bins of chosen duration (change from baseline) of a single condition on the same
%    graph (Fig. 2l & n)
% 3. Plot value of average event feature either near or far from uncaging
%    average pre-stim vs specific post-stim x sec time bin per t-series
%    (mean +- sem) for each condition (Currently for 4th post-stim time
%    bin) (Fig. 2m & o, Ext. Data Fig. 2f)
%
% For Fig. 2l-o & Ext. Data Fig. 2f:
% param_name = {'Number'}; pulse_num = 10; plotchange = 1; GroupBy= {'tseries'};
% STAT = 'nonparametric'; AreaCutOff = 10; ExcludeBaseline = 60; Post_TimeBin_Idx = 4;
%
% Calls the following functions:
% 'PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt.m'
% 'PlotTSMean_PreVPost_ChooseStat.m'
% 'violinplot.m'
% 'fdr_bh.m'

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig2l_o_ExtDataFig2f';

param_name = {'Number'};
pulse_num = 10;
plotchange = 1; %0 to plot all x-sec time bins pre and post, 1 to plot post- avg(pre)
GroupBy= {'tseries'};%{'tseries', 'slice', 'mouse'};
STAT = 'nonparametric'; %'parametric' to run a t-test to compare pre- and post-stim or 'nonparametric' to run a Wilcoxon signed rank test
if strcmp(STAT, 'parametric')
    stat = 'p';
elseif strcmp(STAT, 'nonparametric')
    stat = 'np';
end

if UseAreaThresh == 1
    if ResNonRes == 1
        WTS = 'Res';
    else
        WTS = 'NRes';
    end
else
    WTS = 'AllTS';
end

AreaCutOff = 10; %bands <= AreaCutoff will be combined to quantify events near the uncaging site and > AreaCutOff will be combined for events far from the uncaging site

ExcludeBaseline = 60; %time ins seconds to exclude from the beginning of the baseline period, for activity to stabilize
BLBinsExclude = ExcludeBaseline / time_window; %how many baseline bins should be excluded?

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

for param = 1:length(param_name)
    if strcmp(param_name{param}, 'Number')
        param_pre = 'numevents_Pre'; %Which parameter would you like to aggregate?
        param_post = 'numevents_Post';
        param_unit = '';
        if strcmp(WE, 'UncagingEvents')
            if plotchange == 1
                YLIMS = [-1 5]; %[-1 3] [-1 5]
            else
                YLIMS = [-10 10]; %[-10 10]; [-1 3]
            end
            YLIMSAVG = [0 15];
        else
            YLIMS = [-1 30]; %[-10 10];
            YLIMSAVG = [0 30];
        end
    elseif strcmp(param_name{param}, 'Area')
        param_pre = 'area_Pre'; %Which parameter would you like to aggregate?
        param_post = 'area_Post';
        param_unit = '(um^2)'; %What are the units of this parameter? put in (), e.g. '(s)'
        if plotchange == 1
            YLIMS = [-100 200];
        else
            YLIMS = [-50 200]; %[-10 10]; [-1 3]
        end
        YLIMSAVG = [0 500];
    elseif strcmp(param_name{param}, 'Duration')
        param_pre = 'duration_Pre'; %Which parameter would you like to aggregate?
        param_post = 'duration_Post';
        %         per_change = 'duration_PerChange';
        param_unit = '(s)'; %What are the units of this parameter? put in (), e.g. '(s)'
        if plotchange == 1
            YLIMS = [-10 30];
        else
            YLIMS = [-20 40]; %[-10 10]; [-1 3]
        end
        YLIMSAVG = [0 30];
    elseif strcmp(param_name{param}, 'Amplitude')
        param_pre = 'dffMax_Pre'; %Which parameter would you like to aggregate?
        param_post = 'dffMax_Post';
        param_unit = '(dF/F max)'; %What are the units of this parameter? put in (), e.g. '(s)'
        if plotchange == 1
            YLIMS = [-5 5];
        else
            YLIMS = [-20 20]; %[-10 10]; [-1 3]
        end
        YLIMSAVG = [0 15];
    end
    
    % 1. Aggregate across t-series and combine event features across area bands
    % for averages near vs far from the uncaging site
    FT_XSEC_BINS = struct(); %each field will have the average value for evt features combined across area bands for near vs far events in X sec static bins (per t-series)
    FT_XSEC_BINS_CHANGE = struct(); %each field will have the average value for evt features combined across area bands for near vs far events in post-stim X sec static bins (post - avg pre)
    FT_XSEC_PREPOSTAVG = struct(); %each field will have the average value for evt features combined across area bands for near vs far events average of X sec static bins pre and post-stim (per t-series)
    Responder = struct();
    for cond = 1:length(condition)
        totalXSecBins = size(mydata(1).(condition{cond}).(TB).event_idx_pre, 2) + size(mydata(1).(condition{cond}).(TB).event_idx_post, 2) - BLBinsExclude;
        FT_XSEC_BINS.(condition{cond}).near = nan(length(mydata), totalXSecBins); %rows are t-series, columns are x second time bins (pre followed by post)
        FT_XSEC_BINS.(condition{cond}).far = nan(length(mydata), totalXSecBins); %rows are t-series, columns are x second time bins (pre followed by post)
        FT_XSEC_BINS_CHANGE.(condition{cond}).near = nan(length(mydata),size(mydata(1).(condition{cond}).(TB).event_idx_post, 2)); %rows are t-series, columns are post-stim x second time bins (change from average baseline)
        FT_XSEC_BINS_CHANGE.(condition{cond}).far = nan(length(mydata),size(mydata(1).(condition{cond}).(TB).event_idx_post, 2)); %rows are t-series, columns are post-stim x second time bins (change from average baseline)
        FT_XSEC_PREPOSTAVG.(condition{cond}).near = nan(length(mydata),2); %rows are t-series, column 1: avg pre-stim column 2: avg post-stim
        FT_XSEC_PREPOSTAVG.(condition{cond}).far = nan(length(mydata),2); %rows are t-series, column 1: avg pre-stim column 2: avg post-stim
        
        Responder.(condition{cond}) = []; %keep track of which tseries had at least one post-event above the baseline area threshold
        
        for tseries = TSERIES.(condition{cond}) %1:length(mydata)
            if ~isempty(mydata(tseries).(condition{cond}))
                for bin = 1:2 %bin 1 will be near the uncaging site, bin 2 will be distant from the uncaging site
                    if bin == 1 %bin near the uncaging site
                        Bands2Combine = find(Area_Bands(2:end) <= AreaCutOff); %due to binning of events into area bands in Section A, the near bin includes events < AreaCutoff 
                        DIST = 'near';
                    elseif bin == 2 %bin distant from the uncaging site
                        Bands2Combine = find(Area_Bands(2:end) > AreaCutOff); %due to binning of events into area bands in Section A, the near bin includes events >= AreaCutoff 
                        DIST = 'far';
                    end
                    
                    ft_pre = mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_pre)(Bands2Combine, BLBinsExclude+1:end);
                    ft_post = mydata(tseries).(condition{cond}).(TB).(AB)(lmk_idx).(WE).(param_post)(Bands2Combine, :);
                    if strcmp(param_name{param}, 'Number')
                        ft_pre_comb = sum(ft_pre, 1); %add the number of events across area bands for corresponding time bins
                        ft_post_comb = sum(ft_post, 1); %add the number of events across area bands for corresponding time bins
                    else
                        ft_pre_comb = cell2mat(arrayfun(@(x) nanmean(horzcat(ft_pre{:,x})), 1:size(ft_pre,2), 'UniformOutput', 0)); %combine evt features across area bands for corresponding time bins and average
                        ft_post_comb = cell2mat(arrayfun(@(x) nanmean(horzcat(ft_post{:,x})), 1:size(ft_post,2), 'UniformOutput', 0));
                    end
                    
                    FT_XSEC_BINS.(condition{cond}).(DIST)(tseries, :) = [ft_pre_comb, ft_post_comb]; %all pre-stim x sec bins, followed by all post-stim x sec bins
                    FT_XSEC_BINS_CHANGE.(condition{cond}).(DIST)(tseries, :) = ft_post_comb - nanmean(ft_pre_comb); %post-stim x sec bins (change from average baseline)
                    FT_XSEC_PREPOSTAVG.(condition{cond}).(DIST)(tseries, :) = [nanmean(ft_pre_comb), nanmean(ft_post_comb)]; %column1: avgerage of all pre-stim x sec bins column2: avgerage of all post-stim x sec bins
                    
                    clear ft_pre ft_post ft_pre_comb ft_post_comb DIST Bands2Combine
                end
                clear bin
            end
        end
        clear tseries totalXSecBins
    end
    
    % 2. Plot value of average event feature near vs far from uncaging
    %    in x sec static bins of a single condition on the same graph
    %    (Fig. 2l & n)
    for g = 1:length(GroupBy)
        for cond = 1:length(condition)
            if plotchange == 1
                if strcmp(GroupBy{g}, 'tseries')
                    struct2plot = {FT_XSEC_BINS_CHANGE.(condition{cond}).near, FT_XSEC_BINS_CHANGE.(condition{cond}).far} ;
                    PLOT_TS = TSERIES.(condition{cond});
                end
                XTICKLABELS = {'0-30', '30-60','60-90','90-120','120-150'};
                YLABEL = sprintf('change in mean %s of events %s', lower(param_name{param}), param_unit);
                pchange = 'Change';
            elseif plotchange == 0
                if strcmp(GroupBy{g}, 'tseries')
                    struct2plot = {FT_XSEC_BINS.(condition{cond}).near, FT_XSEC_BINS.(condition{cond}).far} ;
                    PLOT_TS = TSERIES.(condition{cond});
                elseif strcmp(GroupBy{g}, 'slice')
                    struct2plot = {FT_XSEC_BINS_BySlice.(condition{cond}).near, FT_XSEC_BINS_BySlice.(condition{cond}).far} ;
                    PLOT_TS = Slice_IDX.(condition{cond});
                elseif strcmp(GroupBy{g}, 'mouse')
                    struct2plot = {FT_XSEC_BIN_ByMouse.(condition{cond}).near, FT_XSEC_BINS_ByMouse.(condition{cond}).far} ;
                    PLOT_TS = Mouse_IDX.(condition{cond});
                end
                XTICKLABELS = {'-150-120', '-120-90', '-90-60', '-60-30', '-30-0', '0-30', '30-60','60-90','90-120','120-150'};
                XTICKLABELS = XTICKLABELS(BLBinsExclude+1:end);
                YLABEL = sprintf('mean %s of events %s', lower(param_name{param}), param_unit);
                pchange = '';
            end
            
            DATALABELS = {sprintf('events < %dum from uncaging', AreaCutOff), sprintf('events >= %dum from uncaging', AreaCutOff) };
            
            
            PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt( struct2plot,...
                repelem({{PLOT_TS}},size(struct2plot,2)), 1,...
                {sprintf('%s-100ms pulses: CytoGCaMP %s', extractAfter(condition{cond}, '_'), extractBefore(condition{cond}, '_')), sprintf('event %s per %dS bins compared 30S avg baseline', lower(param_name{param}), time_window),...
                sprintf('near vs far from uncaging: avg by %s', GroupBy{g}), sprintf('cell-independent: %s %s', WE, WTS)}, YLABEL, YLIMS, sprintf('%dS time bins', time_window),...
                [0 size(struct2plot{1,1},2)+1], [1:1:size(struct2plot{1,1},2)],XTICKLABELS,...
                DATALABELS, {'-ro','-ko'}, 'sem')
                xtickangle(45)
                
                cd(save_dir)
                saveas(gcf, sprintf('Cyto_CellInd%s_%s2P%dumNearVFar_%sEvt%s%dSbinsAvg%s_%s_MeanSEM.tif', WE, (condition{cond}), AreaCutOff, pchange, param_name{param}, time_window, GroupBy{g}, WTS))
                saveas(gcf, sprintf('Cyto_CellInd%s_%s2P%dumNearVFar_%sEvt%s%dSbinsAvg%s_%s_MeanSEM.svg', WE, (condition{cond}), AreaCutOff, pchange, param_name{param}, time_window, GroupBy{g}, WTS))
                close
                cd(start_dir)
                
            clear DATALABELS struct2plot XTICKLABELS XTICKLABELS YLABEL pchange ThreshDir TS2Plot PLOT_TS
        end
    end
    clear g
    
    % 3. Plot value of average event feature either near or far from uncaging
    %    average pre-stim vs specific post-stim x sec time bin per t-series
    %    (mean +- sem) for each condition (Ext. Data Fig. 2f) & 
    %    change from average baseline as violin plots: near and far
    %    together (Fig. 2m & o)
    for g = 1:length(GroupBy)
        for cond = 1:length(condition)
            Post_TimeBin_Idx = 4; %Chose the post-stim time bin to plot
            totalXSecBins_pre = size(mydata(1).(condition{cond}).(TB).event_idx_pre, 2) - BLBinsExclude;
            SPEC_TB = totalXSecBins_pre + Post_TimeBin_Idx;
            for ab = 1:2
                if ab == 1
                    DIST = 'near';
                    YLIMSAVGMod = YLIMSAVG;
                elseif ab == 2
                    DIST = 'far';
                    if strcmp(param_name{param}, 'Number')
                        YLIMSAVGMod = YLIMSAVG * 3;
                    else
                        YLIMSAVGMod = YLIMSAVG;
                    end
                end
                
                if strcmp(GroupBy{g}, 'tseries')
                    Param_Pre = FT_XSEC_PREPOSTAVG.(condition{cond}).(DIST)(:,1);
                    Param_Post = FT_XSEC_BINS.(condition{cond}).(DIST)(:, SPEC_TB);
                    PLOT_TS = TSERIES.(condition{cond});
                elseif strcmp(GroupBy{g}, 'slice')
                    Param_Pre = FT_XSEC_PREPOSTAVG_BySlice.(condition{cond}).(DIST)(:,1);
                    Param_Post = FT_XSEC_BINS_BySlice.(condition{cond}).(DIST)(:, SPEC_TB);
                    PLOT_TS = Slice_IDX.(condition{cond});
                elseif strcmp(GroupBy{g}, 'mouse')
                    Param_Pre = FT_XSEC_PREPOSTAVG_ByMouse.(condition{cond}).(DIST)(:,1);
                    Param_Post = FT_XSEC_BINS_ByMouse.(condition{cond}).(DIST)(:, SPEC_TB);
                    PLOT_TS = Mouse_IDX.(condition{cond});
                end
                
                Param_PreVPost = [Param_Pre, Param_Post];
                YLABEL = sprintf('event %s %s', lower(param_name{param}), param_unit);
                POST_TIMEBIN_START = (time_window * Post_TimeBin_Idx) - time_window;
                POST_TIMEBIN_END = time_window * Post_TimeBin_Idx;
                
                
                PlotTSMean_PreVPost_ChooseStat( Param_PreVPost, PLOT_TS,...
                    {sprintf('Event %s by %s: %s %s %s uncaging site %s', param_name{param}, GroupBy{g}, strrep(condition{cond}, '_', ' '), WE, DIST, WTS),...
                    sprintf('avg %ds time bins pre vs %d-%dS post', time_window, POST_TIMEBIN_START, POST_TIMEBIN_END),...
                    sprintf('exclude ramping cells: %d, excluding baseline: %dS',ExcludeRampingCells, ExcludeBaseline)},...
                    YLABEL, YLIMSAVGMod, 'sem', STAT )
                
                cd(save_dir)
                saveas(gcf, sprintf('Cyto_%s2P%s%dum%sUncSite_Evt%s_AVG%s%dSTBPre_%d_%dSPost_%s_%s.tif', condition{cond}, WE, AreaCutOff, DIST, upper(param_name{param}), GroupBy{g}, time_window, POST_TIMEBIN_START, POST_TIMEBIN_END, WTS, stat))
                saveas(gcf, sprintf('Cyto_%s2P%s%dum%sUncSite_Evt%s_AVG%s%dSTBPre_%d_%dSPost_%s_%s.svg', condition{cond}, WE, AreaCutOff, DIST, upper(param_name{param}), GroupBy{g}, time_window, POST_TIMEBIN_START, POST_TIMEBIN_END, WTS, stat))
                close
                cd(start_dir)
                
                clear Param_PreVPost YLABEL DIST YLIMSAVGMod Param_Pre Param_Post PLOT_TS
            end
            clear ab totalXSecBins_pre SPEC_TB
            
            % violin plot: change from baseline near and far for each
            % t-series
            Struct2Plot = struct();
            if strcmp(GroupBy{g}, 'tseries')
                PLOT_TS = TSERIES.(condition{cond});
                Struct2Plot.near = FT_XSEC_BINS_CHANGE.(condition{cond}).near(PLOT_TS, Post_TimeBin_Idx);
                Struct2Plot.far = FT_XSEC_BINS_CHANGE.(condition{cond}).far(PLOT_TS, Post_TimeBin_Idx);
            elseif strcmp(GroupBy{g}, 'slice')
                PLOT_TS = Slice_IDX.(condition{cond});
                Struct2Plot.near = FT_XSEC_BINS_CHANGE_BySlice.(condition{cond}).near(PLOT_TS, Post_TimeBin_Idx);
                Struct2Plot.far = FT_XSEC_BINS_CHANGE_BySlice.(condition{cond}).far(PLOT_TS, Post_TimeBin_Idx);
            elseif strcmp(GroupBy{g}, 'mouse')
                PLOT_TS = Mouse_IDX.(condition{cond});
                Struct2Plot.near = FT_XSEC_BINS_CHANGE_ByMouse.(condition{cond}).near(PLOT_TS, Post_TimeBin_Idx);
                Struct2Plot.far = FT_XSEC_BINS_CHANGE_ByMouse.(condition{cond}).far(PLOT_TS, Post_TimeBin_Idx);
            end
            
            PVALS = NaN(1,2);
            if strcmp(STAT, 'parametric')
                [~, PVALS(1,1)] = ttest(Struct2Plot.near);
                [~, PVALS(1,2)] = ttest(Struct2Plot.far);
            elseif strcmp(STAT, 'nonparametric')
                PVALS(1,1) = signrank(Struct2Plot.near);
                PVALS(1,2) = signrank(Struct2Plot.far);
            end
            % correct for multiple comparisons using Benjamini Yekutieli FDR
            % method with FDR of 5%
            [~, ~, ~, adj_pval] = fdr_bh(PVALS, 0.05, 'dep');
            
            DIST = fieldnames(Struct2Plot);
            DataLabels = cellfun(@(x) strrep(x, '_', ' '), DIST, 'UniformOutput', 0);
            DataLabels_PVAL = cell(3, length(DataLabels));
            for f = 1:length(DataLabels)
                %             DataLabels_PVAL{1, cond} = sprintf('p=%.4f : %s ',adj_pval(cond), DataLabels{cond});
                DataLabels_PVAL{1, f} = sprintf('%s ',DataLabels{f});
                DataLabels_PVAL{2, f} = sprintf('p=%.4f ',PVALS(f));
                DataLabels_PVAL{3, f} = sprintf('adj p=%.4f ',adj_pval(f));
                clear idx
            end
            tickLabels = MultiLineLabels(DataLabels_PVAL); %allows for multiline xlabels
            clear f

            TITLE = {sprintf('Event %s by %s: %s-100ms pul %s %s', param_name{param}, GroupBy{g}, strrep(condition{cond}, '_', ' '), WE, WTS),...
                sprintf('avg %ds time bins pre vs %d-%dS post', time_window, POST_TIMEBIN_START, POST_TIMEBIN_END),...
                sprintf('exclude ramping cells: %d, excluding baseline: %dS',ExcludeRampingCells, ExcludeBaseline)};
            YLABEL = sprintf('change in event %s (%s) (post peak - avg pre)', lower(param_name{param}), param_unit);
            
            figure();
            v = violinplot(Struct2Plot);
            title(TITLE)
            xticklabels(tickLabels)
            ylabel(YLABEL)
            hline(0, '--k')
            ylim([-12 12])
            
            cd(save_dir)
            saveas(gcf, sprintf('%s_%s2P%s%dumNearFarViolin_Evt%s_AVGCh%s%dSTBPre_%d_%dSPost_%s_%s.tif', indicator, condition{cond}, WE, AreaCutOff, upper(param_name{param}), GroupBy{g}, time_window, POST_TIMEBIN_START, POST_TIMEBIN_END, WTS, stat))
            saveas(gcf, sprintf('%s_%s2P%s%dumNearFarViolin_Evt%s_AVGCh%s%dSTBPre_%d_%dSPost_%s_%s.svg', indicator, condition{cond}, WE, AreaCutOff, upper(param_name{param}), GroupBy{g}, time_window, POST_TIMEBIN_START, POST_TIMEBIN_END, WTS, stat))
            close
            cd(start_dir)
            
            clear Post_TimeBin_Idx POST_TIMEBIN_START POST_TIMEBIN_END TITLE YLABEL DataLabels_PVAL DataLabels v PVALS adj_pval DIST PLOT_TS Struct2Plot tickLabels
        end
    end
    clear cond g param_pre param_post param_unit YLIMS YLIMSAVG per_change cp ab FT_XSEC_BINS FT_XSEC_BINS_CHANGE FT_XSEC_PREPOSTAVG...
        pchange struct2plot Responder FT_XSEC_BINS_BySlice FT_XSEC_BINS_CHANGE_BySlice FT_XSEC_PREPOSTAVG_BySlice...
        FT_XSEC_BINS_ByMouse FT_XSEC_BINS_CHANGE_ByMouse FT_XSEC_PREPOSTAVG_ByMouse
end

clear param_name plotchange start_dir save_dir AreaCutOff BLBinsExclude cond condition DATALABELS EBL ExcludeBaseline param...
    condition2exclude condition2plot nearVfar GroupBy TSERIESAB area_band Area_Bands bin condition...
    event_idx_post event_idx_pre lmk_idx param param_name...
    param_post param_pre param_unit per_change pulse_num ResNonRes Responder save_dir start_dir TB ThreshDir time_window...
    TS_combined TS_combined_PrePost TS_combined_change tseries UseAreaThresh WE YLIMS condition2plot condition2exclude plotchange...
    STAT stat WTS AB ExcludeRampingCells PostStimTime RampingCellThresh ResBasedOn