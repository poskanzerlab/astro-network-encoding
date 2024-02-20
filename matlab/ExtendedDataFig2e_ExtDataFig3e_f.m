%% Extended Data Fig. 2e & Extended Data Fig. 3e-f:
%  Histograms of event area and duration pre- and post-uncaging
%  Choose which cells to pull events from (uncaging cells or neighboring
%  cells)
%
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%  3. Run 'Fig2_dFFResponders.m': **only necessary for Extended Data Fig. 2e**
%  4. Run the following code:
%
%     Section A: Bin events into stationary bins of a chosen duration for each condition
%               & calculate event features per cell/region using the res file
%               for each tseries:
%               time_window = 120;
%     Section B: Plot histograms of chosen event features for the chosen cells/regions
%              within the designated time window pre- and post-stim for each condition
%              Use the parameters specified at the beginning of the section
%              for Ext. Data Fig. 2e & Ext. Data Fig. 3e-f
%
% Michelle Cahill 20240113
%% A. Bin events into stationary bins of a chosen duration for multiple conditions
% & calculate event features per cell/region using the res file for each tseries
% From 'Aggregate_MC_ftsFilter_MultCon_XSecBins_20220524.m' (this is
% similar to section E of 'Fig2_3_PreppingDataStruct.m', but ensures all time
% bins are of equal length (necessary when using 150 is not a multiple of 'time_window'). 

time_window = 120; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);
condition = fieldnames(mydata);

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            clear event_idx_pre event_idx_post
            uncaging_frame = mydata(tseries).(condition{cond}).uncagingframe;
            total_frames = mydata(tseries).(condition{cond}).totalframes;
            num_frames = round(time_window / mydata(tseries).(condition{cond}).SecPerFrame); % convert the time window from seconds to frames

            FramesCroppedBeginning = mod(uncaging_frame - 1, num_frames);
            pre_frame_cutoffs = FramesCroppedBeginning:num_frames:mydata(tseries).(condition{cond}).uncagingframe-1;
            num_bins_pre = numel(pre_frame_cutoffs) - 1; %How many time bins are there pre-stim?
            
            FramesCroppedEnd = mod(total_frames-(uncaging_frame-1), num_frames);
            post_frame_cutoffs = mydata(tseries).(condition{cond}).uncagingframe-1:num_frames:total_frames-FramesCroppedEnd;
            num_bins_post = numel(post_frame_cutoffs) - 1; %How many time bins are there post-stim?
             
            frame_cutoffs = [pre_frame_cutoffs,post_frame_cutoffs(2:end)];
            
            %find the indices of events occuring pre vs post stimulation for each
            %individual cell
            event_idx_pre = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), num_bins_pre);
            event_idx_post = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), num_bins_post);
            
            for region = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2)
                clear event_idx_0 event_idx event_time
                event_idx_0 = ~isnan(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx(:, region));
                event_idx = find(event_idx_0); %pulls out the indices of events in that region
                
                if isempty(event_idx) %if there aren't any events in this region, fill in that row with empty matrices
                    
                else
                    event_time = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event_idx)';
                    
                    for event = 1:length(event_idx)
                        event_placed = 0;
                        
                        for bin = 1:(num_bins_post + num_bins_pre)
                            if event_placed == 0
                                if frame_cutoffs(bin) < event_time(event,1) && event_time(event,1) <= frame_cutoffs(bin+1)
                                    if bin <= num_bins_pre
                                        event_idx_pre{region, bin} = [ event_idx_pre{region, bin}, event_idx(event)];
                                    else
                                        event_idx_post{region, bin-num_bins_pre} = [event_idx_post{region, bin-num_bins_pre}, event_idx(event)];
                                    end
                                    event_placed = 1;
                                end
                            end
                        end
                    end
                end
                clear event_idx_0 event_idx event_time event bin event_placed
            end
            clear region uncaging_frame frame_cutoffs num_bins_post num_bins_pre num_frames...
                post_frame_cutoffs pre_frame_cutoffs total_frames FramesCroppedBeginning FramesCroppedEnd
            
            mydata(tseries).(condition{cond}).(TB).event_idx_pre = event_idx_pre;
            mydata(tseries).(condition{cond}).(TB).event_idx_post = event_idx_post;
            
            %calculate events per bin for each
            %cell/region
            EventCt_PerBinPre = cell2mat(cellfun(@(x) numel(x), event_idx_pre, 'UniformOutput', 0));
            EventCt_PerBinPost = cell2mat(cellfun(@(x) numel(x), event_idx_post, 'UniformOutput', 0));
            
            %calculating percent change from the average baseline bins
            AvgEvtCt_Pre = mean(EventCt_PerBinPre,2);
            AvgEvtCt_PerChange = ((EventCt_PerBinPost - AvgEvtCt_Pre) ./ AvgEvtCt_Pre) * 100;
            
            %calculating percent change from the baseline bin immediately
            %before uncaging
            EvtCt_PerChange = ((EventCt_PerBinPost - EventCt_PerBinPre(:,end)) ./ EventCt_PerBinPre(:,end)) * 100;
            
            mydata(tseries).(condition{cond}).(TB).numevents_Pre = EventCt_PerBinPre;
            mydata(tseries).(condition{cond}).(TB).numevents_Post = EventCt_PerBinPost;
            mydata(tseries).(condition{cond}).(TB).AvgEvtCt_PerChange = AvgEvtCt_PerChange;
            mydata(tseries).(condition{cond}).(TB).EvtCt_PerChange = EvtCt_PerChange;
            clear EventCt_PerBinPre EventCt_PerBinPost AvgEvtCt_PerChange AvgEvtCt_Pre EvtCt_PerChange
            
            %area of each event divided between pre and post stim for each cell
            Area_PerBinPre = cellfun(@(x) mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(x), event_idx_pre, 'UniformOutput', 0);
            Area_PerBinPost = cellfun(@(x) mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(x), event_idx_post, 'UniformOutput', 0);
            AvgArea_PerBinPre = cell2mat(cellfun(@(x) mean(x), Area_PerBinPre, 'UniformOutput', 0));
            AvgArea_PerBinPost = cell2mat(cellfun(@(x) mean(x), Area_PerBinPost, 'UniformOutput', 0));
            
            %calculating percent change from the average baseline bins
            AvgArea_Pre = mean(AvgArea_PerBinPre,2);
            AvgArea_PerChange = ((AvgArea_PerBinPost - AvgArea_Pre) ./ AvgArea_Pre) * 100;
            
            %calculating percent change from the baseline bin immediately
            %before uncaging
            Area_PerChange = ((AvgArea_PerBinPost - AvgArea_PerBinPre(:,end)) ./ AvgArea_PerBinPre(:,end)) * 100;
            
            mydata(tseries).(condition{cond}).(TB).area_Pre = Area_PerBinPre;
            mydata(tseries).(condition{cond}).(TB).area_Post = Area_PerBinPost;
            mydata(tseries).(condition{cond}).(TB).avgarea_Pre = AvgArea_PerBinPre;
            mydata(tseries).(condition{cond}).(TB).avgarea_Post = AvgArea_PerBinPost;
            mydata(tseries).(condition{cond}).(TB).avgarea_PerChange = AvgArea_PerChange;
            mydata(tseries).(condition{cond}).(TB).area_PerChange = Area_PerChange;
            clear Area_PerBinPre Area_PerBinPost AvgArea_PerBinPre AvgArea_PerBinPost AvgArea_Pre AvgArea_PerChange Area_PerChange
            
            %duration (in seconds) of each event divided between pre and post stim for each cell
            duration = (mydata(tseries).(condition{cond}).res.ftsFilter.loc.t1 - mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0) * mydata(tseries).(condition{cond}).SecPerFrame;
            Duration_PerBinPre = cellfun(@(x) duration(x), event_idx_pre, 'UniformOutput', 0);
            Duration_PerBinPost = cellfun(@(x) duration(x), event_idx_post, 'UniformOutput', 0);
            AvgDuration_PerBinPre = cell2mat(cellfun(@(x) mean(x), Duration_PerBinPre, 'UniformOutput', 0));
            AvgDuration_PerBinPost = cell2mat(cellfun(@(x) mean(x), Duration_PerBinPost, 'UniformOutput', 0));
            
            %calculating percent change from the average baseline bins
            AvgDuration_Pre = mean(AvgDuration_PerBinPre,2);
            AvgDuration_PerChange = ((AvgDuration_PerBinPost - AvgDuration_Pre) ./ AvgDuration_Pre) * 100;
            
            %calculating percent change from the baseline bin immediately
            %before uncaging
            Duration_PerChange = ((AvgDuration_PerBinPost - AvgDuration_PerBinPre(:,end)) ./ AvgDuration_PerBinPre(:,end)) * 100;
            
            mydata(tseries).(condition{cond}).(TB).duration_Pre = Duration_PerBinPre;
            mydata(tseries).(condition{cond}).(TB).duration_Post = Duration_PerBinPost;
            mydata(tseries).(condition{cond}).(TB).avgduration_Pre = AvgDuration_PerBinPre;
            mydata(tseries).(condition{cond}).(TB).avgduration_Post = AvgDuration_PerBinPost;
            mydata(tseries).(condition{cond}).(TB).avgduration_PerChange = AvgDuration_PerChange;
            mydata(tseries).(condition{cond}).(TB).duration_PerChange = Duration_PerChange;
            clear Duration_PerBinPre Duration_PerBinPost AvgDuration_PerBinPre AvgDuration_PerBinPost duration AvgDuration_Pre AvgDuration_PerChange Duration_PerChange
            
            %dffMax (amplitude) of each event divided between pre and post stim for each cell
            dffMax_PerBinPre = cellfun(@(x) mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(x), event_idx_pre, 'UniformOutput', 0);
            dffMax_PerBinPost = cellfun(@(x) mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(x), event_idx_post, 'UniformOutput', 0);
            AvgdffMax_PerBinPre = cell2mat(cellfun(@(x) mean(x), dffMax_PerBinPre, 'UniformOutput', 0));
            AvgdffMax_PerBinPost = cell2mat(cellfun(@(x) mean(x), dffMax_PerBinPost, 'UniformOutput', 0));
            
            %calculating percent change from the average baseline bins
            AvgdffMax_Pre = mean(AvgdffMax_PerBinPre,2);
            AvgdffMax_PerChange = ((AvgdffMax_PerBinPost - AvgdffMax_Pre) ./ AvgdffMax_Pre) * 100;
            
            %calculating percent change from the baseline bin immediately
            %before uncaging
            dffMax_PerChange = ((AvgdffMax_PerBinPost - AvgdffMax_PerBinPre(:,end)) ./ AvgdffMax_PerBinPre(:,end)) * 100;
            
            mydata(tseries).(condition{cond}).(TB).dffMax_Pre = dffMax_PerBinPre;
            mydata(tseries).(condition{cond}).(TB).dffMax_Post = dffMax_PerBinPost;
            mydata(tseries).(condition{cond}).(TB).avgdffMax_Pre = AvgdffMax_PerBinPre;
            mydata(tseries).(condition{cond}).(TB).avgdffMax_Post = AvgdffMax_PerBinPost;
            mydata(tseries).(condition{cond}).(TB).avgdffMax_PerChange = AvgdffMax_PerChange;
            mydata(tseries).(condition{cond}).(TB).dffMax_PerChange = dffMax_PerChange;
            clear dffMax_PerBinPre dffMax_PerBinPost AvgdffMax_PerBinPre AvgdffMax_PerBinPost AvgdffMax_Pre AvgdffMax_PerChange dffMax_PerChange
            
            %growing propagation summed across all directions of each event divided between pre and post stim for each cell
            Sum_GrowingProp = sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall,2)';
            GrowingProp_PerBinPre = cellfun(@(x) Sum_GrowingProp(x), event_idx_pre, 'UniformOutput', 0);
            GrowingProp_PerBinPost = cellfun(@(x) Sum_GrowingProp(x), event_idx_post, 'UniformOutput', 0);
  
            mydata(tseries).(condition{cond}).(TB).GrowingProp_Pre = GrowingProp_PerBinPre;
            mydata(tseries).(condition{cond}).(TB).GrowingProp_Post = GrowingProp_PerBinPost;
            
            % Frequency of propagative events (>1um) 
            FreqProp_pre = cell2mat(cellfun(@(x) sum(x>1), GrowingProp_PerBinPre, 'UniformOutput', 0));
            FreqProp_post = cell2mat(cellfun(@(x) sum(x>1), GrowingProp_PerBinPost, 'UniformOutput', 0));
            
            mydata(tseries).(condition{cond}).(TB).FreqGrowingProp_Pre = FreqProp_pre;
            mydata(tseries).(condition{cond}).(TB).FreqGrowingProp_Post = FreqProp_post;
            
            % Indices of propagative events (>1um)
            internal_propIdx_pre = cellfun(@(x) find(x>1), GrowingProp_PerBinPre, 'UniformOutput', 0);
            internal_propIdx_post = cellfun(@(x) find(x>1), GrowingProp_PerBinPost, 'UniformOutput', 0);
            global_propIdx_pre = cellfun(@(x, y) x(y), event_idx_pre, internal_propIdx_pre, 'UniformOutput', 0);
            global_propIdx_post = cellfun(@(x, y) x(y), event_idx_post, internal_propIdx_post, 'UniformOutput', 0);
            
            mydata(tseries).(condition{cond}).(TB).event_idx_pre_prop = global_propIdx_pre;
            mydata(tseries).(condition{cond}).(TB).event_idx_post_prop =  global_propIdx_post;
            clear internal_propIdx_pre internal_propIdx_post global_propIdx_pre global_propIdx_post
            
            % Indices of static events (<= 1um)
            internal_statIdx_pre = cellfun(@(x) find(x<=1), GrowingProp_PerBinPre, 'UniformOutput', 0);
            internal_statIdx_post = cellfun(@(x) find(x<=1), GrowingProp_PerBinPost, 'UniformOutput', 0);
            global_statIdx_pre = cellfun(@(x, y) x(y), event_idx_pre, internal_statIdx_pre, 'UniformOutput', 0);
            global_statIdx_post = cellfun(@(x, y) x(y), event_idx_post, internal_statIdx_post, 'UniformOutput', 0);
            
            mydata(tseries).(condition{cond}).(TB).event_idx_pre_stat = global_statIdx_pre;
            mydata(tseries).(condition{cond}).(TB).event_idx_post_stat =  global_statIdx_post;
            clear internal_statIdx_pre internal_statIdx_post global_statIdx_pre global_statIdx_post...
                Sum_GrowingProp GrowingProp_PerBinPre GrowingProp_PerBinPost FreqProp_pre FreqProp_post
            
            if ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.region.landmarkDir)
                %propagation of events towards the uncaging site pre and post stim for each cell
                PropToward = mydata(tseries).(condition{cond}).res.ftsFilter.region.landmarkDir.chgToward';
                PropToward_PerBinPre = cellfun(@(x) PropToward(x), event_idx_pre, 'UniformOutput', 0);
                PropToward_PerBinPost = cellfun(@(x) PropToward(x), event_idx_post, 'UniformOutput', 0);
                
                mydata(tseries).(condition{cond}).(TB).PropToward_Pre = PropToward_PerBinPre;
                mydata(tseries).(condition{cond}).(TB).PropToward_Post = PropToward_PerBinPost;
                clear PropToward PropToward_PerBinPre PropToward_PerBinPost
                
                %propagation of events away the uncaging site pre and post stim for each cell
                PropAway = mydata(tseries).(condition{cond}).res.ftsFilter.region.landmarkDir.chgAway';
                PropAway_PerBinPre = cellfun(@(x) PropAway(x), event_idx_pre, 'UniformOutput', 0);
                PropAway_PerBinPost = cellfun(@(x) PropAway(x), event_idx_post, 'UniformOutput', 0);
                
                mydata(tseries).(condition{cond}).(TB).PropAway_Pre = PropAway_PerBinPre;
                mydata(tseries).(condition{cond}).(TB).PropAway_Post = PropAway_PerBinPost;
                clear PropAway PropAway_PerBinPre PropAway_PerBinPost
            end
            clear event_idx_pre event_idx_post
        end
    end
end

clear tseries time_window TB tseries cond
%% B. Histograms of chosen event features for the chosen cells/regions
%   within the designated time window pre- and post-stim for each condition
%  1. Aggregate together event features for each condition pre- and post-stim
%  2. Histograms of chosen event features for the chosen cells/regions
%     within the designated time window pre- and post-stim for each condition
%     Running a rank sum test between pre- and post-stim event features
%
%  For Extended Data Fig. 2e: Event area and duration 120-0s pre- and
%  0-120s post-stim in cells directly stimulated with NTs that increase
%  fluorescence post-stim (responders)
%       time_window = 120; WhichCells = 1; UseAreaThresh = 1; ResNonRes = 1; 
%       PostStimTime = []; ResBasedOn = 'UncagingEvents';
%       ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; ExcludeBaseline = 0;
%       Exclude0 = 1; param_name = {'area', 'duration'};
%
%  For Extended Data Fig. 3e-f: Event area and duration 120-0s pre- and
%  0-120s post-stim in all neighboring cells from all t-series/recordings
%  in each condition
%       time_window = 120; WhichCells = 2; UseAreaThresh = 0; 
%       ResNonRes = 1; PostStimTime = []; ResBasedOn = 'UncagingEvents'; <-these parameters aren't relevant when UseAreaThresh = 0;
%       ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; ExcludeBaseline = 0;
%       Exclude0 = 1; param_name = {'area', 'duration'};

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig3e_f';

condition = fieldnames(mydata);
newCondition = condition; 
BaselineComp = struct();
pulse_num = 10;

time_window = 120; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

WhichCells = 2; %1 to plot uncaging cells, 2 to plot neighboring cells with some activity
if WhichCells == 1
    WE = 'UncagingEvents';
    WC = 'uncaging';
elseif WhichCells == 2
    WE = 'NeighboringEvents';
    WC = 'AvgNeighbor';
end

UseAreaThresh = 0; %0= take all tseries from that pulse number, 1= take only tseries exceeding a specific threshold
ResNonRes = 1; %0 to plot only non-responders, 1 to plot responders (if UseAreaThresh == 1)
PostStimTime = []; %time period post-stim for determining responders
ResBasedOn = 'UncagingEvents'; % If UseAreaThresh == 1, should responder TS be based on response of uncaging cell or neighboring network? 'NeighboringEvents' or 'UncagingEvents'

ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 0; %number of seconds at the beginning of baseline to exclude (must be a multiple of time_window)
Exclude0 = 1; %value of 1 will exclude 0 values from the histograms

if UseAreaThresh == 1
    if ResNonRes == 1
        WTS = 'Responders';
    else
        WTS = 'NonResponders';
    end
else
    WTS = 'AllTS';
end

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

% 1. Aggregate together event features for each condition pre- and post-stim
param_name = {'area', 'duration'}; %'numevents',  'GrowingProp', 'numevents', 'propGrowOverall_PreVPost', 'propShrinkOverall_PreVPost', 'dffMax'};
for param = 1:length(param_name)
    for cond = 1:length(condition)
        temp_pre = [];
        temp_post = [];
        
        for tseries = TSERIES.(condition{cond}) %select_TS %1:length(mydata)
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
                
                if WhichCells == 1
                    Nonzero_cells = find(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk == 1);
                    Nonzero_cells = setdiff(Nonzero_cells, ramping_cell_idx);
                elseif WhichCells == 2
                    Cell_idx = find(isnan(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk)); %Identifying the neighboring cells
                    sum_evts = sum(mydata(tseries).(condition{cond}).numevents_PreVPost,2); %Sum the events pre & post in each cell
                    %                     change_inf = mean(mydata(tseries).(condition{cond}).(TB).(per_change),2);
                    %                     Nonzero_cells = Cell_idx(sum_evts(Cell_idx) > 0 & isfinite(change_inf(Cell_idx))); %Find the indices of neighboring cells with at least one event and without infinity as a value
                    Nonzero_cells = Cell_idx(sum_evts(Cell_idx) > 0); %Find the indices of neighboring cells with at least one event
                    Nonzero_cells = setdiff(Nonzero_cells, ramping_cell_idx);
                end
                ExcludeBaselineBins = ExcludeBaseline / time_window;
                
                param_Pre = mydata(tseries).(condition{cond}).(TB).(strcat(param_name{param}, '_Pre'))(Nonzero_cells, ExcludeBaselineBins+1:end);
                param_Post = mydata(tseries).(condition{cond}).(TB).(strcat(param_name{param}, '_Post'))(Nonzero_cells, :);
                
                param_Pre_vert = reshape(param_Pre, [], 1);
                param_Post_vert = reshape(param_Post, [], 1);
                
                if iscell(param_Pre_vert)
                    param_Pre_vert = cell2mat(param_Pre_vert')';
                    param_Post_vert = cell2mat(param_Post_vert')';
                end
                
                temp_pre = [temp_pre; param_Pre_vert(~isinf(param_Pre_vert)&~isnan(param_Pre_vert))]; %take only baseline events
                temp_post = [temp_post;  param_Post_vert(~isinf(param_Post_vert)&~isnan(param_Post_vert))]; %take only post-stim events
                clear param_Pre param_Post param_Pre_vert param_Post_vert Cell_idx sum_evts change_inf Nonzero_cells ramping_cell_idx ExcludeBaselineBins
            end
        end
        BaselineComp.Pre.(param_name{param}).(newCondition{cond}) = temp_pre;
        BaselineComp.Post.(param_name{param}).(newCondition{cond}) = temp_post;
        clear temp_pre temp_post pulse_cond select_TS
    end
end
clear cond tseries

% 2. Histograms of chosen event features for the chosen cells/regions
% within the designated time window pre- and post-stim for each condition
% Running a rank sum test between pre- and post-stim event features
for param = 1:length(param_name)
    if strcmp((param_name{param}), 'area')
        BINWIDTH = 1;
        BW = '1';
        XLIMS = [0 200];
    elseif strcmp((param_name{param}), 'duration')
        BINWIDTH = 1;
        BW = '1';
        XLIMS = [0 60];
    elseif strcmp((param_name{param}), 'dffMax')
        BINWIDTH = 0.25;
        BW = 'pt25';
        XLIMS = [0 25];
    end

    for cond = 1:length(condition)
        PRE = BaselineComp.Pre.(param_name{param}).(newCondition{cond});
        POST = BaselineComp.Post.(param_name{param}).(newCondition{cond});
        
        if Exclude0 == 1
            PRE = PRE(PRE>0);
            POST = POST(POST>0);
        end
        
        if contains(newCondition{cond}, 'GABA')
            POST_COLOR = [1 0 1];
        elseif contains(newCondition{cond}, 'Glu')
            POST_COLOR = [0 1 1];
        end
        
        ranksum_pval = ranksum(PRE, POST);
         
        TITLE = {sprintf('%s %s %s pre- and post-stim', indicator, strrep(newCondition{cond}, '_', ' '), param_name{param}),...
             sprintf('%s cells %s (%dS poststim) Exclude 0 vals(%d)', WC, WTS, PostStimTime, Exclude0),...
             sprintf('Exclude Ramping Cells %d, Exclude Baseline %dS', ExcludeRampingCells, ExcludeBaseline),...
             sprintf('Rank sum: p = %.15f', ranksum_pval)};
         YLABEL = 'probability';
         XLABEL = sprintf('%s (um)', param_name{param});
        
        %histograms
        figure()
        h_pre = histogram(PRE);
        hold on
        h_post = histogram(POST);
        h_pre.Normalization = 'probability';
        h_pre.BinWidth = BINWIDTH;
        h_pre.FaceColor = [0.5 0.5 0.5];
        h_pre.EdgeColor = 'none';
        
        h_post.Normalization = 'probability';
        h_post.BinWidth = BINWIDTH; 
        h_post.FaceColor = POST_COLOR;
        h_post.EdgeColor = 'none';
        title(TITLE)
        ylabel(YLABEL)
        xlabel(XLABEL)
        xlim(XLIMS)
        legend({'pre-stim', 'post-stim'})
        
        cd(save_dir)
        saveas(gcf, sprintf('%s_%s%dS_HistPreVPost_%sBins_ExRamp%dExBl%dSExZeros%d_%s%s%dSpost.tif', newCondition{cond}, param_name{param}, time_window, BW, ExcludeRampingCells, ExcludeBaseline, Exclude0, WC, WTS, PostStimTime))
        saveas(gcf, sprintf('%s_%s%dS_HistPreVPost_%sBins_ExRamp%dExBl%dSExZeros%d_%s%s%dSpost.svg', newCondition{cond}, param_name{param}, time_window, BW, ExcludeRampingCells, ExcludeBaseline, Exclude0, WC, WTS, PostStimTime))
        close
        cd(start_dir)

        clear h_pre h_post ranksum_pval TITLE STAT struct2plot POST_COLOR...
            PRE POST XLABEL YLABEL MAX EDGES cdf2plot PrePost p TITLE_pvals
    end
    clear XLIMS BINWIDTH BW
end

clear param cond WTS WC param_name Exclude0 save_dir start_dir...
    condition newCondition param_name BaselineComp PlotBaseline param pulse_num pulse_cond...
    ResNonRes UseAreaThresh time_window TB WhichCells ExcludeBaseline ExcludeRampingCells...
    XLABEL YLABEL PostStimTime RampingCellThresh ResBasedOn TSERIES WE