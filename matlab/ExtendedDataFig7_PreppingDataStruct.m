%% Extended Data Fig. 7: Prepare the data struct for calculating the change in event frequency post-uncaging 
% Bin events into static time bins of chosen duration 
% From Aggregate_MC_ftsFiler_MultCond_XSecBins_20220524.m

% This is similar to Section E of Fig2_3_PreppingDataStruct.m
% modified 20220524 to change pre_frame_cutoffs and post_frame_cutoffs,
% cutting off some frames at the beginning and end of the recording to have
% the same number of frames accounted for in each bin 
%% Using event indices within each static time bin
% Calculate parameters per cell/region using the res file for each tseries
time_window = 30; %time window (in seconds) taken before the stim and after the stim
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
