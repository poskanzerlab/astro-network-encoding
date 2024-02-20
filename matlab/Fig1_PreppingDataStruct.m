%% Preparing data structs with AQuA res files from 2P Ca2+ imaging experiments 
%  during receptor agonist bath application
%  Used to prepare data struct for Fig. 1: bath application of receptor agonist (Baclofen & t-ACPD)
%
% Michelle Cahill 20231208
%% Section A: preparing mydata struct
% From Aggregate_MC_ftsFilter_uncMod_MultCond_20210302.m
% Calculate the frame the stimulus occurs ('uncagingframe')
% Note down imaging parameters (pulled from AQuA res file for each
% recording)
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig1C';

experiment_type = 'BathApp'; % 'uncaging', 'BathApp'
if strcmp(experiment_type, 'uncaging')
    MultUncagingReps = 1; %0 if there was only one train of uncaging pulses in the recording, 1 if there were multiple trains of uncaging pulses (as in GluSnFR uncaging control data)
elseif strcmp(experiment_type, 'BathApp')
    method_AgonistEntry = 'MaxCurve'; %What method should be used to indicate the frame agonist enters imaging chamber 'MaxCurve', 'Threshold'? Run on Alexa594 trace
end

condition = fieldnames(mydata);

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            if strcmp(experiment_type, 'uncaging')
                r = mydata(tseries).(condition{cond}).UncLaserVoltageRecording; %load the voltage recording for this tseries
                if ~isempty(r)
                    if MultUncagingReps == 1
                        [interpulse_duration_ms, num_pulses, num_bursts, pulseduration_ms, uncaging_start_ms, uncaging_end_ms, uncaging_duration_ms, uncaging_frame_start, uncaging_frame_end, max_voltage ] = square_pulse_detection_MultiTrain_20200723( r, mydata(tseries).(condition{cond}).res.opts.frameRate );
                    else
                        [interpulse_duration_ms, num_pulses, pulseduration_ms, uncaging_start_ms, uncaging_end_ms, uncaging_duration_ms, uncaging_frame_start, uncaging_frame_end, max_voltage ] = square_pulse_detection_20200408( r, mydata(tseries).(condition{cond}).res.opts.frameRate );
                    end
                end
            elseif strcmp(experiment_type, 'BathApp')
                ALEXA_Trace = mydata(tseries).(condition{cond}).Alexa594.RawTraceXY;
                mydata(tseries).(condition{cond}).Alexa594.LAremoved = RemovingLightArtifacts(ALEXA_Trace(:,2), 50, 10);
                LA_removed = mydata(tseries).(condition{cond}).Alexa594.LAremoved;
                
                if strcmp(method_AgonistEntry, 'MaxCurve')
                    [IncDec, fitresult, gof, SigFitTrace] = createFit_sigmoid20210816(ALEXA_Trace(1:600,1), LA_removed(1:600));
                    mydata(tseries).(condition{cond}).Alexa594.CoE = coeffvalues(fitresult);
                    mydata(tseries).(condition{cond}).Alexa594.SigFitTrace = SigFitTrace; %y-vals of the fitted sigmoid curve
                    [t1, t2, t3] = PointsMaxCurvature_4CoE(mydata(tseries).(condition{cond}).Alexa594.CoE); %frame of increase as frame of maximum curvature in the sigmoid curve
                    mydata(tseries).(condition{cond}).Alexa594.FrameMaxCurvature = min([t1 t2 t3]);
                    
                    %                         % plot point of max curvature and save for later
                    %                         % inspection
                    %                         hold on
                    %                         if mydata(tseries).(condition{cond}).Alexa594.FrameMaxCurvature > 0
                    %                             plot(mydata(tseries).(condition{cond}).Alexa594.FrameMaxCurvature, SigFitTrace(round(mydata(tseries).(condition{cond}).Alexa594.FrameMaxCurvature)), 'mx', 'MarkerSize', 15)
                    %                             legend('filtered fluorescence trace', 'Fit Sigmoid Curve', 'Max Curve', 'Location', 'NorthWest' );
                    %                         else
                    %                             legend('filtered fluorescence trace', 'Fit Sigmoid Curve', 'Location', 'NorthWest' );
                    %                         end
                    %                         title(sprintf('Alexa 594 %s: tseries %d', condition{cond}, tseries))
                    %                         cd(save_dir)
                    %                         saveas(gcf, sprintf('Alexa594SigFitMaxCurve_%stseries%d.tif', condition{cond}, tseries))
                    %                         saveas(gcf, sprintf('Alexa594SigFitMaxCurve_%stseries%d.svg', condition{cond}, tseries))
                    %                         close
                    %                         cd(start_dir)
                elseif strcmp(method_AgonistEntry, 'Threshold')
                    Z = 3;
                    VARI = 'std';
                    mydata(tseries).(condition{cond}).Alexa594.FirstFrameAboveBL = FirstFrame_AboveBaselineMean(LA_removed, [1 300], Z, VARI, 375);
                    
                    %                         figure()
                    %                         plot(LA_removed)
                    %                         hold on
                    %                         plot(mydata(tseries).(condition{cond}).Alexa594.FirstFrameAboveBL, LA_removed(round(mydata(tseries).(condition{cond}).Alexa594.FirstFrameAboveBL)), 'mx', 'MarkerSize', 15)
                    %                         legend('fluorescence trace', sprintf('%d%s above BL', Z, VARI), 'Location', 'NorthWest');
                    %                         title(sprintf('Alexa 594 %s: tseries %d', condition{cond}, tseries))
                    %                         cd(save_dir)
                    %                         saveas(gcf, sprintf('Alexa594SigFit_%d%sAboveBL_%stseries%d.tif', Z, VARI, condition{cond}, tseries))
                    %                         saveas(gcf, sprintf('Alexa594SigFit_%d%sAboveBL_%stseries%d.svg', Z, VARI, condition{cond}, tseries))
                    %                         close
                    %                         cd(start_dir)
                end
                clear IncDec fitresult gof SigFitTrace t1 t2 t3 LA_removed Z VARI ALEXA_Trace
            end
            
            % Fill in the structured array of mydata with imaging
            % parameters from the AQuA res file
            mydata(tseries).(condition{cond}).totalframes = mydata(tseries).(condition{cond}).res.opts.sz(1,3);
            mydata(tseries).(condition{cond}).SecPerFrame = mydata(tseries).(condition{cond}).res.opts.frameRate;
            
            %Change region indices to match any manual re-naming (added
            %20230807)
            if ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.region)
                [mydata(tseries).(condition{cond}).res.ftsFilter.region.cellOrgIdx, mydata(tseries).(condition{cond}).res.ftsFilter.region.cell] = MatchRegionIdxToCellName(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell);
            end
            
            if strcmp(experiment_type, 'uncaging')
                mydata(tseries).(condition{cond}).uncagingframe = uncaging_frame_start;
                mydata(tseries).(condition{cond}).uncagingframe_end = uncaging_frame_end;
                mydata(tseries).(condition{cond}).uncaging_start_ms = uncaging_start_ms;
                mydata(tseries).(condition{cond}).uncaging_end_ms = uncaging_end_ms;
                mydata(tseries).(condition{cond}).uncaging_duration_ms = uncaging_duration_ms;
                mydata(tseries).(condition{cond}).pulseduration_ms = pulseduration_ms;
                mydata(tseries).(condition{cond}).num_pulses = num_pulses;
                mydata(tseries).(condition{cond}).interpulse_duration_ms = interpulse_duration_ms;
                mydata(tseries).(condition{cond}).pulsepower = max_voltage;
                if MultUncagingReps == 1
                    mydata(tseries).(condition{cond}).num_bursts = num_bursts; %To be used with square_pulse_detection_MultiTrain_20200723
                end
            elseif strcmp(experiment_type, 'BathApp')
                if strcmp(method_AgonistEntry, 'MaxCurve')
                    mydata(tseries).(condition{cond}).uncagingframe = floor(mydata(tseries).(condition{cond}).Alexa594.FrameMaxCurvature);
                elseif strcmp(method_AgonistEntry, 'Threshold')
                    mydata(tseries).(condition{cond}).uncagingframe = floor(mydata(tseries).(condition{cond}).Alexa594.FirstFrameAboveBL);
                end
            end
            clear f curr_files curr_dir r uncaging_recording interpulse_duration_ms num_pulses pulseduration_ms uncaging_start_ms uncaging_end_ms ...
                uncaging_duration_ms uncaging_frame_start uncaging_frame_end max_voltage num_bursts
        end
    end
    clear tseries
end

fprintf ('Done!')

cd(start_dir);
clear dirs fields cond
%% Section B: Calculate parameters for each condition for population-wide activity the res
% file for each tseries (to be used when no regions are defined)
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond})) && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
            clear event_idx_pre event_idx_post event_time_sec event_time_frames
            uncaging_frame = mydata(tseries).(condition{cond}).uncagingframe;
            
            %find the indices of events occuring pre vs post simulation
            pre_idx = [];
            post_idx = [];
            for event = 1:length(mydata(tseries).(condition{cond}).res.ftsFilter.basic.area) %loop through each event, assigning it to pre or post stim
                if mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event) <= uncaging_frame %break the indices into pre vs post stimulus
                    pre_idx = [pre_idx, event];
                else
                    post_idx = [post_idx, event];
                end
            end
            event_idx_pre = pre_idx;
            event_idx_post = post_idx;
            
            mydata(tseries).(condition{cond}).event_idx_pre = event_idx_pre;
            mydata(tseries).(condition{cond}).event_idx_post = event_idx_post;
            
            %calculate events per minute pre and post stimulation for each
            %cell/region
            clear eventnumber_PreVPost
            eventnumber_PreVPost = nan(1,2); %rows are cells/regions, column 1 pre-stim column 2 post-stim
            eventnumber_PreVPost(1, 1) = numel(event_idx_pre) / ((mydata(tseries).(condition{cond}).uncagingframe * mydata(tseries).(condition{cond}).SecPerFrame)/60); %1st column is # of events/min pre-stim
            eventnumber_PreVPost(1, 2) = numel(event_idx_post) / (((mydata(tseries).(condition{cond}).totalframes - mydata(tseries).(condition{cond}).uncagingframe)*mydata(tseries).(condition{cond}).SecPerFrame) / 60);%# of events/min post-stim
            
            mydata(tseries).(condition{cond}).numevents_PreVPost = eventnumber_PreVPost;
            
            %area of each event divided between pre and post stim for each cell
            clear area_PreVPost avgarea_PreVPost
            area_PreVPost = cell(1, 2);
            avgarea_PreVPost = nan(1, 2);
            
            area_PreVPost(1, 1) = {mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_pre)};
            area_PreVPost(1, 2) = {mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_post)};
            avgarea_PreVPost(1, 1) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_pre));
            avgarea_PreVPost (1, 2) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_post));
            
            mydata(tseries).(condition{cond}).area_PreVPost = area_PreVPost;
            mydata(tseries).(condition{cond}).avgarea_PreVPost = avgarea_PreVPost;
            
            %duration (in seconds) of each event divided between pre and post stim for each cell
            clear duration_PreVPost avgduration_PreVPost
            duration_PreVPost = cell(1, 2);
            avgduration_PreVPost = nan(1, 2);
            frames_pre = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t1(event_idx_pre) - mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event_idx_pre);
            duration_PreVPost(1, 1) = {frames_pre * mydata(tseries).(condition{cond}).SecPerFrame};
            avgduration_PreVPost(1, 1) = mean(frames_pre * mydata(tseries).(condition{cond}).SecPerFrame);
            
            frames_post = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t1(event_idx_post) - mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event_idx_post);
            duration_PreVPost(1, 2) = {frames_post * mydata(tseries).(condition{cond}).SecPerFrame};
            avgduration_PreVPost(1, 2) = mean(frames_post * mydata(tseries).(condition{cond}).SecPerFrame);
            mydata(tseries).(condition{cond}).duration_PreVPost = duration_PreVPost;
            mydata(tseries).(condition{cond}).avgduration_PreVPost = avgduration_PreVPost;
            
            %dffMax (amplitude) of each event divided between pre and post stim for each cell
            clear dffMax_PreVPost avgdffMax_PreVPost
            dffMax_PreVPost = cell(1, 2);
            avgdffMax_PreVPost = nan(1, 2);
            dffMax_PreVPost(1, 1) = {mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_pre)};
            dffMax_PreVPost(1, 2) = {mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_post)};
            avgdffMax_PreVPost(1, 1) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_pre));
            avgdffMax_PreVPost (1, 2) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_post));
            mydata(tseries).(condition{cond}).dffMax_PreVPost = dffMax_PreVPost;
            mydata(tseries).(condition{cond}).avgdffMax_PreVPost = avgdffMax_PreVPost;
            
            %propGrowOverall of each event divided between pre and post stim for each cell
            clear propGrowOverall_PreVPost avgpropGrowOverall_PreVPost
            propGrowOverall_PreVPost = cell(1, 2);
            avgpropGrowOverall_PreVPost = nan(1, 2);
            
            propGrowOverall_PreVPost(1, 1) = {sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall(event_idx_pre,:),2)'};
            propGrowOverall_PreVPost(1, 2) = {sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall(event_idx_post,:),2)'};
            avgpropGrowOverall_PreVPost(1, 1) = mean(sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall(event_idx_pre,:),2));
            avgpropGrowOverall_PreVPost (1, 2) = mean(sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall(event_idx_post,:),2));
            
            mydata(tseries).(condition{cond}).propGrowOverall_PreVPost = propGrowOverall_PreVPost;
            mydata(tseries).(condition{cond}).avgpropGrowOverall_PreVPost = avgpropGrowOverall_PreVPost;
            
            %propShrinkOverall of each event divided between pre and post stim for each cell
            clear propShrinkOverall_PreVPost avgpropShrinkOverall_PreVPost
            propShrinkOverall_PreVPost = cell(1, 2);
            avgpropShrinkOverall_PreVPost = nan(1, 2);
            
            propShrinkOverall_PreVPost(1, 1) = {sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propShrinkOverall(event_idx_pre,:),2)'};
            propShrinkOverall_PreVPost(1, 2) = {sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propShrinkOverall(event_idx_post,:),2)'};
            avgpropShrinkOverall_PreVPost(1, 1) = mean(sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propShrinkOverall(event_idx_pre,:),2));
            avgpropShrinkOverall_PreVPost (1, 2) = mean(sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propShrinkOverall(event_idx_post,:),2));
            
            mydata(tseries).(condition{cond}).propShrinkOverall_PreVPost = propShrinkOverall_PreVPost;
            mydata(tseries).(condition{cond}).avgpropShrinkOverall_PreVPost = avgpropShrinkOverall_PreVPost;
            
            %spatiotemporal information across all regions (adapted from Trisha's
            %code)
            clear xSize ySize zSize
            xSize = mydata(tseries).(condition{cond}).res.opts.sz(1); %pixels in rows (minus what AQuA takes off of borders) this confusingly is what we think of as usual Y-coordinates
            ySize = mydata(tseries).(condition{cond}).res.opts.sz(2); %pixels in columns (minus what AQuA takes off of borders) this confusingly is what we think of as usual X-coordinates
            zSize = mydata(tseries).(condition{cond}).res.opts.sz(3); %number of frames
            
            %make cell arrays for the x,y and z info for each event
            %xLocs 1xn n=#events where each element has the active x pixels
            %corresponding to the active frame in tLocs
            for event = 1:length(mydata(tseries).(condition{cond}).res.ftsFilter.loc.x3D) %loop through each event
                clear currLocs
                currLocs = mydata(tseries).(condition{cond}).res.ftsFilter.loc.x3D{event}; %getting the x3D file for a single event
                %             [mydata(tseries).(condition{cond}).xLocs{event}, mydata(tseries).(condition{cond}).yLocs{event}, mydata(tseries).(condition{cond}).tLocs{event}] = ind2sub([xSize ySize zSize], currLocs); %transform x3D into active x and y pixels per frame
                %             mydata(tseries).(condition{cond}).xLocs{event} = xSize - mydata(tseries).(condition{cond}).xLocs{event}; %flipped x idx starting
                [mydata(tseries).(condition{cond}).yLocs{event}, mydata(tseries).(condition{cond}).xLocs{event}, mydata(tseries).(condition{cond}).tLocs{event}] = ind2sub([xSize ySize zSize], currLocs); %transform x3D into active x and y pixels per frame
                mydata(tseries).(condition{cond}).yLocs{event} = xSize + 1 - mydata(tseries).(condition{cond}).yLocs{event}; %flipped y idx starting
            end
        end
    end
    
    clear currLocs xSize ySize zSize dffMax_PreVPost avgdffMax_PreVPost duration_PreVPost avgduration_PreVPost area_PreVPost avgarea_PreVPost eventnumber_PreVPost...
        event_idx_0 event_idx pre_idx post_idx event_time event_idx_pre event_idx_post event_time_sec event_time_frames uncaging_frame tseries t_series_idx region...
        frames_post frames_pre event propShrinkOverall_PreVPost avgpropShrinkOverall_PreVPost propGrowOverall_PreVPost avgpropGrowOverall_PreVPost
end

clear cond
%% Section C: Bin events into static time bins of chosen duration
% From Aggregate_MC_ftsFilter_MultCond_XSecBins_20220524.m
% cutting off some frames at the beginning and end of the recording to have
% the same number of frames accounted for in each bin 

% Using event indices within each static time bin
% Calculate parameters per cell/region using the res file for each tseries
time_window = 60; %time window (in seconds) taken before the stim and after the stim
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