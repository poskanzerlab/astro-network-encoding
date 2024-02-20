%% Preparing data structs with AQuA res files from 2P uncaging experiments
%  Figure 2 & 3: 2P uncaging of GABA and gluamate and 2P imaging of CytoGCaMP6f
%
% Michelle Cahill 20240110
%% Section A: preparing mydata struct
% Modified from Aggregate_MC_ftsFilter_uncMod_MultCond_20210302.m
% 1. Record the frame the stimulus occurs ('uncagingframe') and other uncaging
%    parameters (using 'UncLaserVoltageRecording') 
% 2. Change region indices in each  res file to match any manual re-naming done during AQuA event detection
% 3. Note down imaging parameters (pulled from AQuA res file for each
%    recording)
%
%  For all Cyto-GCaMP datasets:
%    experiment_type = 'uncaging'; MultUncagingReps = 0;
%  For GluSnFR dataset:
%    experiment_type = 'uncaging'; MultUncagingReps = 1;

experiment_type = 'uncaging'; % 'uncaging', 'BathApp'
if strcmp(experiment_type, 'uncaging')
    MultUncagingReps = 1; %0 if there was only one train of uncaging pulses in the recording, 1 if there were multiple trains of uncaging pulses in the same recording (as in GluSnFR uncaging control data)
elseif strcmp(experiment_type, 'BathApp')
    method_AgonistEntry = 'MaxCurve'; %What method should be used to indicate the frame agonist enters imaging chamber 'MaxCurve', 'Threshold'? Run on Alexa594 trace
end

condition = fieldnames(mydata);

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            if strcmp(experiment_type, 'uncaging') %Record timing of uncaging using the voltage recording
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
            
            % Fill in the structured array of mydata with the res file, the time of
            % uncaging and the total time of the t-series
            mydata(tseries).(condition{cond}).totalframes = mydata(tseries).(condition{cond}).res.opts.sz(1,3);
            mydata(tseries).(condition{cond}).SecPerFrame = mydata(tseries).(condition{cond}).res.opts.frameRate;
            
            % Change region indices to match any manual re-naming done during AQuA event detection
            if ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.region)
                [mydata(tseries).(condition{cond}).res.ftsFilter.region.cellOrgIdx, mydata(tseries).(condition{cond}).res.ftsFilter.region.cell] = MatchRegionIdxToCellName(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell);
            end
            
            % Add in details about stimulation calculated above
            if strcmp(experiment_type, 'uncaging')
                if ~isempty(r)
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

clear cond MultUncagingReps
%% Section B. Identify pulse numbers in dataset and which recordings/tseries are included in each group, excluding any designated files 
% From Aggregate_MC_ftsFilter_uncMod_MultCond_20210302.m &
% ExcludeTS_FileNames_MC20230803.m
% create pulse_x_idx struct

if ~exist('RecordingsToExclude', 'var') %if all files are to be included
    data_sets = {}; %{'WT'}; %{'WT', 'CGP', 'LY', 'LaserUncagingControl'}; {'WT', 'Cx43', 'CBX'}; {'RuBiGlu_R1', 'RuBiGlu_R2', 'RuBiGlu_R3'};
    pulse_x_idx_Full = struct();
    pulse_x_idx = struct();
    
    %Which pulse condition are included in the datasets?
    NUM_PULSES = [];
    for cond = 1:length(condition)
        for tseries = 1:length(mydata)
            if ~isempty(mydata(tseries).(condition{cond}))
                NUM_PULSES = [NUM_PULSES; mydata(tseries).(condition{cond}).num_pulses];
            end
        end
    end
    NUM_PULSES = unique(NUM_PULSES);
    
    %For each individual condition and pulse number, find tseries indices
    for cond = 1:length(condition)
        for n = 1:length(NUM_PULSES)
            temp_ts_idx = [];
            for tseries = 1:length(mydata)
                if ~isempty(mydata(tseries).(condition{cond})) && mydata(tseries).(condition{cond}).num_pulses == NUM_PULSES(n)
                    temp_ts_idx = [temp_ts_idx, tseries];
                end
            end
            temp_cond = sprintf('%s_%d',condition{cond}, NUM_PULSES(n));
            pulse_x_idx_Full.(temp_cond) = temp_ts_idx;
            clear temp_cond temp_ts_idx
        end
    end
    clear cond n tseries
    
    %Between conditions (RuBiGABA_WT and RuBiGlu_WT, for example), if all
    %tseries indices are the same, combine the indices into one (WT_10, for
    %example); only used for WT, Cx43 and CBX dataset
    FN = fieldnames(pulse_x_idx_Full);
    if ~isempty(data_sets) && contains(data_sets{1,2}, 'Cx43')
        for d = 1:length(data_sets)
            DATA_SET = data_sets{d};
            for n = 1:length(NUM_PULSES)
                temp_FieldNames = FN(endsWith(FN, sprintf('%s_%d',DATA_SET, NUM_PULSES(n))));
                %         temp_FieldNames = FN(contains(FN, sprintf('%s_%d',DATA_SET, NUM_PULSES(n))));
                TS_cond1 = pulse_x_idx_Full.(temp_FieldNames{1});
                TS_cond2 = pulse_x_idx_Full.(temp_FieldNames{2});
                NotCommon = setxor(TS_cond1, TS_cond2);
                
                %Are the tseries not common between the two conditions simply
                %empty?
                EMPTY = [];
                MEMBER1 = NotCommon(~ismember(NotCommon, TS_cond1));
                MEMBER2 = NotCommon(~ismember(NotCommon, TS_cond2));
                for m = MEMBER1
                    f = extractBefore(temp_FieldNames{1}, sprintf('_%d', NUM_PULSES(n)));
                    if isempty(mydata(m).(f))
                        EMPTY = [EMPTY, m];
                    end
                    clear f
                end
                clear m
                for m = MEMBER2
                    f = extractBefore(temp_FieldNames{2}, sprintf('_%d', NUM_PULSES(n)));
                    if isempty(mydata(m).(f))
                        EMPTY = [EMPTY, m];
                    end
                    clear f
                end
                clear m
                NOT_EMPTY_TS = setxor(EMPTY, NotCommon);
                clear EMPTY MEMBER1 MEMBER2
                
                if isempty(NotCommon) || isempty(NOT_EMPTY_TS)
                    pulse_x_idx.(sprintf('%s_%d',DATA_SET, NUM_PULSES(n))) = TS_cond1;
                else
                    pulse_x_idx.(temp_FieldNames{1}) = TS_cond1;
                    pulse_x_idx.(temp_FieldNames{2}) = TS_cond2;
                end
                clear TS_cond1 TS_cond2 NOT_EMPTY_TS NotCommon temp_FieldNames
            end
            clear n DATA_SET
        end
    else
        pulse_x_idx = pulse_x_idx_Full;
    end
else %if some files are to be excluded
    pulse_x_idx = struct(); %Recording the tseries idx for each concentration or treatment value in each condition, exlucding specific days
    condition = fieldnames(mydata);
    NUM_PULSES = [];
    for cond = 1:length(condition)
        for tseries = 1:length(mydata)
            if ~isempty(mydata(tseries).(condition{cond}))
                NUM_PULSES = [NUM_PULSES; mydata(tseries).(condition{cond}).num_pulses];
            end
        end
    end
    NUM_PULSES = unique(NUM_PULSES);
    
    for cond = 1:length(condition)
        for n = 1:length(NUM_PULSES)
            TS_temp = [];
            temp_cond = sprintf('%s_%d',condition{cond}, NUM_PULSES(n));
            for tseries = 1:length(mydata)
                if ~isempty(mydata(tseries).(condition{cond})) && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
                    EXCLUDE_RECORDING = []; %is this t-series a slice that should be excluded?
                    for r = 1:length(RecordingsToExclude)
                        EXCLUDE_RECORDING = [EXCLUDE_RECORDING, contains(mydata(tseries).(condition{cond}).file, RecordingsToExclude{r})];
                    end
                    EXCLUDE_RECORDING = sum(EXCLUDE_RECORDING);
                    
                    if EXCLUDE_RECORDING == 0
                        TS_temp = [TS_temp, tseries];
                    end
                end
            end
            pulse_x_idx.(temp_cond) = TS_temp;
            clear concentration_TS_idx TS_temp tseries EXCLUDE_RECORDING r temp_cond
        end
    end
    clear cond n ans
end

clear d FN NUM_PULSES
%% Section C: Record the indices of "ramping cells", that should be
% excluded, for each recording. Uses the struct 'ramping_cells'
% Ramping cells significantly increase or decrease event frequency at baseline
% From ExcludeRampingCellsMultThresh_MC20220314.m

ramping_fields = fieldnames(ramping_cells);
for f = 1:length(ramping_fields)
    TS = unique(ramping_cells.(ramping_fields{f}).filename); %pull out name for each tseries
    for tseries = 1:length(TS)
        TSCells = cellfun(@(x) strcmp(x, TS{tseries}), ramping_cells.(ramping_fields{f}).filename); %Find all of the cells from that t-series
        pvals = ramping_cells.(ramping_fields{f}).p(TSCells); %record the ramping p-values
        cell_idx = ramping_cells.(ramping_fields{f}).cell_recording(TSCells); %record the AQuA cell id
        num_evt = ramping_cells.(ramping_fields{f}).n_events(TSCells); %record how many events took place in each cell
%         exclude_cell_idx_p05 = cell_idx(pvals <= 0.05 & num_evt > 5); %record the cell ids with ramping p-values <=0.05 and greater than 5 events
        exclude_cell_idx_p1 = cell_idx(pvals <= 0.1 & num_evt > 5); %record the cell ids with ramping p-values <=0.2 and greater than 5 events
%         exclude_cell_idx_p2 = cell_idx(pvals <= 0.2 & num_evt > 5); %record the cell ids with ramping p-values <=0.2 and greater than 5 events
        
        % figure out which tseries in mydata this corresponds to
        condition = fieldnames(mydata);
        foundTS = 0;
        while foundTS == 0
            for cond = 1:length(condition)
                for t = 1:length(mydata)
                    if ~isempty(mydata(t).(condition{cond}))
                        FN = mydata(t).(condition{cond}).file; %extract the file name of the AQuA res file
                        RCD = extractBefore(TS{tseries}, '_'); %extract the date of the tseries of the ramping cells;
                        
                        if sum(contains(condition, 'R1')) || sum(contains(condition, 'R2')) || sum(contains(condition, 'R3'))
                            % To be used on multi-round glutamate uncaging data
                            partialFN = extractAfter(TS{tseries}, 'RuBiG');
                            RCTN = extractBetween(partialFN, '_', '_'); %extract the t-series number
                            RCTN = extractBefore(RCTN, 4); %in case there were extra characters after the tseries number
                            RCTN = strcat('-', RCTN);
                            
                            if isempty(RCTN) %This catches file names in ramping cells with a '-' before the t-series number instead of '_'
                                RCTN = extractBetween(partialFN, '-', '_'); %extact the t-series number
                                RCTN = extractBefore(RCTN, 4); %in case there were extra characters after the tseries number
                                RCTN = strcat('-', RCTN);
                            end   
                        else
                            % To be used on all other recordings
                            RCTN = extractBetween(TS{tseries}, '-', '_');
                            RCTN = extractBefore(RCTN, 4); %in case there were extra characters after the tseries number
                            RCTN = strcat('-', RCTN);
                        end
                        
                        if contains(FN, RCD) && contains(FN, RCTN) %the file name for the AQuA res file must include the date and t-series number of the ramping cell info
                            foundTS = 1;
%                             mydata(t).(condition{cond}).RampingCellIdx.pval_05 = exclude_cell_idx_p05;
                            mydata(t).(condition{cond}).RampingCellIdx.pval_1 = exclude_cell_idx_p1;
%                             mydata(t).(condition{cond}).RampingCellIdx.pval_2 = exclude_cell_idx_p2;
                        end
                        
                        if isempty(RCTN) %if something is wrong and the t-series number is empty, move on (this t-sereis will be identified as one without ramping cells ID'ed
                            foundTS = 1;
                        end
                        clear FN RCD RCTN partialFN
                    end
                end
            end
        end
        clear foundTS condition cond t TSCells pvals cell_idx exclude_cell_idx_p05 exclude_cell_idx_p1 exclude_cell_idx_p2 num_evt
    end
    clear TS tseries
end

clear ramping_fields f
%% Section D: Calculate parameters for each condition per cell/region using the res file for each tseries
% From Aggregate_MC_ftsFilter_uncMod_MultCond_20210302.m
condition = fieldnames(mydata);
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond})) && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
            clear event_idx_pre event_idx_post event_time_sec event_time_frames
            uncaging_frame = mydata(tseries).(condition{cond}).uncagingframe;
            
            %find the indices of events occuring pre vs post stimulation for each
            %individual cell
            event_idx_pre = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 1);
            event_idx_post = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 1);
            event_time_sec = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 1);
            event_time_frames = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 1);
            
            for region = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2)
                clear event_idx_0 event_idx pre_idx post_idx event_time
                event_idx_0 = ~isnan(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx(:, region));
                event_idx = find(event_idx_0); %pulls out the indices of events in that region
                
                if isempty(event_idx) %if there aren't any events in this region, fill in that row with empty matrices
                    
                else
                    event_time = nan(length(event_idx),1); %find the frame in which each event occurred
                    
                    pre_idx = [];
                    post_idx = [];
                    for event = 1:length(event_idx)
                        event_time(event,1) = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event_idx(event));
                        
                        if event_time(event,1) <= uncaging_frame %break the indices into pre vs post stimulus
                            pre_idx = [pre_idx, event_idx(event)];
                        else
                            post_idx = [post_idx, event_idx(event)];
                        end
                    end
                    event_idx_pre(region) = {pre_idx};
                    event_idx_post(region) = {post_idx};
                    event_time_sec(region) = {event_time * mydata(tseries).(condition{cond}).SecPerFrame};
                    event_time_frames(region) = {event_time};
                end
            end
            mydata(tseries).(condition{cond}).event_idx_pre = event_idx_pre;
            mydata(tseries).(condition{cond}).event_idx_post = event_idx_post;
            mydata(tseries).(condition{cond}).event_times_sec = event_time_sec;
            mydata(tseries).(condition{cond}).event_times_frames = event_time_frames;
            
            %calculate events per minute pre and post stimulation for each
            %cell/region
            clear eventnumber_PreVPost
            eventnumber_PreVPost = nan(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2),2); %rows are cells/regions, column 1 pre-stim column 2 post-stim
            for region = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2)
                eventnumber_PreVPost(region, 1) = numel(event_idx_pre{region}) / ((mydata(tseries).(condition{cond}).uncagingframe * mydata(tseries).(condition{cond}).SecPerFrame)/60); %1st column is # of events/min pre-stim
                eventnumber_PreVPost(region, 2) = numel(event_idx_post{region}) / (((mydata(tseries).(condition{cond}).totalframes - mydata(tseries).(condition{cond}).uncagingframe) * mydata(tseries).(condition{cond}).SecPerFrame) / 60);%# of events/min post-stim
            end
            
            mydata(tseries).(condition{cond}).numevents_PreVPost = eventnumber_PreVPost;
            
            %area of each event divided between pre and post stim for each cell
            clear area_PreVPost avgarea_PreVPost
            area_PreVPost = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 2);
            avgarea_PreVPost = nan(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 2);
            for region = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2)
                area_PreVPost(region, 1) = {mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_pre{region})};
                area_PreVPost(region, 2) = {mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_post{region})};
                avgarea_PreVPost(region, 1) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_pre{region}));
                avgarea_PreVPost (region, 2) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_post{region}));
            end
            mydata(tseries).(condition{cond}).area_PreVPost = area_PreVPost;
            mydata(tseries).(condition{cond}).avgarea_PreVPost = avgarea_PreVPost;
            
            %duration (in seconds) of each event divided between pre and post stim for each cell
            clear duration_PreVPost avgduration_PreVPost
            duration_PreVPost = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 2);
            avgduration_PreVPost = nan(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 2);
            for region = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2)
                frames_pre = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t1(event_idx_pre{region}) - mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event_idx_pre{region});
                duration_PreVPost(region, 1) = {frames_pre * mydata(tseries).(condition{cond}).SecPerFrame};
                avgduration_PreVPost(region, 1) = mean(frames_pre * mydata(tseries).(condition{cond}).SecPerFrame);
                
                frames_post = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t1(event_idx_post{region}) - mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event_idx_post{region});
                duration_PreVPost(region, 2) = {frames_post * mydata(tseries).(condition{cond}).SecPerFrame};
                avgduration_PreVPost(region, 2) = mean(frames_post * mydata(tseries).(condition{cond}).SecPerFrame);
            end
            mydata(tseries).(condition{cond}).duration_PreVPost = duration_PreVPost;
            mydata(tseries).(condition{cond}).avgduration_PreVPost = avgduration_PreVPost;
            
            
            %dffMax (amplitude) of each event divided between pre and post stim for each cell
            clear dffMax_PreVPost avgdffMax_PreVPost
            dffMax_PreVPost = cell(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 2);
            avgdffMax_PreVPost = nan(size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2), 2);
            for region = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2)
                dffMax_PreVPost(region, 1) = {mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_pre{region})};
                dffMax_PreVPost(region, 2) = {mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_post{region})};
                avgdffMax_PreVPost(region, 1) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_pre{region}));
                avgdffMax_PreVPost (region, 2) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_post{region}));
            end
            mydata(tseries).(condition{cond}).dffMax_PreVPost = dffMax_PreVPost;
            mydata(tseries).(condition{cond}).avgdffMax_PreVPost = avgdffMax_PreVPost;
            
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
    
    clear currLocs xSize ySize zSize dffMax_PreVPost avgdffMax_PreVPost duration_PreVPost avgduration_PreVPost area_PreVPost avgarea_PreVPost eventnumber_PreVPost event_idx_0 event_idx pre_idx post_idx event_time event_idx_pre event_idx_post event_time_sec event_time_frames uncaging_frame tseries tseries region frames_post frames_pre event
end

clear cond
%% Section E: Bin events into static time bins of chosen duration
% From Aggregate_MC_ftsFilter_MultCond_XSecBins_20210304.m

% Using event indices within each static time bin
% Calculate parameters per cell/region using the res file for each tseries
time_window = 30; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

condition = fieldnames(mydata);

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            clear event_idx_pre event_idx_post
            uncaging_frame = mydata(tseries).(condition{cond}).uncagingframe;
            num_frames = round(time_window / mydata(tseries).(condition{cond}).SecPerFrame); % convert the time window from seconds to frames
            num_bins_pre = ceil((mydata(tseries).(condition{cond}).uncagingframe-1)/num_frames); %How many time bins are needed to cover the entire pre-period
            num_bins_post = ceil((mydata(tseries).(condition{cond}).totalframes - mydata(tseries).(condition{cond}).uncagingframe)/num_frames); %How many time bins are needed to cover the entire post-period
            
            pre_frame_cutoffs = num_frames:num_frames:mydata(tseries).(condition{cond}).uncagingframe-1;
            if ~(pre_frame_cutoffs(end) == mydata(tseries).(condition{cond}).uncagingframe-1)
                pre_frame_cutoffs = [pre_frame_cutoffs, mydata(tseries).(condition{cond}).uncagingframe-1];
            end
            
            post_frame_cutoffs = mydata(tseries).(condition{cond}).uncagingframe:num_frames:mydata(tseries).(condition{cond}).totalframes;
            if ~(post_frame_cutoffs(end) == mydata(tseries).(condition{cond}).totalframes)
                post_frame_cutoffs = [post_frame_cutoffs, mydata(tseries).(condition{cond}).totalframes];
            end
            
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
                                if event_time(event,1) <= frame_cutoffs(bin)
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
            clear region uncaging_frame frame_cutoffs num_bins_post num_bins_pre num_frames post_frame_cutoffs pre_frame_cutoffs
            
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