%% Aggregate together individual AQuA res files 
%  into a structured array called 'mydata', each condition is a field
%  within mydata and rows contain individual AQuA res files & information specific 
%  to that recording
% 
% From Aggregate_MC_ftsFilter_uncMod_MultCond_20210302.m section 1
%
% All MATLAB workspaces provided in Dryad already have AQuA res files aggregated 
% into the structured array, 'mydata'. To create a new 'mydata' for a different 
% set of AQuA res files, run the following code before running either
% Fig1_PreppingDataStruct.m or Fig2_3_PreppingDataStruct.m. This script
% would take the place of running section A of Fig1_PreppingDataStruct.m or Fig2_3_PreppingDataStruct.m.
%
% ***********************************************************************
% Preparing data:
% 1. Run time-series imaging recordings through AQuA event detection. 
%    Each recording will have an individual AQuA output folder containing 
%    a MATLAB AQuA res file named with the *RecordingName*_AQuA.mat
%
% 2. A .csv file, indicating stimulus time(s) is required. The corresponding 
%    .csv file should be manually added to the corresponding AQuA output folder.
%    The 'experiment_type' determines how the stimulus time is calculated 
%    and what values should be recorded in the .csv file:
%
%     If experiment_type = 'uncaging', the .csv file should be a voltage
%     recording, with voltage exceeding a threshold when the uncaging laser
%     was active/pulsed. The filename should end with 'VoltageRecording_001.csv' 
%             Column 1: time from recording start (ms)
%             Column 2: Pockels cell voltage 
%
%     If experiment_type = 'BathApp', the .csv file should contain
%     fluorescence values from a second channel acquired during imaging,
%     corresponding to a fluorescent dye (such as Alexa594) added along
%     with agonist. These values can be measured in FIJI/ImageJ by recording
%     the z-axis profile of the second channel images.
%     The filename should end with 'ZAxisProfile.csv'' 
%             Column 1: frame number from recording start
%             Column 2: fluorescence values from 2nd imaging channel 
%     The time agonist enters the bath can be calculated using two different methods:
%          If method_AgonistEntry = 'MaxCurve': 
%                 1. A sigmoidal curve is fit to the fluorescence trace
%                 2. The points of maximum curvature are calculated 
%                 3. The earliest of these values is marked as the frame of agonist entry. 
%          If method_AgonistEntry = 'Threshold':
%                 1. A threshold is calculated on baseline values (3std above baseline mean) 
%                 2. The first frame >= threshold is marked as the frame of agonist entry.
%                    *This will currently begin searching frames beginning
%                     at frame 375 (final input to
%                    'FirstFrame_AboveBaselineMean') to ensure baseline
%                    frames won't be marked as frame of agonist entry
%
%  Variables to adjust in code below:
%     save_dir: pathname for the folder to save the variable 'dirs' 
%     experiment_type: a string that should be 'uncaging' or 'BathApp' to
%                      match the type of data collected (as described above in 2)
%     MultUncagingReps: 0 or 1 (if experiment_type = 'uncaging') to
%                       determine if a single uncaging train or multiple
%                       uncaging trains should be detected from the voltage
%                       recording.
%     method_AgonistEntry: a string that should be either 'MaxCurve' or
%                   'Threshold' (if experiment_type = 'BathBpp') to
%                   determine the method used to mark the frame of agonist
%                   entry into the bath (as described above in 2)
%     condition: a cell array with each cell containing a string that will
%                be a fieldname within 'mydata'. 
%     dirs: a structured array with a field for each condition. Each field
%           contains a 1xn cell array for n AQuA recordings to include in 
%           that condition. Each cell contains a pathname to an AQuA output
%           folder to be included in that condition. This variable will be
%           saved to the save_dir.
%       **Note, if AQuA recordings are paired across conditions, ensure their
%        pathnames are indexed in the same way within their respective cell
%        arrays (ie. if condition 1 and 2 are paired, dirs.(condition{1}){1,i} and
%        dirs.(condition{2}){1,i} should contain pathnames to the AQuA output
%        folders for the same field-of-view recorded in condition 1 and condition
%        2). Both paired and unpaired conditions can be included**
%     fields: a cell array containing fields within each AQuA res file to
%             delete. Individual AQuA res files are large & removing these
%             fields considerably reduces the size of each AQuA res file to 
%             make the entire workspace more manageable.
% **************************************************************************
% The following functions are called:
%
%  MatchRegionIdxToCellName.m
%
%  If experiment_type = 'uncaging':
%      square_pulse_detection_MultiTrain_20200723.m
%      square_pulse_detection_20200408.m
%  If experiment_type = 'BathApp'
%      RemovingLightArtifacts.m
%      createFit_sigmoid20210816.m
%      PointsMaxCurvature_4CoE.m
%      FirstFrame_AboveBaselineMean.m
%
% Michelle Cahill 20240131
%% Create a structured array called mydata with each condition as a field and each recording within that condition as a row

start_dir = pwd;
save_dir = 'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\PostEventDet'; %specify the folder where items should be saved

experiment_type = 'uncaging'; %Options: 'uncaging', 'BathApp'; This will determine the general method used to specify the stimulus time(s)
if strcmp(experiment_type, 'uncaging')
    MultUncagingReps = 1; %0 if there was only one train of uncaging pulses in the recording, 1 if there were multiple trains of uncaging pulses in the same recording (as in GluSnFR uncaging control data)
elseif strcmp(experiment_type, 'BathApp')
    method_AgonistEntry = 'MaxCurve'; %Options: 'MaxCurve', 'Threshold'. What method should be used to indicate the frame agonist enters imaging chamber? 
end

condition = {'RuBiGABA_WT', 'RuBiGlu_WT', 'RuBiGABA_CGP', 'RuBiGlu_LY', 'NoRuBi_LaserUncagingControl'}; %What should each condition within 'mydata' be named? Alternative names: {'Baclofen', 'tACPD', 'Baclofen_CBX', 'tACPD_CBX'};

dirs = struct(); % Create a structured array with a field for each condition containing pathnames to all AQuA output folders to be included in 'mydata'

% For each condition, list the pathnames for each AQuA res folder to be included
dirs.(condition{1}) = {'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_multipulse\AQuA\20200623_GluSnFRp14f_S1RuBiGlu1_1_100-001Cycle00001_Ch1',... 
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_multipulse\AQuA\20200623_GluSnFRp14f_S1RuBiGlu2_1_100-002Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_multipulse\AQuA\20200623_GluSnFRp14f_S1RuBiGlu3_1_100-003Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_multipulse\AQuA\20200623_GluSnFRp14f_S2RuBiGlu1_1_100-012Cycle00001_Ch1',...
    };

dirs.(condition{2}) = {'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_multipulse\AQuA\20200623_GluSnFRp14f_S1RuBiGlu1_10_100-007Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_multipulse\AQuA\20200623_GluSnFRp14f_S1RuBiGlu2_10_100-008Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_multipulse\AQuA\20200623_GluSnFRp14f_S1RuBiGlu3_10_100-009Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_multipulse\AQuA\20200623_GluSnFRp14f_S2RuBiGlu1_10_100-018Cycle00001_Ch1',...
    };

dirs.(condition{3}) = {'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_70AU_RuBiGlu2_10_100-015Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_70AU_RuBiGlu3_10_100-016Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_70AU_RuBiGlu3_10_100-020Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_70AU_RuBiGlu3_10_100-025Cycle00001_Ch1',...
    };

dirs.(condition{4}) = {'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_25AU_RuBiGlu2_10_100-006Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_25AU_RuBiGlu2_10_100-014Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_25AU_RuBiGlu3_10_100-017Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_25AU_RuBiGlu3_10_100-021Cycle00001_Ch1',...
    };

dirs.(condition{5}) = {'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_25AU_RuBiGlu2_10_100-006Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_25AU_RuBiGlu2_10_100-014Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_25AU_RuBiGlu3_10_100-017Cycle00001_Ch1',...
    'G:\Michelle_Event_Detection\EventDetectionGluSnFR2PUncaging_FollowingLaserAlignment\AQuA\20230608_GluSnFRp15m_S1R1_25AU_RuBiGlu3_10_100-021Cycle00001_Ch1',...
    };

mydata = struct();

fields = {'datOrg', 'evt', 'fts','riseLst','dF','featureTable','userFeatures',...
    'dffMat','dMat', 'dMatFilter', 'dMatFav', 'scl','ov', 'riseLstFilter', 'riseLstFav','stg',...
}; %fields from the AQuA res file to delete.

for cond = 1:length(condition)
    for t_series_idx = 1:length(dirs.(condition{cond}))
        curr_dir = dirs.(condition{cond}){t_series_idx};
        
        if sum(cell2mat(strfind(curr_dir, 'G:', 'ForceCellOutput', 1))) > 0 %Handle cases in which there isn't data for one tseries for one condition
            cd (curr_dir)
            
            f = dir('*.mat'); %load in the AQuA results file
            curr_files = load(f.name);
            
            if strcmp(experiment_type, 'uncaging')
                r = dir('*VoltageRecording_001.csv'); %load in the electrical recording for the uncaging laser if this file exists
                if ~isempty(r)
                    if MultUncagingReps == 1
                        [interpulse_duration_ms, num_pulses, num_bursts, pulseduration_ms, uncaging_start_ms, uncaging_end_ms, uncaging_duration_ms, uncaging_frame_start, uncaging_frame_end, max_voltage ] = square_pulse_detection_MultiTrain_20200723( r, curr_files.res.opts.frameRate );
                    else
                        [interpulse_duration_ms, num_pulses, pulseduration_ms, uncaging_start_ms, uncaging_end_ms, uncaging_duration_ms, uncaging_frame_start, uncaging_frame_end, max_voltage ] = square_pulse_detection_20200408( r, curr_files.res.opts.frameRate );
                    end
                end
            elseif strcmp(experiment_type, 'BathApp')
                r = dir('*ZAxisProfile.csv'); %load in the z-axis profile of Ch2 fluorscence trace (Alexa 594) if this file exists
                if ~isempty(r)
                    r = table2array(readtable(r.name));
                    mydata(t_series_idx).(condition{cond}).Alexa594.RawTraceXY = r;
                    mydata(t_series_idx).(condition{cond}).Alexa594.LAremoved = RemovingLightArtifacts(r(:,2), 50, 10);
                    LA_removed = mydata(t_series_idx).(condition{cond}).Alexa594.LAremoved;
                    
                    if strcmp(method_AgonistEntry, 'MaxCurve')
                        [IncDec, fitresult, gof, SigFitTrace] = createFit_sigmoid20210816(r(1:600,1), LA_removed(1:600));
                        %                     [IncDec, fitresult, gof, SigFitTrace] = createFit_sigmoid20210816(r(1:800,1), LA_removed(1:800));
                        mydata(t_series_idx).(condition{cond}).Alexa594.CoE = coeffvalues(fitresult);
                        mydata(t_series_idx).(condition{cond}).Alexa594.SigFitTrace = SigFitTrace; %y-vals of the fitted sigmoid curve
                        [t1, t2, t3] = PointsMaxCurvature_4CoE(mydata(t_series_idx).(condition{cond}).Alexa594.CoE); %frame of increase as frame of maximum curvature in the sigmoid curve
                        mydata(t_series_idx).(condition{cond}).Alexa594.FrameMaxCurvature = min([t1 t2 t3]);
                        
                        % plot point of max curvature and save for later
                        % inspection
                        hold on
                        if mydata(t_series_idx).(condition{cond}).Alexa594.FrameMaxCurvature > 0
                            plot(mydata(t_series_idx).(condition{cond}).Alexa594.FrameMaxCurvature, SigFitTrace(round(mydata(t_series_idx).(condition{cond}).Alexa594.FrameMaxCurvature)), 'mx', 'MarkerSize', 15)
                            legend('filtered fluorescence trace', 'Fit Sigmoid Curve', 'Max Curve', 'Location', 'NorthWest' );
                        else
                            legend('filtered fluorescence trace', 'Fit Sigmoid Curve', 'Location', 'NorthWest' );
                        end
                        title(sprintf('Alexa 594 %s: tseries %d', condition{cond}, t_series_idx))
                        cd(save_dir)
                        saveas(gcf, sprintf('Alexa594SigFitMaxCurve_%stseries%d.tif', condition{cond}, t_series_idx))
                        saveas(gcf, sprintf('Alexa594SigFitMaxCurve_%stseries%d.svg', condition{cond}, t_series_idx))
                        close
                        cd(start_dir)
                    elseif strcmp(method_AgonistEntry, 'Threshold')
                        Z = 3;
                        VARI = 'std';
                        mydata(t_series_idx).(condition{cond}).Alexa594.FirstFrameAboveBL = FirstFrame_AboveBaselineMean(LA_removed, [1 300], Z, VARI, 375);
                        
                        figure()
                        plot(LA_removed)
                        hold on
                        plot(mydata(t_series_idx).(condition{cond}).Alexa594.FirstFrameAboveBL, LA_removed(round(mydata(t_series_idx).(condition{cond}).Alexa594.FirstFrameAboveBL)), 'mx', 'MarkerSize', 15)
                        legend('fluorescence trace', sprintf('%d%s above BL', Z, VARI), 'Location', 'NorthWest');
                        title(sprintf('Alexa 594 %s: tseries %d', condition{cond}, t_series_idx))
                        cd(save_dir)
                        saveas(gcf, sprintf('Alexa594SigFit_%d%sAboveBL_%stseries%d.tif', Z, VARI, condition{cond}, t_series_idx))
                        saveas(gcf, sprintf('Alexa594SigFit_%d%sAboveBL_%stseries%d.svg', Z, VARI, condition{cond}, t_series_idx))
                        close
                        cd(start_dir)
                    end
                    
                    clear IncDec fitresult gof SigFitTrace t1 t2 t3 LA_removed Z VARI
                end
            end
            
            % Fill in the structured array of mydata with the res file, the time of
            % uncaging and the total time of the t-series
            
            mydata(t_series_idx).(condition{cond}).file = dirs.(condition{cond}){t_series_idx};
            mydata(t_series_idx).(condition{cond}).res = curr_files.res;
            mydata(t_series_idx).(condition{cond}).res = rmfield(mydata(t_series_idx).(condition{cond}).res, fields);
            mydata(t_series_idx).(condition{cond}).totalframes = mydata(t_series_idx).(condition{cond}).res.opts.sz(1,3);
            mydata(t_series_idx).(condition{cond}).SecPerFrame = mydata(t_series_idx).(condition{cond}).res.opts.frameRate;
            %Change region indices to match any manual re-naming (added
            %20230807)
            if ~isempty(curr_files.res.ftsFilter.region)
                [mydata(t_series_idx).(condition{cond}).res.ftsFilter.region.cellOrgIdx, mydata(t_series_idx).(condition{cond}).res.ftsFilter.region.cell] = MatchRegionIdxToCellName(curr_files.res.ftsFilter.region.cell);
            end
            
            if ~isempty(r)
                if strcmp(experiment_type, 'uncaging')
                    mydata(t_series_idx).(condition{cond}).uncagingframe = uncaging_frame_start;
                    mydata(t_series_idx).(condition{cond}).uncagingframe_end = uncaging_frame_end;
                    mydata(t_series_idx).(condition{cond}).uncaging_start_ms = uncaging_start_ms;
                    mydata(t_series_idx).(condition{cond}).uncaging_end_ms = uncaging_end_ms;
                    mydata(t_series_idx).(condition{cond}).uncaging_duration_ms = uncaging_duration_ms;
                    mydata(t_series_idx).(condition{cond}).pulseduration_ms = pulseduration_ms;
                    mydata(t_series_idx).(condition{cond}).num_pulses = num_pulses;
                    mydata(t_series_idx).(condition{cond}).interpulse_duration_ms = interpulse_duration_ms;
                    mydata(t_series_idx).(condition{cond}).pulsepower = max_voltage;
                    if MultUncagingReps == 1
                        mydata(t_series_idx).(condition{cond}).num_bursts = num_bursts; %To be used with square_pulse_detection_MultiTrain_20200723
                    end
                elseif strcmp(experiment_type, 'BathApp')
                    if strcmp(method_AgonistEntry, 'MaxCurve')
                        mydata(t_series_idx).(condition{cond}).uncagingframe = floor(mydata(t_series_idx).(condition{cond}).Alexa594.FrameMaxCurvature);
                    elseif strcmp(method_AgonistEntry, 'Threshold')
                        mydata(t_series_idx).(condition{cond}).uncagingframe = floor(mydata(t_series_idx).(condition{cond}).Alexa594.FirstFrameAboveBL);
                    end
                end
            end
            
            clear f curr_files curr_dir r uncaging_recording interpulse_duration_ms num_pulses pulseduration_ms uncaging_start_ms uncaging_end_ms ...
                uncaging_duration_ms uncaging_frame_start uncaging_frame_end max_voltage num_bursts
        end
    end
    clear t_series_idx
end

cd(save_dir)

save('dirs_2024XXXX.mat', 'dirs'); %save the pathnames of all included AQuA recordings
fprintf ('Done!')

cd(start_dir);
clear dirs fields cond