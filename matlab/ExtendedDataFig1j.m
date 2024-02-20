%% Extended Data Fig. 1 j: ROI-based analysis of fluorescence traces
%  pinkFlamindo (cAMP) responses to bath application of receptor agonist
%
% 1. Load workspace 'FIJI_pinkFlamindo_ReceptorAgonistBathApp_ExtFig1.mat'
% 2. Run the following code
%
% Michelle Cahill 20231206
%% Section A: From PrepPinkFlamindStructForFittingSigmoids_MultCond_MC20211208.m
% Calcuate dF/F for each defined region from raw fluorescence traces

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1j_methIncDecMaxCurve_methAgEntThresh_TraLenEntire_LPFilt_dFFFixedBL__ComOut291'; %
indicator = 'pinkFlamindo'; %'pinkFlamindo', 'CytoGCaMP'
condition = fields(mydata);

% 1. Approximate frame of agonist entry (uncaging frame)
%   average frame of peak light artifact + 90 frames (based on 20200127
%   data)
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        LA_Channel = mydata(tseries).(condition{cond}).RawFluorescence;
        maxPeakFrame = NaN(1,size(LA_Channel, 2));
        for region = 1:size(LA_Channel, 2)
            [AllPeaksHeight, AllPeaksFrame] = findpeaks(LA_Channel(:,region));
            MaxPeakHeight = max(AllPeaksHeight);
            maxPeakFrame(1,region) = AllPeaksFrame(AllPeaksHeight == MaxPeakHeight);
            clear AllPeaksHeight AllPeaksFrame MaxPeakHeight
        end
        clear region
        PeakLAFrame = nanmean(maxPeakFrame);
        mydata(tseries).(condition{cond}).uncagingframe = round(PeakLAFrame) + 90;
        clear LA_Channel maxPeakFrame PeakLAFrame
    end
end
clear cond tseries

% 2. Remove light artifact from traces 20 frames before and 10 frames after
t2remove = 'RawMinusBackground';
tLARemoved = sprintf('%s_LARemoved', t2remove);
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        LA_Channel = mydata(tseries).(condition{cond}).RawFluorescence;
        Trace2Edit = mydata(tseries).(condition{cond}).(t2remove);
        LARemoved = NaN(size(Trace2Edit));
        for region = 1:size(LA_Channel, 2)
            [AllPeaksHeight, AllPeaksFrame] = findpeaks(LA_Channel(:,region));
            MaxPeakHeight = max(AllPeaksHeight);
            maxPeakFrame = AllPeaksFrame(AllPeaksHeight == MaxPeakHeight);
            start_LA = maxPeakFrame - 20;
            end_LA = maxPeakFrame + 10;
            LARemoved(:,region) = Trace2Edit(:,region);
            LARemoved(start_LA:end_LA, region) = mean(Trace2Edit(start_LA-30:start_LA-1,region));
            clear AllPeaksHeight AllPeaksFrame MaxPeakHeight start_LA end_LA
        end
        clear region
        mydata(tseries).(condition{cond}).(tLARemoved) = LARemoved;
        clear LA_Channel maxPeakFrame PeakLAFrame LARemoved
    end
end
clear cond tseries

% 3. Calculate dF/F per region for a fixed baseline (used for ext. data
% fig. 1j)
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        RawTrace = mydata(tseries).(condition{cond}).(tLARemoved);
        SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
        BaselinePeriod = [1,100];
        dFF_T = sprintf('dFF_FixedBL%d_%d', BaselinePeriod(1), BaselinePeriod(2));
        mydata(tseries).(condition{cond}).ImageJByRegion.(dFF_T) =  dFF_FixedBaseline_20211208(RawTrace, BaselinePeriod);
        
%         figure()
%         plot(mydata(tseries).(condition{cond}).ImageJByRegion.(dFF_T))
%         hold on
%         vline(mydata(tseries).(condition{cond}).uncagingframe, 'r')
%         title({sprintf('ImageJ Traces By Region: %s TS %d', condition{cond}, tseries),...
%             sprintf('dF/F fixed baseline: frame %d-%d',BaselinePeriod(1), BaselinePeriod(2))})
%         xlabel('time (frames)')
%         ylabel('dF/F')
%         
%         cd(save_dir)
%         saveas(gcf, sprintf('ImageJTracesByRegion_dFF_CytoG_%s_TS%d.tif', condition{cond}, tseries))
%         saveas(gcf, sprintf('ImageJTracesByRegion_dFF_CytoG_%s_TS%d.svg', condition{cond}, tseries))
%         close
%         cd(start_dir)
        clear RawTrace SecPerFrame BaselinePeriod
    end
end
clear cond tseries 

% 4. Calculate dF/F per region using a moving window and the lower 50% of
% frames
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        RawTrace = mydata(tseries).(condition{cond}).(tLARemoved);
        SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
        MovingWindow = 10;
        dFF_T = sprintf('dFF_lower%ds', MovingWindow);
        mydata(tseries).(condition{cond}).ImageJByRegion.(dFF_T) = dFF_lower_20211202(RawTrace, SecPerFrame, MovingWindow);
        
%         figure()
%         plot(mydata(tseries).(condition{cond}).ImageJByRegion.(dFF_T))
%         hold on
%         vline(mydata(tseries).(condition{cond}).uncagingframe, 'r')
%         title({sprintf('ImageJ Traces By Region: %s TS %d', condition{cond}, tseries),...
%             sprintf('dF/F lower: %dS moving window', MovingWindow)})
%         xlabel('time (frames)')
%         ylabel('dF/F lower')
%         
%         cd(save_dir)
%         saveas(gcf, sprintf('ImageJTracesByRegion_dFF%dSWind_%s_%s_TS%d.tif', MovingWindow, indicator, condition{cond}, tseries))
%         saveas(gcf, sprintf('ImageJTracesByRegion_dFF%dSWind_%s_%s_TS%d.svg', MovingWindow, indicator, condition{cond}, tseries))
%         close
%         cd(start_dir)
        clear RawTrace SecPerFrame MovingWindow
    end
end

% 5. Calculate dF/F per region using a moving window and the lower X% of
% frames (used for ext. data fig. 1k)
LowerXPercent = 100;
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if strcmp(indicator, 'pinkFlamindo')
            RawTrace = mydata(tseries).(condition{cond}).(tLARemoved); %For pinkFlamindo
        elseif strcmp(indicator, 'CytoGCaMP')
            RawTrace = mydata(tseries).(condition{cond}).ImageJByRegion.RawTrace;
        end
        SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
        MovingWindow = 10;
        dFF_T = sprintf('dFF_lower_%dPer_%ds', LowerXPercent, MovingWindow);
        mydata(tseries).(condition{cond}).ImageJByRegion.(dFF_T) = dFF_lowerXPercent_20220127(RawTrace, SecPerFrame, MovingWindow, LowerXPercent);
        
%         figure()
%         plot(mydata(tseries).(condition{cond}).ImageJByRegion.(dFF_T))
%         hold on
%         vline(mydata(tseries).(condition{cond}).uncagingframe, 'r')
%         title({sprintf('ImageJ Traces By Region: %s TS %d', condition{cond}, tseries),...
%             sprintf('dF/F lower %d percent: %dS moving window', LowerXPercent, MovingWindow)})
%         xlabel('time (frames)')
%         ylabel('dF/F lower')
%         
%         cd(save_dir)
%         saveas(gcf, sprintf('ImageJTracesByRegion_dFFLower%dPer_%dSWind_%s_%s_TS%d.tif', LowerXPercent, MovingWindow, indicator, condition{cond}, tseries))
%         saveas(gcf, sprintf('ImageJTracesByRegion_dFFLower%dPer_%dSWind_%s_%s_TS%d.svg', LowerXPercent, MovingWindow, indicator, condition{cond}, tseries))
%         close
%         cd(start_dir)
        clear RawTrace SecPerFrame MovingWindow
    end
end
clear cond LowerXPercent dFF_T tseries t2remove tLARemoved Trace2Edit 
%% Section B: Extended Data Figure 1j
% From FittingSigmoids_ImageJTracesByRegion_MultCond_MC20211202.m
% 1. Fit sigmoids to fluorescence traces from each region & calculate frame of response
%    as frame of max curvature, when slope is above certain threshold and
%    whether activity increased, decreased or stayed the same
%    For each condition, record the average latency to response for all
%    cells in a t-series that show an increase (Latency) 
%    Record the response direction of each cell (inc, dec or no change)
%    a. Choose how to determine a response (method_IncDec):
%       -Fitting sigmoid and finding frame of max curvature ('MaxCurve')
%       -Finding first frame above or below baseline thresh ('Threshold')
%    b. Choose raw or dF/F trace (TRACE_TYPE)
%    c. Choose low-pass filter or no filter (LowPassFilt)
% 2. Plot average traces for each condition +- errorbars of choice
% 3. Plot pie charts indicating percent of cells in each condition showing
%    increase or no change following agonist
%
% For Ext. Data Fig. 1j: method_IncDec = 'MaxCurve'; method_AgonistEntry = 'Threshold';
%        ExcludeBaseline = 0; ExcludeRampingCells = 0; TRACE = 'ImageJByRegion';
%        TRACE_TYPE = 'dFF'; dFF_T = 'dFF_FixedBL1_100'; LowPassFilt = 1;

method_IncDec = 'MaxCurve'; %'Threshold', 'MaxCurve'
method_AgonistEntry = 'Threshold'; %'Threshold', 'MaxCurve'

ExcludeBaseline = 0; %first x seconds to exclude
ExcludeRampingCells = 0; %1 to exclude the list of ramping cells, 0 to include all cells
TRACE = 'ImageJByRegion';
TRACE_TYPE = 'dFF'; %'RawTrace', 'dFF'
dFF_T = 'dFF_FixedBL1_100'; % Choose the dF/F normalization 'dFF_FixedBL1_100'
LowPassFilt = 1; %1 to use the lowpass filter; 0 to not
if strcmp(method_IncDec, 'MaxCurve')
    TraceLength = 'EntireTrace'; %'FirstPeak' or 'EntireTrace'
end

Latency = struct(); %average latency of cells in each t-series showing an increase
IncDecNoChange = struct();
TracesAligned = struct();
condition = fieldnames(mydata);

if strcmp(TRACE_TYPE, 'dFF')
    TRACE_TYPE = dFF_T;
end
if LowPassFilt == 1
    FILT = 'LPFilt';
elseif LowPassFilt == 0
    FILT = 'NoFilt';
end

% 1. Fit sigmoids to fluorescence traces from each to determine category of
% activity change
for cond = 1:length(condition)
    temp_latency = [];
    temp_IncDecNoChange = [];
    mean_AllTraces = [];
    mean_IncTraces = [];
    
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            clear fluorescence_trace TS_latency TS_IncDecNoChange
            TS_latency = [];
            TS_IncDecNoChange = [];
            TS_TRACES = [];
            if strcmp(indicator, 'CytoGCaMP')
                if strcmp(method_AgonistEntry, 'MaxCurve')
                    frame_agonistEntry = floor(mydata(tseries).(condition{cond}).Alexa594.FrameMaxCurvature);
                elseif strcmp(method_AgonistEntry, 'Threshold')
                    frame_agonistEntry = mydata(tseries).(condition{cond}).Alexa594.FirstFrameAboveBL;
                end
            elseif strcmp(indicator, 'pinkFlamindo')
                %                 if strcmp(method_AgonistEntry, 'MaxCurve')
                %                     frame_agonistEntry = mydata(tseries).(condition{cond}).EndLightArtifact + 175; %Based on frames between light artifact and max curve of Alexa 488 20200127_Baclofen1
                %                 elseif strcmp(method_AgonistEntry, 'Threshold')
                %                     frame_agonistEntry = mydata(tseries).(condition{cond}).EndLightArtifact + 90; %Based on frames between light artifact and first frame 3std above BL mean of Alexa 488 20200127_Baclofen1
                %                 end
                frame_agonistEntry = mydata(tseries).(condition{cond}).uncagingframe;
            end
%             frame_agonistEntry = mydata(tseries).(condition{cond}).uncagingframe;
            SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
            fluorescence_trace = mydata(tseries).(condition{cond}).(TRACE).(TRACE_TYPE);
            
            XVals = (1:1:size(fluorescence_trace,1));
            ExcludeBaselineFrames = round(ExcludeBaseline / mydata(tseries).(condition{cond}).SecPerFrame);
            %             EndFrame = size(fluorescence_trace, 2);
            %             EndFrame = 300;
            if ExcludeBaselineFrames < 0
                BaselinePeriod = [ExcludeBaselineFrames, frame_agonistEntry-10];
            else
                BaselinePeriod = [1, frame_agonistEntry-10];
            end
            
            %low-pass filter fluorescence signal
            lowpassfilt = designfilt('lowpassfir', 'FilterOrder', 20, 'CutoffFrequency', 0.01, 'SampleRate', 1/SecPerFrame); % 'CutoffFrequency', .05
            
            if ExcludeRampingCells
                if isfield(mydata(tseries).(condition{cond}), 'RampingCellIdx')
                    ramping_cell_idx = mydata(tseries).(condition{cond}).RampingCellIdx;
                else
                    ramping_cell_idx = [];
                end
            else
                ramping_cell_idx = [];
            end
            
            Nonzero_cells = 1:size(fluorescence_trace,2);
            Nonzero_cells = Nonzero_cells(~ismember(Nonzero_cells, ramping_cell_idx));
            
            for region = Nonzero_cells
                clear fluorescenceFilt trace2plot
                if LowPassFilt == 1
                    fluorescenceFilt = filtfilt(lowpassfilt, fluorescence_trace(:, region));
                    trace2plot = fluorescenceFilt;
                elseif LowPassFilt == 0
                    trace2plot = fluorescence_trace(:, region);
                end
                
                if strcmp(method_IncDec, 'MaxCurve')
                    if strcmp(TraceLength, 'FirstPeak')
                        PostStim_PeakFrame = Peak3STD(trace2plot, frame_agonistEntry, BaselinePeriod);
                        if ~isnan(PostStim_PeakFrame)
                            if (PostStim_PeakFrame + 30) <= size(trace2plot, 1)
%                                 EndFrame = PostStim_PeakFrame + 10;
                                EndFrame = PostStim_PeakFrame + 30;
                            else
                                EndFrame = size(trace2plot, 1);
                            end
                        else
                            EndFrame = size(trace2plot, 1);
                        end
                    elseif strcmp(TraceLength, 'EntireTrace')
                        EndFrame = size(trace2plot, 1);
                    end
                    
                    [IncDec, fitresult, gof, SigFitTrace] = createFit_sigmoidDefinedStartEnd(XVals', trace2plot,...
                        BaselinePeriod(1), EndFrame, BaselinePeriod);
                    ylabel(sprintf('%s', TRACE_TYPE))

                    mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).IncDec = IncDec; %categorizing the traces as increase(1), decrease(-1) or no change(0)
                    mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).CoE = coeffvalues(fitresult);
                    mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).GOF = gof;
                    mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).SigFitTrace = SigFitTrace; %y-vals of the fitted sigmoid curve
                    
                    if ~(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).IncDec == 0) &&...
                            mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).CoE(2) > frame_agonistEntry-15 %If the trace changed after the agonist entered the bath (around frame 400), calculate the frame of increase
                        [t1, t2, t3] = PointsMaxCurvature_4CoE_V2(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).CoE); %frame of increase as frame of maximum curvature in the sigmoid curve
                        mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature = min([t1 t2 t3]);
                        if mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature < size(trace2plot, 1) &&...
                                mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature >= mydata(tseries).(condition{cond}).uncagingframe &&...
                                mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature < length(SigFitTrace)
                            mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyMaxCurveSec = (mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature - frame_agonistEntry) * SecPerFrame;
                        else
                            mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyMaxCurveSec = nan;
%                             IncDec = 0;
%                             mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).IncDec = IncDec;
                        end
                        clear t1 t2 t3
                        
                        temp_slope = diff(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).SigFitTrace); %frame of increase as the first frame where the slope of the sigmoid curve is no longer 0
                        if mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).IncDec == 1
                            mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise = find(temp_slope > 0.01, 1)  + ExcludeBaselineFrames;
                            mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyInitRiseSec = (mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise - frame_agonistEntry) * SecPerFrame;
                        elseif mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).IncDec == -1
                            mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise = find(temp_slope < -0.01, 1)  + ExcludeBaselineFrames;
                            mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyInitRiseSec = (mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise - frame_agonistEntry) * SecPerFrame;
                        end
                        clear temp_slope
                    else
                        mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature = nan;
                        mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyMaxCurveSec = nan;
                        mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise = nan;
                        mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyInitRiseSec = nan;
                    end
                    
                    TS_latency = [TS_latency; mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyMaxCurveSec];
                    TS_IncDecNoChange = [TS_IncDecNoChange; IncDec];
                    temp_IncDecNoChange = [temp_IncDecNoChange; IncDec];
                    TS_TRACES = [TS_TRACES, trace2plot];
                    
                    if ismember(region, Nonzero_cells) > 0
                        figure(1)
                        hold on
                        if ~isnan(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature)&&...
                                ~isnan(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyMaxCurveSec)&&...
                                mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature <= EndFrame &&...
                                mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature > 0
                            plot(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature, SigFitTrace(round(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature - ExcludeBaselineFrames)), 'mx', 'MarkerSize', 15)
                            plot(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise, SigFitTrace(round(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise - ExcludeBaselineFrames)), 'cx', 'MarkerSize', 15)
                        end
                        title({sprintf('%s %s: %s TS %d region %d', indicator, TRACE, strrep(condition{cond}, '_', ' '), tseries, region),...
                            sprintf('%s', strrep(TRACE_TYPE, '_', ''))})
                        if LowPassFilt == 1
                            legend('filtered fluorescence trace', 'Fit Sigmoid Curve', 'Max Curve', 'Non-0 slope (>0.01)', 'Location', 'NorthWest' );
                        else
                            legend('raw fluorescence trace', 'Fit Sigmoid Curve', 'Max Curve', 'Non-0 slope (>0.01)', 'Location', 'NorthWest' );
                        end
                        cd(save_dir)
                        saveas(gcf, sprintf('%s_%s_%s_SigFit%sTrace%s_TS%dReg%d.tif', indicator, condition{cond}, TRACE, TRACE_TYPE, FILT, tseries, region))
                        saveas(gcf, sprintf('%s_%s_%s_SigFit%sTrace%s_TS%dReg%d.svg', indicator, condition{cond}, TRACE, TRACE_TYPE, FILT, tseries, region))
                        close
                        cd(start_dir)
                        
                        if LowPassFilt == 1
                            figure()
                            hold on
                            plot(fluorescence_trace(:, region), 'k')
                            plot(BaselinePeriod(1):EndFrame, SigFitTrace, 'r')
                            if ~isnan(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature)&&...
                                    ~isnan(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).LatencyMaxCurveSec)&&...
                                    mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature <= EndFrame &&...
                                mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature > 0
                                plot(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature, SigFitTrace(round(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameMaxCurvature - ExcludeBaselineFrames)), 'mx', 'MarkerSize', 15)
                                plot(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise, SigFitTrace(round(mydata(tseries).(condition{cond}).SigmoidEntireTrace.(TRACE).(FILT).(TRACE_TYPE)(region).FrameInitialRise - ExcludeBaselineFrames)), 'cx', 'MarkerSize', 15)
                            end
                            title(sprintf('%s %s: %s TS %d region %d', indicator, TRACE, strrep(condition{cond}, '_', ' '), tseries, region))
                            legend('raw trace', 'fitted sigmoid', 'Max Curve', 'Non-0 slope (>0.01)','Location', 'Best');
                            cd(save_dir)
                            saveas(gcf, sprintf('%s_%s_%s_%sSigFitNoFilt_TS%dReg%d.tif',indicator, condition{cond}, TRACE, TRACE_TYPE, tseries, region))
                            saveas(gcf, sprintf('%s_%s_%s_%sSigFitNoFilt_TS%dReg%d.svg',indicator, condition{cond}, TRACE, TRACE_TYPE, tseries, region))
                            close
                            cd(start_dir)
                        end
                        
                    end
                    close
                    clear IncDec fitresult gof t1 t2 t3 SigFitTrace fluorescenceFilt trace2plot
                elseif strcmp(method_IncDec, 'Threshold')
                    Z = 3;
                    VARI = 'std';

                    [FirstFrameAboveBelowBL, IncDec] = FirstFrame_AboveBelowBaselineMean(trace2plot, BaselinePeriod, Z, VARI, frame_agonistEntry);
                    mydata(tseries).(condition{cond}).(TRACE).(FILT).(TRACE_TYPE).FirstFrameAboveBelowBL(region) = FirstFrameAboveBelowBL;
                    mydata(tseries).(condition{cond}).(TRACE).(FILT).(TRACE_TYPE).IncDec(region) = IncDec;
                    
                    figure()
                    hold on
                    if LowPassFilt == 1
                        plot(fluorescence_trace(:, region), '-k')
                        plot(trace2plot, '-r')
                    else
                        plot(trace2plot, '-k')
                    end

                    if IncDec ~= 0 && ~isnan(IncDec)
                        plot(FirstFrameAboveBelowBL, trace2plot(round(FirstFrameAboveBelowBL)), 'mx', 'MarkerSize', 15)
                        if LowPassFilt == 1
                            legend('raw fluorescence trace', 'low-pass filtered trace', sprintf('%d%s above/below BL: %d', Z, VARI, IncDec), 'Location', 'NorthWest');
                        else
                            legend('raw fluorescence trace', sprintf('%d%s above/below BL: %d', Z, VARI, IncDec), 'Location', 'NorthWest');
                        end
                    else
                        if LowPassFilt == 1
                            legend('raw fluorescence trace', 'low-pass filtered trace', 'Location', 'NorthWest');
                        else
                            legend('raw fluorescence trace', 'Location', 'NorthWest');
                        end
                    end
                    title({sprintf('%s %s %s:', TRACE, TRACE_TYPE, condition{cond}), sprintf('tseries %d region %d, %s',tseries, region, FILT)})
                    cd(save_dir)
                    saveas(gcf, sprintf('%s%s%s_%d%sAboveBL_%sTS%dReg%d.tif', TRACE, TRACE_TYPE, FILT, Z, VARI, condition{cond}, tseries, region))
                    saveas(gcf, sprintf('%s%s%s_%d%sAboveBL_%sTS%dReg%d.svg', TRACE, TRACE_TYPE, FILT, Z, VARI, condition{cond}, tseries, region))
                    close
                    cd(start_dir)
                    
                    reg_latency = (FirstFrameAboveBelowBL - frame_agonistEntry) * SecPerFrame;
                    mydata(tseries).(condition{cond}).(TRACE).(FILT).(TRACE_TYPE).LatencyAboveBelowThreshSec(region) = reg_latency;
                    
                    TS_latency = [TS_latency; reg_latency];
                    TS_IncDecNoChange = [TS_IncDecNoChange; IncDec];
                    temp_IncDecNoChange = [temp_IncDecNoChange; IncDec];
                    TS_TRACES = [TS_TRACES, trace2plot];
                    clear Z VARI reg_latency fluorescenceFilt trace2plot
                end
            end
            
            mean_AT = nanmean(TS_TRACES,2);
            mean_IT = nanmean(TS_TRACES(:,TS_IncDecNoChange == 1),2);
            mean_AllTraces = [mean_AllTraces; mean_AT'];
            mean_IncTraces = [mean_IncTraces; mean_IT'];
            
            mean_latency = nanmean(TS_latency(TS_IncDecNoChange == 1));
            temp_latency = [temp_latency; mean_latency];
            
            clear frame_agonistEntry SecPerFrame fluorescence_trace XVals region ExcludeBaselineFrames...
                EndFrame BaselinePeriod PostStim_PeakFrame ramping_cell_idx Nonzero_cells...
                lowpassfilt x f TS_latency TS_IncDecNoChange mean_AT mean_IT mean_latency...
                TS_TRACES
        end
    end
    Latency.(condition{cond}) = temp_latency;
    IncDecNoChange.(condition{cond}) = temp_IncDecNoChange;
    
    OnsetFrame = [];
    for tseries = 1:length(mydata)
        OnsetFrame = [OnsetFrame, mydata(tseries).(condition{cond}).uncagingframe];
    end
    
    [TracesAligned.AllTraces.(condition{cond}), AlignedOnsetFrame] = AlignTraces_20211203(mean_AllTraces,OnsetFrame);
    [TracesAligned.IncTraces.(condition{cond}), AlignedOnsetFrame] = AlignTraces_20211203(mean_IncTraces,OnsetFrame);

    clear tseries TSERIES pulse_cond temp_latency temp_IncDecNoChange OnsetFrame...
        mean_AllTraces mean_IncTraces
end
clear cond

% 2. Plot mean+-error bars for multiple conditions on the same graph
%    Extended Data Fig. 1j, right
avg2plot = {'IncTraces'}; %{'AllTraces', 'IncTraces'}
condition2plot = {'AllConditions'}; % Which conditions should be plotted
COLORS = {'MgCy'}; %'RdPu', 'Blues', 'YlGnBu', 'MgCy'
COLOR_Scheme = {''}; %'seq' if COLORS are 'RdPu' or 'Blues'; '' if COLORS are 'MgCy'
ERRORBARS = 'sem'; %'sem' or 'std'
for cp = 1:length(condition2plot)
    %     if strcmp(condition2plot{cp}, 'AllConditions')
    %         condition = fieldnames(mydata);
    %     else
    %         temp_condition = fieldnames(mydata);
    %         if ~contains(condition2plot{cp}, '_')
    %             idx = cell2mat(cellfun(@(x) contains(x, condition2plot{cp}), temp_condition, 'UniformOutput', 0));
    %         else %to pull out conditions with different strings (for example: All RuBiGABA conditions and NoRuBi_LaserUncagingControl would be used with 'RuBiGABA_NoRuBi'
    %             idx_1 = cell2mat(cellfun(@(x) contains(x, extractBefore(condition2plot{cp}, '_')), temp_condition, 'UniformOutput', 0));
    %             idx_2 = cell2mat(cellfun(@(x) contains(x, extractAfter(condition2plot{cp}, '_')), temp_condition, 'UniformOutput', 0));
    %             idx = idx_1 | idx_2;
    %         end
    %         condition = temp_condition(idx, 1);
    %         clear temp_condition idx idx_1 idx_2
    %     end
    
    condition = sort(condition);
    if ~isempty(COLOR_Scheme{cp})
        cmap = cbrewer(COLOR_Scheme{cp}, COLORS{cp}, 8,'PCHIP');
        if contains(condition2plot{cp}, '_')
            cmap(1,:) = [0.5 0.5 0.5];
        end
        if contains(condition2plot{cp}, '_') || strcmp(COLOR_Scheme{cp}, 'div')
            cmap_num = floor(linspace(1, 8, length(condition)));
        else
            cmap_num = floor(linspace(3, 8, length(condition)));
        end
    else
        cmap = [1 0 1; 0 1 1]; %magneta and cyan
        cmap_num = [1, 2];
    end
    
    for TA = 1:length(avg2plot)
        first_bin_All = [];
        last_bin_All = [];
        for cond = 1:length(condition)
            mat2plot = TracesAligned.(avg2plot{TA}).(condition{cond});
            ValBin = find(sum(~isnan(mat2plot),1)>0); %Find only the bins with values, not nan
            first_bin_All = [first_bin_All, ValBin(1)];
            last_bin_All = [last_bin_All, ValBin(end)];
        end
        first_bin = max(first_bin_All);
        last_bin = min(last_bin_All);
        
        % Shifting the Onset Frame to 0 and converting from frames to seconds 
        XVALS = 1:numel(first_bin:last_bin);
        XVALS = (XVALS - numel(first_bin:AlignedOnsetFrame)) * mydata(1).(condition{1}).SecPerFrame;
        
        figure()
        set(gcf,'position', [100, 200, 1500, 600])
        hold on
        
        for cond = 1:length(condition)
            mat2plot = TracesAligned.(avg2plot{TA}).(condition{cond});
%             XLIMS = [1 (last_bin - first_bin+1)];
            mat2plot_cropped = mat2plot(:, first_bin:last_bin);
            
            if strcmp(ERRORBARS, 'std')
                boundedline( XVALS, nanmean(mat2plot_cropped, 1), nanstd(mat2plot_cropped,1), 'alpha', 'cmap', cmap(cmap_num(cond),:))
            elseif strcmp(ERRORBARS, 'sem')
                SEM = nanstd(mat2plot_cropped,1) / (sqrt(size(mat2plot_cropped, 1)-1));
                boundedline( XVALS, nanmean(mat2plot_cropped, 1), SEM, 'alpha', 'cmap', cmap(cmap_num(cond),:))
            end
            
            MEAN_Latency = nanmean(Latency.(condition{cond}));
            [~, temp_frame] = min(abs(XVALS - MEAN_Latency));
            temp_AvgTrace = nanmean(mat2plot_cropped, 1);
            plot(MEAN_Latency, temp_AvgTrace(temp_frame), 'xk', 'MarkerSize',15)
            clear mat2plot mat2plot_cropped Pval_placement MEAN_Latency temp_frame temp_AvgTrace
        end
        
        vline(0, 'r')
%         vline(numel(first_bin:AlignedOnsetFrame), 'r')
%         xlim(XLIMS)
            xlabel('time (s)')
            ylabel(sprintf('mean %s +- %s', strrep(TRACE_TYPE, '_', ' '), ERRORBARS))

        title(sprintf('%s %s %s %s: %s', indicator, TRACE, strrep(TRACE_TYPE, '_', ' '), FILT, avg2plot{TA}))
        DATA_Labels = strrep(condition, '_', ' ');
        DATA_Labels_exp = cell(length(DATA_Labels) *2, 1);
        
        track_temp = 1;
        for d = 1:length(DATA_Labels)
            DATA_Labels_exp{track_temp,1} = '';
            DATA_Labels_exp{track_temp+1,1} = DATA_Labels{d};
            track_temp = track_temp + 2;
        end
        %     end
        legend(DATA_Labels_exp, 'Location', 'best')
        legend('boxoff')
        
        cd(save_dir)
        saveas(gcf, sprintf('%s_%s_%s%s%s_Mean%s_Sec.tif', indicator, TRACE, TRACE_TYPE, FILT, avg2plot{TA}, ERRORBARS))
        saveas(gcf, sprintf('%s_%s_%s%s%s_Mean%s_Sec.svg', indicator, TRACE, TRACE_TYPE, FILT, avg2plot{TA}, ERRORBARS))
        close
        cd(start_dir)
        
        clear first_bin first_bin_All last_bin last_bin_All track_temp DATA_Labels DATA_Labels_exp ValBin XLIMS d SEM XVALS
    end
    clear TA
end
clear cp COLOR_Scheme COLORS condition2plot avg2plot cmap cmap_num cond

% 3. Pie charts of cells that show increase, decrease or no change for each condition
%    Extended Data Fig. 1j, left
condition = fieldnames(IncDecNoChange);
for cond = 1:length(condition)
    temp_IncDecNoChange = IncDecNoChange.(condition{cond});
    PercentInc = (sum(temp_IncDecNoChange == 1))/numel(temp_IncDecNoChange);
    PercentDec = (sum(temp_IncDecNoChange == -1))/numel(temp_IncDecNoChange);
    PercentNoChange = (sum(temp_IncDecNoChange == 0))/numel(temp_IncDecNoChange);
    figure()
    p = pie([PercentNoChange, PercentInc, PercentDec]);
    if PercentNoChange > 0 && PercentInc > 0 && PercentDec > 0
        p(1).FaceColor = [ 0.5 0.5 0.5];
        p(3).FaceColor = 'green';
        p(5).FaceColor = 'red';
        legend('No Change', 'Increase', 'Decrease', 'Location', 'Best')
    elseif PercentNoChange > 0 && PercentInc > 0
        p(1).FaceColor = [ 0.5 0.5 0.5];
        p(3).FaceColor = 'green';
        legend('No Change', 'Increase', 'Location', 'Best')
    elseif PercentNoChange > 0 && PercentDec > 0
        p(1).FaceColor = [ 0.5 0.5 0.5];
        p(3).FaceColor = 'red';
        legend('No Change', 'Decrease', 'Location', 'Best')
    elseif PercentNoChange > 0
        p(1).FaceColor = [ 0.5 0.5 0.5];
        legend('No Change', 'Location', 'Best')
    elseif PercentInc > 0 && PercentDec > 0
        p(1).FaceColor = 'green';
        p(3).FaceColor = 'red';
        legend('Increase', 'Decrease', 'Location', 'Best')
    elseif PercentInc > 0
        p(1).FaceColor = 'green';
        legend('Increase', 'Location', 'Best')
    elseif PercentDec > 0
        p(1).FaceColor = 'red';
        legend('Decrease', 'Location', 'Best')
    end
    title({sprintf('%s %s %s %s: %s', indicator, TRACE, TRACE_TYPE, FILT, strrep(condition{cond}, '_', ' ')),...
        sprintf('n=%d', numel(temp_IncDecNoChange))})
    cd(save_dir)
    saveas(gcf, sprintf('%s_%s_%s%s%s_IncDecNoChange_Pie.tif', indicator, condition{cond}, TRACE, TRACE_TYPE, FILT))
    saveas(gcf, sprintf('%s_%s_%s%s%s_IncDecNoChange_Pie.svg', indicator, condition{cond}, TRACE, TRACE_TYPE, FILT))
    close
    cd(start_dir)
    clear p temp_IncDecNoChange PercentInc PercentDec PercentNoChange
end
clear cond

clear AlignedOnsetFrame dFF_T ERRORBARS ExcludeBaseline ExcludeRampingCells...
    FILT IncDecNoChange Latency LowPassFilt method_IncDec method_AgonistEntry...
    save_dir start_dir TRACE TRACE_TYPE TraceLength TracesAligned