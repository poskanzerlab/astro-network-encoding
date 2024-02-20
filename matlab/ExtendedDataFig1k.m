%% Extended Data Fig. 1 k: ROI-based analysis of fluorescence traces
% pinkFlamindo (cAMP) and GCaMP (Ca2+) responses to bath application of receptor agonist
% Peak counting: comparing between cAMP and Ca2+ datasets
%
% 1. Load workspace 'FIJI_AQuA_CytoGCaMP_BathAppBacLY_ExtFig1k.mat' and
%    run Section A & B with indicator = 'CytoGCaMP'
% 2. Load workspace 'FIJI_pinkFlamindo_ReceptorAgonistBathApp_ExtFig1.mat' and
%    re-run Section A & B with indicator = 'pinkFlamindo'
% 3. Run section C after running A&B on both datasets to generate plots in Ext. Data Fig. 1k
%
% Michelle Cahill 20231207
%% Section A: From PrepPinkFlamindStructForFittingSigmoids_MultCond_MC20211208.m
% Calcuate dF/F for each defined region from raw fluorescence traces
indicator = 'pinkFlamindo'; %'pinkFlamindo', 'CytoGCaMP'
condition = fields(mydata);

% 1. Approximate frame of agonist entry (uncaging frame)
%   average frame of peak light artifact + 90 frames (based on 20200127
%   data)
if strcmp(indicator, 'pinkFlamindo')
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
end

% 2. Remove light artifact from traces 20 frames before and 10 frames after
if strcmp(indicator, 'pinkFlamindo')
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
end

if strcmp(indicator, 'pinkFlamindo')
    TRACE = tLARemoved; %For pinkFlamindo, use the raw trace minus background for dF/F calculations
elseif strcmp(indicator, 'CytoGCaMP')
    TRACE = 'RawTrace'; %For CytoGCaMP, use the raw trace for dF/F calculations
end

% 3. Calculate dF/F per region for a fixed baseline (used for ext. data
% fig. 1j)
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if strcmp(indicator, 'pinkFlamindo')
            RawTrace = mydata(tseries).(condition{cond}).(TRACE);
        elseif strcmp(indicator, 'CytoGCaMP')
            RawTrace = mydata(tseries).(condition{cond}).ImageJByRegion.(TRACE);
        end
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
        if strcmp(indicator, 'pinkFlamindo')
            RawTrace = mydata(tseries).(condition{cond}).(TRACE);
        elseif strcmp(indicator, 'CytoGCaMP')
            RawTrace = mydata(tseries).(condition{cond}).ImageJByRegion.(TRACE);
        end
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
            RawTrace = mydata(tseries).(condition{cond}).(TRACE);
        elseif strcmp(indicator, 'CytoGCaMP')
            RawTrace = mydata(tseries).(condition{cond}).ImageJByRegion.(TRACE);
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
clear cond LowerXPercent dFF_T tseries t2remove tLARemoved Trace2Edit TRACE
%% Section B: Counting peaks
% From PeaksPerMinByActiveRegion_MultCond_MC20220126.m
% Detect peaks xSTD above baseline mean & calculate peaks per min pre vs post-stim
% plot pre vs post by region or by t-series
% Michelle Cahill 20220126
%
% For Ext. Data Fig. 1k: Trace2Plot = 'ImageJByRegion'; dFF_T = 'dFF_lower_100Per_10s';
%         ERRORBARS = 'sem'; STD_THRESH = 3; LowPassFilt = 0; AVERAGE_ACROSS = 'TSERIES';
%         TreatmentVals = 100;

% 1. Calculating peaks per minute pre vs post by active region from ImageJ Traces
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1k'; 

Trace2Plot = 'ImageJByRegion'; %'dff', 'ActivePx', 'EventCountPerFrame', 'EventCountPerFrame_NormAllEvts', 'ImageJTraces'
dFF_T = 'dFF_lower_100Per_10s';
ERRORBARS = 'sem';
STD_THRESH = 3; %Number of standard deviations above baseline for detecting peaks
LowPassFilt = 0; %1 to use the lowpass filter; 0 to not
AVERAGE_ACROSS = 'TSERIES'; % 'TSERIES', 'REGIONS'
TreatmentVals = 100;
indicator = 'pinkFlamindo'; % 'pinkFlamindo', 'CytoGCaMP'

if ~exist('Change_PeakPerMin')
    Change_PeakPerMin = struct();
end
if ~exist('PeakPerMinPreVPost')
    PeakPerMinPreVPost = struct();
end
if ~exist('TimeToFirstPeak')
    TimeToFirstPeak = struct();
end

if strcmp(indicator, 'CytoGCaMP') %For CytoGCaMP, only use regions with >= 1 AQuA detected event during the recording
    for cond = 1:length(condition)
        for TV = TreatmentVals
            concentration_TS_idx = sprintf('%s_%s%d',indicator, condition{cond}, TV);
            
            %         temp_TS_idx = pulse_x_idx.(concentration_TS_idx);
            %         temp_TS_idx = pulse_x_idx.WT_100;
            temp_TS_idx = 1:length(mydata);
            
            temp_preVpost = [];
            temp_TimeToFirstPeak = [];
            temp_Change = [];
            for tseries = temp_TS_idx
                if ~isempty(mydata(tseries).(condition{cond}))  && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
                    % Finding cells with at lest 1 AQuA event within across all
                    % conditions in that FOV
                    %                     CellsWithActivityCond1 = sum(mydata(tseries).(condition{1}).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                    %                     CellsWithActivityCond2 = sum(mydata(tseries).(condition{2}).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                    %                     CellsWithActivity = sum([CellsWithActivityCond1; CellsWithActivityCond2], 1);
                    % Finding cells with at lest 1 AQuA event within this
                    % condition
                    CellsWithActivity = sum(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                    ACTIVE_REGIONS = find(CellsWithActivity);
                    PerFieldActive = mydata(tseries).(condition{cond}).(Trace2Plot).(dFF_T)(:,ACTIVE_REGIONS);
                    clear CellsWithActivity ACTIVE_REGIONS CellsWithActivityCond1 CellsWithActivityCond2
                    
                    FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe;
                    SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
                    MinPre = ((FrameOfIncrease -1) * SecPerFrame) / 60;
                    MinPost = ((size(PerFieldActive, 1) - FrameOfIncrease) * SecPerFrame) / 60;
                    ts_preVpost = [];
                    ts_TimeToFirstPeak = [];
                    for region = 1:size(PerFieldActive, 2)
                        if LowPassFilt == 1
                            %low-pass filter fluorescence signal
                            lowpassfilt = designfilt('lowpassfir', 'FilterOrder', 20, 'CutoffFrequency', 0.1, 'SampleRate', 1/SecPerFrame); % 'CutoffFrequency', .05
                            trace2plot = filtfilt(lowpassfilt, PerFieldActive(:, region));
                            peakthresh = 1;
                        else
                            trace2plot = PerFieldActive(:, region);
                            peakthresh = 3;
                        end
                        peakthresh = 0;
                        MinPeakWidth = 5;
                        
                        [peak_frames, peak_dFFMax, peak_latency] = Peaks_3STD(trace2plot, FrameOfIncrease, SecPerFrame, STD_THRESH, peakthresh, MinPeakWidth);
                        
                        cd(save_dir)
                        saveas(gcf, sprintf('%s_%s%duMPeaksTS%dReg%d_ImageJFilt%d_Peaks%dSTD_peakthresh%d_MinPeakWidth%d.tif', indicator, condition{cond}, TV, tseries, region, LowPassFilt, STD_THRESH, peakthresh, MinPeakWidth))
                        close
                        cd(start_dir)
                        
                        NumPeaksPre = numel(peak_frames{1});
                        NumPeaksPost = numel(peak_frames{2});
                        PeakPerMinPre = NumPeaksPre / MinPre;
                        PeakPerMinPost = NumPeaksPost / MinPost;
                        ts_preVpost = [ts_preVpost; PeakPerMinPre, PeakPerMinPost];
                        if ~isempty(peak_latency{2})
                            ts_TimeToFirstPeak = [ts_TimeToFirstPeak; min(peak_latency{2})];
                        else
                            ts_TimeToFirstPeak = [ts_TimeToFirstPeak; nan];
                        end
                        
                        clear peak_frames NumPeaksPre NumPeaksPost PeakPerMinPre PeakPerMinPost...
                            peak_dFFMax peak_latency trace2plot lowpassfilt peakthresh MinPeakWidth
                    end
                    clear region MinPre MinPost SecPerFrame FrameOfIncrease PerFieldActive
                    
                    ts_change = ts_preVpost(:,2) - ts_preVpost(:,1);
                    
                    if strcmp(AVERAGE_ACROSS, 'REGIONS')
                        temp_preVpost = [temp_preVpost; ts_preVpost];
                        temp_TimeToFirstPeak = [temp_TimeToFirstPeak; ts_TimeToFirstPeak];
                        temp_Change = [temp_Change; ts_change];
                    elseif strcmp(AVERAGE_ACROSS, 'TSERIES')
                        TS_preVpostMean = mean(ts_preVpost, 1, 'omitnan');
                        temp_preVpost = [temp_preVpost; TS_preVpostMean];
                        ts_TimeToFirstPeakMean = mean(ts_TimeToFirstPeak, 1, 'omitnan');
                        temp_TimeToFirstPeak = [temp_TimeToFirstPeak; ts_TimeToFirstPeakMean];
                        ts_changeMean = mean(ts_change, 1, 'omitnan');
                        temp_Change = [temp_Change; ts_changeMean];
                    end
                    clear TS_preVpostMean ts_preVpost ts_TimeToFirstPeakMean...
                        ts_TimeToFirstPeak ts_changeMean ts_change
                end
            end
            PeakPerMinPreVPost.(concentration_TS_idx) = temp_preVpost;
            TimeToFirstPeak.(concentration_TS_idx) = temp_TimeToFirstPeak;
            Change_PeakPerMin.(concentration_TS_idx) = temp_Change;
            clear concentration_TS_idx temp_TS_idx tseries temp_preVpost temp_TimeToFirstPeak temp_Change
        end
    end
    clear TV cond
elseif strcmp(indicator, 'pinkFlamindo') %For pinkFlamindo, use all identified regions
    for cond = 1:length(condition)
        for TV = TreatmentVals
            concentration_TS_idx = sprintf('%s_%s%d',indicator, condition{cond}, TV);
            
            %         temp_TS_idx = pulse_x_idx.(concentration_TS_idx);
            %         temp_TS_idx = pulse_x_idx.WT_100;
            temp_TS_idx = 1:length(mydata);
            
            temp_preVpost = [];
            temp_TimeToFirstPeak = [];
            temp_Change = [];
            for tseries = temp_TS_idx
                if ~isempty(mydata(tseries).(condition{cond}))
                    PerFieldActive = mydata(tseries).(condition{cond}).(Trace2Plot).(dFF_T);
                    
                    FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe;
                    SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
                    MinPre = ((FrameOfIncrease -1) * SecPerFrame) / 60;
                    MinPost = ((size(PerFieldActive, 1) - FrameOfIncrease) * SecPerFrame) / 60;
                    ts_preVpost = [];
                    ts_TimeToFirstPeak = [];
                    for region = 1:size(PerFieldActive, 2)
                        if LowPassFilt == 1
                            %low-pass filter fluorescence signal
                            lowpassfilt = designfilt('lowpassfir', 'FilterOrder', 20, 'CutoffFrequency', 0.1, 'SampleRate', 1/SecPerFrame); % 'CutoffFrequency', .05
                            trace2plot = filtfilt(lowpassfilt, PerFieldActive(:, region));
                            peakthresh = 1;
                        else
                            trace2plot = PerFieldActive(:, region);
                            peakthresh = 3;
                        end
                        peakthresh = 0;
                        MinPeakWidth = 5;
                        
                        [peak_frames, peak_dFFMax, peak_latency] = Peaks_3STD(trace2plot, FrameOfIncrease, SecPerFrame, STD_THRESH, peakthresh, MinPeakWidth);
                        
                        cd(save_dir)
                        saveas(gcf, sprintf('%s_%s%duMPeaksTS%dReg%d_ImageJFilt%d_Peaks%dSTD_peakthresh%d_MinPeakWidth%d.tif', indicator, condition{cond}, TV, tseries, region, LowPassFilt, STD_THRESH, peakthresh, MinPeakWidth))
                        close
                        cd(start_dir)
                        
                        NumPeaksPre = numel(peak_frames{1});
                        NumPeaksPost = numel(peak_frames{2});
                        PeakPerMinPre = NumPeaksPre / MinPre;
                        PeakPerMinPost = NumPeaksPost / MinPost;
                        ts_preVpost = [ts_preVpost; PeakPerMinPre, PeakPerMinPost];
                        if ~isempty(peak_latency{2})
                            ts_TimeToFirstPeak = [ts_TimeToFirstPeak; min(peak_latency{2})];
                        else
                            ts_TimeToFirstPeak = [ts_TimeToFirstPeak; nan];
                        end
                        
                        clear peak_frames NumPeaksPre NumPeaksPost PeakPerMinPre PeakPerMinPost...
                            peak_dFFMax peak_latency trace2plot lowpassfilt peakthresh MinPeakWidth
                    end
                    clear region MinPre MinPost SecPerFrame FrameOfIncrease PerFieldActive
                    
                    ts_change = ts_preVpost(:,2) - ts_preVpost(:,1);
                    
                    if strcmp(AVERAGE_ACROSS, 'REGIONS')
                        temp_preVpost = [temp_preVpost; ts_preVpost];
                        temp_TimeToFirstPeak = [temp_TimeToFirstPeak; ts_TimeToFirstPeak];
                        temp_Change = [temp_Change; ts_change];
                    elseif strcmp(AVERAGE_ACROSS, 'TSERIES')
                        TS_preVpostMean = mean(ts_preVpost, 1, 'omitnan');
                        temp_preVpost = [temp_preVpost; TS_preVpostMean];
                        ts_TimeToFirstPeakMean = mean(ts_TimeToFirstPeak, 1, 'omitnan');
                        temp_TimeToFirstPeak = [temp_TimeToFirstPeak; ts_TimeToFirstPeakMean];
                        ts_changeMean = mean(ts_change, 1, 'omitnan');
                        temp_Change = [temp_Change; ts_changeMean];
                    end
                    clear TS_preVpostMean ts_preVpost ts_TimeToFirstPeakMean...
                        ts_TimeToFirstPeak ts_changeMean ts_change
                end
            end
            PeakPerMinPreVPost.(concentration_TS_idx) = temp_preVpost;
            TimeToFirstPeak.(concentration_TS_idx) = temp_TimeToFirstPeak;
            Change_PeakPerMin.(concentration_TS_idx) = temp_Change;
            clear concentration_TS_idx temp_TS_idx tseries temp_preVpost temp_TimeToFirstPeak temp_Change
        end
    end
    clear TV cond
end

clear start_dir save_dir Trace2Plot dFF_T indicator mydata condition
%% After section A & B have been run on both Cyto-GCaMP and pinkFlamindo datasets, run the following
% Section C: plots for Extended Data Fig. 1k
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1k';

% Re-order fieldnames in Change_PeakPerMin and PeakPerMinPreVPost
% (Cyto-Baclofen, pinkFlamindo-Baclofen, Cyto-LY, pinkFlamindo-LY)
FIELDS = fieldnames(Change_PeakPerMin);
ORDER_OF_FIELDS = {FIELDS{1}; FIELDS{3}; FIELDS{2}; FIELDS{4}};
condition = ORDER_OF_FIELDS;
clear FIELDS

% Plot peaks per minute pre- vs post-stim by region or t-series, all
% conditions on the same axes (Extended Data Fig. 1k, top)
for TV = TreatmentVals
%     Param_PreVPost = struct();
%     for cond = 1:length(condition)
%         concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
%         Param_PreVPost.(concentration_TS_idx) = PeakPerMinPreVPost.(concentration_TS_idx);
%     end
    Param_PreVPost = orderfields(PeakPerMinPreVPost, ORDER_OF_FIELDS);
    
    %plotting each t-series mean pre and mean post and the group mean all
    %conditions of this treatment value on the same axes
    CONDITION = fieldnames(Param_PreVPost);
    TITLE = {sprintf('Peak frequency: %duM low pass filt (%d)', TV, LowPassFilt),...
        sprintf('peaks > baseline mean + %dSTD, by %s', STD_THRESH, lower(AVERAGE_ACROSS))};
    COLORS = [0.75 0.75 0.75;... %light grey for pre-stim
        0.4660 0.6740 0.1880]; % green for post-stim
    YLABEL = 'peaks/min';
    XLABEL = 'conditions';
    
    BarPlotGrouped_PreVPostTSMeanERROR_IndTrials(Param_PreVPost, CONDITION, COLORS, TITLE, YLABEL, XLABEL, ERRORBARS)
    
    cd(save_dir)
    saveas(gcf, sprintf('AllCond%duMPkPerMinBy%s_Peaks%dSTD_LPFILT%d.tif', TV, lower(AVERAGE_ACROSS), STD_THRESH, LowPassFilt))
    saveas(gcf, sprintf('AllCond%duMPkPerMinBy%s_Peaks%dSTD_LPFILT%d.svg', TV, lower(AVERAGE_ACROSS), STD_THRESH, LowPassFilt))
    close
    cd(start_dir)
    
    clear Param_PreVPost cond CONDITION TITLE COLORS YLABEL XLABEL YLABEL NUM_TS 
end
clear TV

% 3. Plot change in peak frequency (post-pre) comparing across all conditions
for TV = TreatmentVals
    %     condition = sort(condition);
    %     if ~isempty(COLOR_Scheme{cp})
    %         cmap = cbrewer(COLOR_Scheme{cp}, COLORS{cp}, 8,'PCHIP');
    %         if contains(condition2plot{cp}, '_')
    %             cmap(1,:) = [0.5 0.5 0.5];
    %         end
    %         if contains(condition2plot{cp}, '_') || strcmp(COLOR_Scheme{cp}, 'div')
    %             cmap_num = flip(floor(linspace(1, 8, length(condition))));
    %         else
    %             cmap_num = flip(floor(linspace(3, 8, length(condition))));
    %         end
    %     else
    if length(condition) == 4
        cmap = [1 0 1; 1 0 1; 0 1 1; 0 1 1]; %magneta and cyan
        cmap_num = [1, 2, 3, 4];
    elseif length(condition) == 2
        cmap = [1 0 1;  0 1 1]; %magneta and cyan
        cmap_num = [1, 2];
    end
    %     end
    
    %     latency_dist = struct();
    %     for cond = 1:length(condition)
    %         concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
    %         temp_latency = Change_PeakPerMin.(concentration_TS_idx);
    %         dc_temp_latency = temp_latency;
    %         latency_dist.(concentration_TS_idx) = dc_temp_latency;
    %
    %         clear temp_latency temp_IncDecNoChange dc_temp_latency
    %     end
    %     clear cond
    Param_PreVPost = orderfields(Change_PeakPerMin, ORDER_OF_FIELDS);
    
%     concentration_TS_idx = sprintf('%s%d',condition{1}, TV);
    TITLE = {sprintf('Change in peak frequency ImageJTrace: %duM low pass filt (%d)', TV, LowPassFilt),...
        sprintf('peaks > baseline mean + %dSTD, by %s', STD_THRESH, lower(AVERAGE_ACROSS))};
    YLABEL = 'change in peak frequency (post-pre)';
    
    % mean +- errorbars & individual datapoints offset
    if length(condition) == 2
        PlotMean_IndDataPts_MultCond( Param_PreVPost, TITLE, YLABEL, cmap, cmap_num,...
            ERRORBARS, 'RankSum')
    else
        PlotMean_IndDataPts_MultCond( Param_PreVPost, TITLE, YLABEL, cmap, cmap_num,...
            ERRORBARS, 'ANOVA')
    end
    cd(save_dir)
    saveas(gcf, sprintf('AllConditions_ChangeInPeakFrequencyBy%s_%dSTD_MeanIndPts_LowPassFilt%d.tif', lower(AVERAGE_ACROSS), STD_THRESH, LowPassFilt))
    saveas(gcf, sprintf('AllConditions_ChangeInPeakFrequencyBy%s_%dSTD_MeanIndPts_LowPassFilt%d.svg', lower(AVERAGE_ACROSS), STD_THRESH, LowPassFilt))
    close
    cd(start_dir)
    
    % violin plot
    if length(condition) == 2
        violinplot_STATS(Param_PreVPost, [], TITLE, YLABEL, cmap, cmap_num, 'RankSum')
    else
        violinplot_STATS(Param_PreVPost, [], TITLE, YLABEL, cmap, cmap_num, 'ANOVA')
    end
    cd(save_dir)
    saveas(gcf, sprintf('AllConditions_ChangeInPeakFrequencyBy%s_%dSTD_Violins_LowPassFilt%d.tif', lower(AVERAGE_ACROSS), STD_THRESH, LowPassFilt))
    saveas(gcf, sprintf('AllConditions_ChangeInPeakFrequencyBy%s_%dSTD_Violins_LowPassFilt%d.svg', lower(AVERAGE_ACROSS), STD_THRESH, LowPassFilt))
    close
    cd(start_dir)
    
    RANKSUM_PVAL.Baclofen = ranksum(Param_PreVPost.(condition{1}), Param_PreVPost.(condition{2}));
    RANKSUM_PVAL.LY379268 = ranksum(Param_PreVPost.(condition{3}), Param_PreVPost.(condition{4}));
    
    clear Param_PreVPost YLABEL TITLE concentration_TS_idx Param_PreVPost
end
clear TV

%%
clear start_dir save_dir ERRORBARS STD_THRESH PeakPerMinPreVPost AVERAGE_ACROSS TimeToFirstPeak LowPassFilt...
    Change_PeakPerMin cmap cmap_num condition LowPassFilt ORDER_OF_FIELDS RANDSUM_PVAL...
    TreatmentVals