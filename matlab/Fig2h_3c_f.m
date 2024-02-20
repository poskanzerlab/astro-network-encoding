%% Fig. 2h & 3c-f: dF/F raster plots individual recordings in each row
% recordings are sorted by the latency to the first post-stim frame with
% dF/F values > threshold and related panels
% For a flexible number of conditions and the ability to choose which 
% cells to include for mean dF/F calculation for each recording.
% modified from FluorescenceTraceFromAQuAEvtsNoRegion_MC20211117.m
%
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%  3. Run the following code
%
% A: Calculate the mean dF/F from identified events per t-series, 
%    choosing which cells to include in mean dF/F calculation. These
%    calculations are used in the following sections
% B: Fig. 2h & Fig. 3c: Plot raster plots/spectrograms, sorting t-series 
%    in each condition by the time of the first post-stim peak >= threshold. 
%    Choose which cells & tseries to include & what method of normalization to use
% C: Fig. 3c stacked bar graphs showing proportion of recordings 
%    exhibiting an initial fluorescence increase following uncaging
%   (responder). Fisher's exact test compares proportion of responders
%    across WT and Cx43 conditions.
% D: Fig. 3d&e Latency to post-stim frames reaching above threshold for dF/F values
% E. Fig. 3f: Response persistence. Proportion of post-stim frames
%    above threshold following first post-stim frame above threshold. 
%
% Michelle Cahill 20240111
%% Section A: Calculate the mean dF/F from identified events per t-series, 
%  choosing which cells to include in mean dF/F calculation 
%    Currently uses dFF traces from each identified event after removing 
%    the contribution from other events wiht an overlapping footprint: dffMatFilter(:,:,2)
%
%   For Fig. 2h: WE = 'UncagingEvents';
%   For Fig. 3c-f: WE = 'NeighboringEvents';
%   *run both WE = 'UncagingEvents'; and WE = 'NeighboringEvents'; in this
%   section before moving following sections.

condition = fieldnames(mydata);

ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 60; %number of seconds at the beginning of baseline to exclude, 0 to include all baseline period
EBL = sprintf('Exclude%dSBL', ExcludeBaseline);

WE = 'UncagingEvents'; %'NeighboringEvents', 'UncagingEvents', 'AllEvents' which events should be included?

% 1. Calculate the average fluorescence from identified events from chosen regions in
%    each t-series. 
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            clear ActivePxPerFrame_NormAllActPx ActivePxPerFrame_NormAllPxReg
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
            
            ExcludeBaselineFrames = round(ExcludeBaseline / mydata(tseries).(condition{cond}).SecPerFrame);
            FluorescenceTraces_AllRegions = [];
            for r = regions
                evt_idx = [mydata(tseries).(condition{cond}).event_idx_pre{r,1},mydata(tseries).(condition{cond}).event_idx_post{r,1}];
                temp_FluorescenceTraces = NaN(numel(evt_idx), mydata(tseries).(condition{cond}).totalframes);
                e = 1;
                for evt = evt_idx
                    if ~isempty(mydata(tseries).(condition{cond}).tLocs{1,evt}) % && mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(evt) > ExcludeBaselineFrames
                        temp_FluorescenceTraces(e,:) = mydata(tseries).(condition{cond}).res.dffMatFilter(evt, :, 2);
                        e = e+1;
                    end
                end
                FluorescenceTraces_AllRegions = [FluorescenceTraces_AllRegions; temp_FluorescenceTraces(1:e-1, :)];
                clear evt_idx e evt temp_FluorescenceTraces
            end
            clear r
            
            mydata(tseries).(condition{cond}).(EBL).(WE).dffMatFilter2 = FluorescenceTraces_AllRegions; %average fluorescence traces of events from each region from (dffMatFilter(:,:,2))
            
            %Find frames that are >= 3std+baseline mean for avg dF/F across
            %all events
            AVG_TRACE = mean(FluorescenceTraces_AllRegions, 1, 'omitnan');
            if ~isempty(AVG_TRACE)
                BL = AVG_TRACE(ExcludeBaselineFrames+1:mydata(tseries).(condition{cond}).uncagingframe-1);
                BL_Mean = nanmean(BL);
                BL_STD = nanstd(BL);
                PeakThreshold = BL_Mean + (3*BL_STD);
                FramesAboveThresh = AVG_TRACE >= PeakThreshold;
            else
                FramesAboveThresh = false(1, mydata(tseries).(condition{cond}).totalframes);
                BL_Mean = NaN;
                BL_STD = NaN;
                PeakThreshold = NaN;
            end
            mydata(tseries).(condition{cond}).(EBL).(WE).dffMatFilter2_Thresh3std = FramesAboveThresh;
            mydata(tseries).(condition{cond}).(EBL).(WE).Thresh.Threshold = PeakThreshold;
            mydata(tseries).(condition{cond}).(EBL).(WE).Thresh.BL_mean = BL_Mean;
            mydata(tseries).(condition{cond}).(EBL).(WE).Thresh.BL_STD = BL_STD;
            clear BL BL_Mean BL_STD PeakThreshold FramesAboveThresh AVG_TRACE... 
             ax1 ax2 XVALS Traces2Plot TITLE ExcludeBaselineFrames FluorescenceTraces_AllRegions regions ramping_cell_idx
        end
        clear num_evts ActivePxPerFrame evt region region ExcludeBaselineFrames FluorescenceTraces_AllRegions
    end
    clear tseries
end
clear cond
%% Section B. Plot raster plots/spectrograms, sorting t-series in each condition by the time of the first post-stim peak >= threshold 
% Choose which cells & tseries to include & what method of normalization to
% use. Run Section A before running this section.
% 1. Combine together dF/F values from events from chosen regions in each tseries
%    and record the mean dF/F values for each tseries/recording as a row in
%    XsecBin_dffMat. Centering recordings around uncaging.
% 2. Fig. 2h: dF/F raster plots (normalized between 0-1 for each recording)
%    Sorting t-series (rows) in each condition by the time of the first post-stim
%    peak 3std above the baseline mean
%    Calls the function Peak3STD.m
%      For Fig. 2h, plotting events from the stimulated cell: 
%            WE = 'UncagingEvents'; SortOn = 'UncagingEvents';
%            ExcludeBaseline = 60; UseAreaThresh = 0; pulse_num = 10; Xsec = 1;
% 3. Fig. 3c: Sort t-series in each condition by the time of the first post-stim
%    peak, represent each t-series/recording (row) as z-score
%    Calls function Peak3STD.m
%      For Fig. 3c, plotting events from the neighboring cells: 
%            WE = 'NeighboringEvents'; SortOn = 'NeighboringEvents';
%            ExcludeBaseline = 60; UseAreaThresh = 0; pulse_num = 10; Xsec = 1;

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig3c';

ExcludeBaseline = 60; %number of seconds at the beginning of baseline to exclude, 0 to include all baseline period
EBL = sprintf('Exclude%dSBL', ExcludeBaseline);

WE = 'NeighboringEvents'; %'NeighboringEvents', 'UncagingEvents', 'AllEvents' which events should be included?
SortOn = 'NeighboringEvents'; %which events should the t-series be sorted by? 'NeighboringEvents', 'UncagingEvents'

UseAreaThresh = 0; %1 to plot Res or NonRes, 0 to plot all TS in a condition
ResNonRes = 1;% If UseAreaThresh == 1, 1: plot responders, 0: plot non-responders
PostStimTime = []; %time period post-stim for determining responders

pulse_num = 10;

XsecBin_dffMat = struct();
Xsec = 1; %size of desired bins (in seconds)
BinSize = sprintf('TimeBins%dS', Xsec);

if UseAreaThresh == 1
    if ResNonRes == 1
        WTS = 'Res';
    else
        WTS = 'NRes';
    end
else
    WTS = 'AllTS';
end

condition = fieldnames(mydata);
TSERIES = struct();
for cond = 1:length(condition)
    if exist('pulse_num', 'var')
        pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num));
    else
        pulse_cond = extractAfter(condition{cond}, '_');
    end
    if UseAreaThresh
        if isempty(PostStimTime)
            RESPONDER_TS = dFFResponders.poststim;
        else
            temptime = sprintf('poststim%ds', PostStimTime);
            RESPONDER_TS = dFFResponders.(temptime);
            clear temptime
        end
        
        if ResNonRes == 1
            TSERIES.(condition{cond}) = intersect(RESPONDER_TS.UncagingEvents.(condition{cond}), pulse_x_idx.(pulse_cond)); %return the t-series of responders within that pulse condition
        elseif ResNonRes == 0
            TSERIES.(condition{cond}) = setdiff(pulse_x_idx.(pulse_cond), RESPONDER_TS.UncagingEvents.(condition{cond})); %return the t-series of non-responders within that pulse condition
        end
    else
        TSERIES.(condition{cond}) = intersect(pulse_x_idx.(pulse_cond), 1:length(mydata));
    end
    clear pulse_cond RESPONDER_TS
end
clear cond temp_idx

% 1. Combine together dF/F values from events from chosen regions in each tseries
%    and record the mean dF/F values for each tseries/recording as a row in
%    XsecBin_dffMat. Centering recordings around uncaging.
condition = fieldnames(mydata);
EVENTS = {'UncagingEvents', 'NeighboringEvents'}; %'UncagingEvents': average together the dF/F traces from all events within the directly stimulated cell; 'NeighboringEvents': average together the dF/F traces from all events within cells not directly stimulated by uncaging;
for e = 1:length(EVENTS)
    for cond = 1:length(condition)
        %create a matrix of nans with a row for each tseries to store the trace
        %of choice in XSec bins
        XsecBin_dffMat.(EVENTS{e}).(condition{cond}) = NaN(max(TSERIES.(condition{cond})),1000);
        ts = 1;
        for tseries = TSERIES.(condition{cond})
            if ~isempty(mydata(tseries).(condition{cond}))
                SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
                ExcludeBaselineFrames = round(ExcludeBaseline / SecPerFrame);
                FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe - ExcludeBaselineFrames;
                tempFT = mydata(tseries).(condition{cond}).(EBL).(EVENTS{e}).dffMatFilter2(:, ExcludeBaselineFrames+1:end);
                FT = mean(tempFT, 1, 'omitnan');
                
                if ~isempty(FT)
                    NumFrames_Xsec = round(Xsec/SecPerFrame); %The number of frames in x seconds, listed above as 5sec (Xsec = 5)
                    FramesCroppedBeginning = mod(FrameOfIncrease-1, NumFrames_Xsec);
                    FramesCroppedEnd = mod((length(FT)- FrameOfIncrease + 1), NumFrames_Xsec);
                    
                    PerFieldActive_regPre = FT(FramesCroppedBeginning+1 : FrameOfIncrease-1);
                    PerFieldActivePerXSec_regPre = mean(reshape(PerFieldActive_regPre, NumFrames_Xsec, []),1, 'omitnan');
                    
                    PerFieldActive_regPost = FT(FrameOfIncrease : end-FramesCroppedEnd);
                    PerFieldActivePerXSec_regPost = mean(reshape(PerFieldActive_regPost, NumFrames_Xsec, []),1, 'omitnan');
                    
                    XsecBin_dffMat.(EVENTS{e}).(condition{cond})(tseries, (500-length(PerFieldActivePerXSec_regPre)+1):500) = PerFieldActivePerXSec_regPre;
                    XsecBin_dffMat.(EVENTS{e}).(condition{cond})(tseries, (500+1):(500+length(PerFieldActivePerXSec_regPost))) = PerFieldActivePerXSec_regPost;
                end
                ts = ts+1;
                clear SecPerFrame ExcludeBaselineFrames FrameOfIncrease tempFT FT...
                    NumFrames_XSec FramesCroppedBeginning FramesCroppedEnd PerFieldActive_regPre...
                    PerFieldActivePerXSec_regPre PerFieldActive_regPost PerFieldActivePerXSec_regPost
            end
        end
        clear ts tseries
    end
end
clear cond EVENTS e

% 2. Fig. 2h: dF/F raster plots (normalized between 0-1 for each recording)
%    Sorting t-series (rows) in each condition by the time of the first post-stim
%    peak 3std above the baseline mean. Greyed out rows at the bottom are
%    recordings with no detected events in the chosen cells
for cond = 1:length(condition)
    for NORM = 1
        mat2plot = XsecBin_dffMat.(WE).(condition{cond})(TSERIES.(condition{cond}),:);
        mat2sort = XsecBin_dffMat.(SortOn).(condition{cond})(TSERIES.(condition{cond}),:);
        ValBin = find(sum(~isnan(mat2plot),1)>0); %Find only the bins with values, not nan
        first_bin = ValBin(1);
        last_bin = ValBin(end);
        SecPerFrame = mydata(1).(condition{cond}).SecPerFrame;
        ExcludeBaselineFrames = round(ExcludeBaseline /SecPerFrame);
        mat2plot_cropped = mat2plot(:, first_bin:last_bin);
        mat2sort_cropped = mat2sort(:, first_bin:last_bin);
        if Xsec == 1
            XVALS = 1:(size(mat2plot_cropped,2)+ExcludeBaselineFrames);
            XVALS = (XVALS - mydata(1).(condition{cond}).uncagingframe) * SecPerFrame;
            XVALS = XVALS(ExcludeBaselineFrames+1:end);
        else
            
        end
        
        %finding the first peak 3STD above baseline mean and all frames
        %3STD above baseline mean
        First_Peak = NaN(1, size(mat2plot_cropped,1));
        FramesAboveThresh = NaN(size(mat2plot_cropped));
        for tseries = 1:size(mat2plot_cropped,1)
            First_Peak(1, tseries) = Peak3STD(mat2sort_cropped(tseries,:),mydata(1).(condition{cond}).uncagingframe-ExcludeBaselineFrames, [1 mydata(1).(condition{cond}).uncagingframe-ExcludeBaselineFrames]);
            
            BL = mat2plot_cropped(tseries, 1:mydata(1).(condition{cond}).uncagingframe-(ExcludeBaselineFrames+1));
            BL_Mean = nanmean(BL);
            BL_STD = nanstd(BL);
            PeakThreshold = BL_Mean + (3*BL_STD);
            FramesAboveThresh(tseries, :) = mat2plot_cropped(tseries,:) >= PeakThreshold;
            
            if isnan(First_Peak(1, tseries)) && sum(FramesAboveThresh(tseries,:)) > 0
                First_Peak(1, tseries) = size(mat2plot_cropped, 2); %If there are no peaks detected, but there is 1 post-stim frame above threshold, set the onset frame as the final frame
            elseif isnan(sum(mat2plot_cropped(tseries,:))) %If there are no values for that t-series, place the tseries at the bottom of the sorted tseries
                First_Peak(1, tseries) = size(mat2plot_cropped, 2) + 5;
            end
            clear BL BL_Mean BL_STD PeakThreshold
        end
        
        First_Peak(isnan(First_Peak)) = size(mat2plot_cropped, 2)+1; %any t-series without a peak should put the time of the peak as the final frame
        [~, FirstPeak_SortedIDX] = sort(First_Peak);
        
        mat2plot_cropped_sorted = mat2plot_cropped(FirstPeak_SortedIDX, :);
        FramesAboveThresh_sorted = FramesAboveThresh(FirstPeak_SortedIDX, :);
        clear First_Peak tseries FirstPeak_SortedIDX 
        
        if NORM == 0 %to plot un-normalized data
            MAT2PLOT = mat2plot_cropped_sorted;
            NORM_BY = '';
            clims = [0 2]; %[0 15];
        elseif NORM == 1
            % Normalizing each trace between 0 and 1
            Norm01 = NaN(size(mat2plot_cropped_sorted));
            for tseries = 1:size(mat2plot_cropped_sorted,1)
                temp_trace = mat2plot_cropped_sorted(tseries,:);
                MIN_X = min(temp_trace);
                MAX_X = max(temp_trace);
                Norm01(tseries,:) = (temp_trace - MIN_X) ./ (MAX_X - MIN_X);
                clear temp_trace MIN_X MAX_X
            end
            MAT2PLOT = Norm01;
            NORM_BY = 'Norm 0-1';
            clims = [0 1]; %[0 15];
        end
        
        if exist('pulse_num', 'var')
            TITLE = {sprintf('%s %s %s %s %d-100ms pulses', indicator, experiment_type, strrep(condition{cond}, '_', ' '), WE, pulse_num),...
                sprintf('Avg %s %s sorted on %s', 'Fluorescence of Identified Events', NORM_BY, SortOn),...
                sprintf('%s %s (%dS PostStim)', WE, WTS, PostStimTime), sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
        else
            TITLE = {sprintf('%s %s %s-100ms pulses %s', indicator, experiment_type, strrep(condition{cond}, '_', ' '), WE),...
                sprintf('Avg %s %s sorted on %s', 'Fluorescence of Identified Events', NORM_BY, SortOn),...
                sprintf('%s %s (%dS PostStim)', WE, WTS, PostStimTime), sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
        end
        
        % Plot spectrograms with each row as a t-series, sorted by the onset
        % time of post-stim peaks & binarized raster plotting frames in
        % each t-series 3STD above baseline mean
        figure() %For each condition plot the % of active px per frame with the mean across regions/t-series beneath
        ax1 = subplot(2,1,1);
        I = imagesc([XVALS(1) XVALS(end)], [1 size(MAT2PLOT, 1)], MAT2PLOT, clims);
        TRANSPARENCY = double(~isnan(MAT2PLOT));
        TRANSPARENCY(TRANSPARENCY==0) = 0.1;
        set(I, 'AlphaData', TRANSPARENCY)
        colorbar
        hold on
        vline(0, 'w--')
        ylabel('t-series')
        if Xsec == 1
            xlabel('time (s)')
        else
            xlabel(sprintf('time (%dS bins)', Xsec))
        end
        title(TITLE)
        
        ax2 = subplot(2,1,2);
        imagesc([XVALS(1) XVALS(end)], [1 size(FramesAboveThresh_sorted, 1)], FramesAboveThresh_sorted, [0 1])
        colormap(ax2, flipud(gray))
        colorbar
        hold on
        vline(0, 'r--')
        ylabel('t-series')
        if Xsec == 1
            xlabel('time (s)')
        else
            xlabel(sprintf('time (%dS bins)', Xsec))
        end
        title(TITLE)
        linkaxes([ax1, ax2], 'x')
        
        cd(save_dir)
        saveas(gcf, sprintf('%s_AvgFluorSpect_3std_%dSBins_ExRamp%dExBl%dS_%s%s_%dSpost_SORTED%s_NORM%d.tif', condition{cond}, Xsec, ExcludeRampingCells, ExcludeBaseline, WE, WTS, PostStimTime, SortOn, NORM))
        saveas(gcf, sprintf('%s_AvgFluorSpect_3std_%dSBins_ExRamp%dExBl%dS_%s%s_%dSpost_SORTED%s_NORM%d.svg', condition{cond}, Xsec, ExcludeRampingCells, ExcludeBaseline, WE, WTS, PostStimTime, SortOn, NORM))
        close
        cd(start_dir)

        clear I TRANSPARENCY mat2plot MAT2PLOT first_bin last_bin mat2plot_cropped ValBin...
            clims ax1 ax2 SEM mean_plot exclude_bins SecPerFrame XLIMS...
            MAT Norm01 tseries Norm01_LPFILT MAT YSPACER TITLE mat2plot_cropped_sorted...
            NORM_BY XVALS lowpassfilt UNC_BIN YLIMS YLABEL AVG_Pre AVG_Post FramesAboveThresh_sorted...
            mat2sort mat2sort_cropped
    end
    clear NORM
end

% 3. Fig. 3c: dF/F raster plots (represented as the z-score for each 
%    recording/tseries [standard deviations from baseline mean]) &
%    binarized raster plots marking frames in each recording 3STD >=
%    baseline mean (z-score >= 3).
%    Sorting t-series (rows) in each condition by the time of the first post-stim
%    peak 3std above the baseline mean. Greyed out rows at the bottom are
%    recordings with no detected events in the chosen cells 
condition = fieldnames(mydata);
for cond = 1:length(condition)
        mat2plot = XsecBin_dffMat.(WE).(condition{cond})(TSERIES.(condition{cond}),:);
        mat2sort = XsecBin_dffMat.(SortOn).(condition{cond})(TSERIES.(condition{cond}),:);
        ValBin = find(sum(~isnan(mat2plot),1)>0); %Find only the bins with values, not nan
        first_bin = ValBin(1);
        last_bin = ValBin(end);
        SecPerFrame = mydata(1).(condition{cond}).SecPerFrame;
        ExcludeBaselineFrames = round(ExcludeBaseline /SecPerFrame);
        mat2plot_cropped = mat2plot(:, first_bin:last_bin);
        mat2sort_cropped = mat2sort(:, first_bin:last_bin);
        if Xsec == 1
            XVALS = 1:(size(mat2plot_cropped,2)+ExcludeBaselineFrames);
            XVALS = (XVALS - mydata(1).(condition{cond}).uncagingframe) * SecPerFrame;
            XVALS = XVALS(ExcludeBaselineFrames+1:end);
        else
            
        end
        
        % 1. finding the first peak 3STD above baseline mean and all frames
        %    3STD above baseline mean
        % 2. representing each frame as standard deviation from baseline
        %    mean
        First_Peak = NaN(1, size(mat2plot_cropped,1));
        FramesAboveThresh = NaN(size(mat2plot_cropped));
        SDFromMean = NaN(size(mat2plot_cropped));
        for tseries = 1:size(mat2plot_cropped,1)
            %1.
            First_Peak(1, tseries) = Peak3STD(mat2sort_cropped(tseries,:),mydata(1).(condition{cond}).uncagingframe-ExcludeBaselineFrames, [1 mydata(1).(condition{cond}).uncagingframe-ExcludeBaselineFrames]);
            
            BL = mat2plot_cropped(tseries, 1:mydata(1).(condition{cond}).uncagingframe-(ExcludeBaselineFrames+1));
            BL_Mean = nanmean(BL);
            BL_STD = nanstd(BL);
            PeakThreshold = BL_Mean + (3*BL_STD);
            FramesAboveThresh(tseries, :) = mat2plot_cropped(tseries,:) >= PeakThreshold;
            
            if isnan(First_Peak(1, tseries)) && sum(FramesAboveThresh(tseries,:)) > 0
                First_Peak(1, tseries) = size(mat2plot_cropped, 2); %If there are no peaks detected, but there is 1 post-stim frame above threshold, set the onset frame as the final frame
            elseif isnan(sum(mat2plot_cropped(tseries,:))) %If there are no values for that t-series, place the tseries at the bottom of the sorted tseries
                First_Peak(1, tseries) = size(mat2plot_cropped, 2) + 5; 
            end
            
            %2.
            SDFromMean(tseries,:) = (mat2plot_cropped(tseries,:) - BL_Mean) ./ BL_STD; 
            clear BL BL_Mean BL_STD PeakThreshold
        end
        
        First_Peak(isnan(First_Peak)) = size(mat2plot_cropped, 2)+1; %any t-series without a peak should put the time of the peak as the final frame + 1
        [~, FirstPeak_SortedIDX] = sort(First_Peak);
        
        MAT2PLOT = SDFromMean(FirstPeak_SortedIDX, :);
        FramesAboveThresh_sorted = FramesAboveThresh(FirstPeak_SortedIDX, :);
        clear First_Peak tseries FirstPeak_SortedIDX SDFromMean
        
        NORM_BY = 'SD from BL mean';
        clims = [-4 4];
        
        if exist('pulse_num', 'var')
            TITLE = {sprintf('%s %s %s %s %d-100ms pulses', indicator, experiment_type, strrep(condition{cond}, '_', ' '), WE, pulse_num),...
                sprintf('Avg %s %s sorted on %s', 'Fluorescence of Identified Events', NORM_BY, SortOn),...
                sprintf('%s %s (%dS PostStim)', WE, WTS, PostStimTime), sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
        else
            TITLE = {sprintf('%s %s %s-100ms pulses %s', indicator, experiment_type, strrep(condition{cond}, '_', ' '), WE),...
                sprintf('Avg %s %s sorted on %s', 'Fluorescence of Identified Events', NORM_BY, SortOn),...
                sprintf('%s %s (%dS PostStim)', WE, WTS, PostStimTime), sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
        end
        
        COLOR = flipud(cbrewer('div', 'RdBu', (numel(clims(1):clims(2))-1)*1+1));
        
        % Plot spectrograms with each row as a t-series, sorted by the onset
        % time of post-stim peaks & binarized raster plotting frames in
        % each t-series 3STD above baseline mean
        figure() %For each condition plot the % of active px per frame with the mean across regions/t-series beneath
        ax1 = subplot(2,1,1);
        I = imagesc([XVALS(1) XVALS(end)], [1 size(MAT2PLOT, 1)], MAT2PLOT, clims);
        TRANSPARENCY = double(~isnan(MAT2PLOT));
        TRANSPARENCY(TRANSPARENCY==0) = 0.5; %0.1;
        set(I, 'AlphaData', TRANSPARENCY)
        if exist('COLOR', 'var')
            colormap(COLOR)
        end
        c = colorbar;
        c.Label.String = 'Z-score';
        hold on
        if exist('COLOR', 'var')
            vline(0, 'k--')
        else
            vline(0, 'w--')
        end
        ylabel('t-series')
        if Xsec == 1
            xlabel('time (s)')
        else
            xlabel(sprintf('time (%dS bins)', Xsec))
        end
        title(TITLE)
        
        ax2 = subplot(2,1,2);
        imagesc([XVALS(1) XVALS(end)], [1 size(FramesAboveThresh_sorted, 1)], FramesAboveThresh_sorted, [0 1])
        colormap(ax2, flipud(gray))
        colorbar
        hold on
        vline(0, 'r--')
        ylabel('t-series')
        if Xsec == 1
            xlabel('time (s)')
        else
            xlabel(sprintf('time (%dS bins)', Xsec))
        end
        title(TITLE)
        linkaxes([ax1, ax2], 'x')
        
        cd(save_dir)
        saveas(gcf, sprintf('%s_AvgFluorSpect_3std_%dSBins_ExRamp%dExBl%dS_%s%s_%dSpost_SORTED%s_SDofBLMean.tif', condition{cond}, Xsec, ExcludeRampingCells, ExcludeBaseline, WE, WTS, PostStimTime, SortOn))
        saveas(gcf, sprintf('%s_AvgFluorSpect_3std_%dSBins_ExRamp%dExBl%dS_%s%s_%dSpost_SORTED%s_SDofBLMean.svg', condition{cond}, Xsec, ExcludeRampingCells, ExcludeBaseline, WE, WTS, PostStimTime, SortOn))
        close
        cd(start_dir)
        
        clear I TRANSPARENCY c mat2plot MAT2PLOT first_bin last_bin mat2plot_cropped ValBin...
            clims ax1 ax2 SEM mean_plot exclude_bins SecPerFrame XLIMS...
            MAT Norm01 tseries Norm01_LPFILT MAT YSPACER TITLE mat2plot_cropped_sorted...
            NORM_BY XVALS lowpassfilt UNC_BIN YLIMS YLABEL AVG_Pre AVG_Post FramesAboveThresh_sorted...
            mat2sort mat2sort_cropped COLOR
end

clear e cond EVENTS AnyDifferences CheckThresholds cond SortOn
%% Section C. Fig. 3c stacked bar graphs showing proportion of recordings 
%  exhibiting an initial fluorescence increase following uncaging
%  (responder). Fisher's exact test compares proportion of responders
%  across WT and Cx43 conditions. 
%  Run section A & B before running this section. 
%  Modified from FluorescenceTraceFromAQuAEvtsNoRegion_MC20211117.m
%  Calls the function cbrewer
%
%  For Fig. 3c: pulse_num = 10; Xsec = 1; EVENTS = {'NeighboringEvents'};

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig3c\StackedBarCharts';

pulse_num = 10;

Xsec = 1; 
FramesAboveThresh_XSecBins = struct();
TotalFramesAboveThresh = struct();

condition = fieldnames(mydata);
TSERIES = struct();
for cond = 1:length(condition)
    if exist('pulse_num', 'var')
        pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num));
    else
        pulse_cond = extractAfter(condition{cond}, '_');
    end
    if UseAreaThresh
        if isempty(PostStimTime)
            RESPONDER_TS = dFFResponders.poststim;
        else
            temptime = sprintf('poststim%ds', PostStimTime);
            RESPONDER_TS = dFFResponders.(temptime);
            clear temptime
        end
        
        if ResNonRes == 1
            TSERIES.(condition{cond}) = intersect(RESPONDER_TS.UncagingEvents.(condition{cond}), pulse_x_idx.(pulse_cond)); %return the t-series of responders within that pulse condition
        elseif ResNonRes == 0
            TSERIES.(condition{cond}) = setdiff(pulse_x_idx.(pulse_cond), RESPONDER_TS.UncagingEvents.(condition{cond})); %return the t-series of non-responders within that pulse condition
        end
    else
        TSERIES.(condition{cond}) = intersect(pulse_x_idx.(pulse_cond), 1:length(mydata));
    end
    clear pulse_cond RESPONDER_TS
end
clear cond temp_idx

% For each recording, identify the frames with dF/F values above threshold
% and the recordings with any post-stim frames above threshold (responder)
% For either dF/F values from uncaging cell events or neighboring cell
% events
EVENTS = {'NeighboringEvents'}; %{'UncagingEvents', 'NeighboringEvents'};
for e = 1:length(EVENTS)
    for cond = 1:length(condition)
        Approx_Bins = ceil((mydata(1).(condition{cond}).totalframes * mydata(1).(condition{cond}).SecPerFrame) / Xsec);
        tempFramesAboveThresh = NaN(max(TSERIES.(condition{cond})), Approx_Bins);
        tempTotalFramesAboveThresh = NaN(max(TSERIES.(condition{cond})), 1);
        for tseries = TSERIES.(condition{cond})
            if ~isempty(mydata(tseries).(condition{cond})) && ~isempty(mydata(tseries).(condition{cond}).(EBL).(EVENTS{e}).dffMatFilter2)
                SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
                if Xsec == 1
                    ExcludeBaselineBins = round(ExcludeBaseline / SecPerFrame);
                else
                    ExcludeBaselineBins = round(ExcludeBaseline / Xsec);
                end
                FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe;
                NumFrames_Xsec = round(Xsec/SecPerFrame); %The number of frames in x seconds, listed above as 5sec (Xsec = 5)
                
                tempFT = mydata(tseries).(condition{cond}).(EBL).(EVENTS{e}).dffMatFilter2_Thresh3std;
                
                FramesCroppedBeginning = mod(FrameOfIncrease-1, NumFrames_Xsec);
                FramesCroppedEnd = mod((length(tempFT)- FrameOfIncrease + 1), NumFrames_Xsec);
                
                tempPre = tempFT(FramesCroppedBeginning+1 : FrameOfIncrease-1);
                XSec_Pre = sum(reshape(tempPre, NumFrames_Xsec, []),1);
                
                tempPost = tempFT(FrameOfIncrease : end-FramesCroppedEnd);
                XSec_Post = sum(reshape(tempPost, NumFrames_Xsec, []),1);
                
                NumBins.pre = length(XSec_Pre);
                NumBins.post = length(XSec_Post);
                
                tempFramesAboveThresh(tseries, 1:NumBins.pre+NumBins.post) = [XSec_Pre, XSec_Post];
                tempTotalFramesAboveThresh(tseries, 1) = sum(tempFT(FrameOfIncrease:end));
                clear ExcludeBaselineBins FrameOfIncrease NumFrames_Xsec tempFT...
                    FramesCroppedBeginning FramesCroppedEnd tempPre XSec_Pre...
                    tempPost XSec_Post SecPerFrame
            end
        end
        FramesAboveThresh_XSecBins.(EVENTS{e}).(condition{cond}) = tempFramesAboveThresh(:,1:NumBins.pre+NumBins.post);
        TotalFramesAboveThresh.(EVENTS{e}).(condition{cond}) = tempTotalFramesAboveThresh;
        clear tseries tempFramesAboveThresh Approx_Bins tempTotalFramesAboveThresh
    end
end
clear cond e EVENTS

% Plot stacked bar graphs for each condition showing the proportion of 
% recordings exhibiting an initial fluorescence increase 
% (out of all recordings in that condition with detected events in the chosen cells)
% Choose which conditions to compare in Fisher's exact test and which
% conditions to exclude (recorded in FishersExact_pval)
condition2plot = {'RuBiGABA', 'RuBiGlu', 'WT'}; % Which conditions should be plotted and compared? Each cell/string in condition2plot will be run separately and will include all conditions containing the string. 'RuBiGABA_NoRuBi', 'RuBiGlu_NoRuBi', 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' 
EXCLUDE = {'CBX', 'CBX', ''}; %Within that dataset, should any condition be left out?
COLORS = {'RdPu', 'Blues', 'MgCy'}; %Color corresponding to cbrewer; 'MgCy' is not a cbrewer color scheme and can be used for datasets with 2 conditions, coloring one magenta and one cyan
COLOR_Scheme = {'seq', 'seq', ''}; %Color schemes corresponding to cbrewer
WE = 'NeighboringEvents'; %'UncagingEvents', 'NeighboringEvents'
FishersExact_pval = struct(); %keep track of whether the number of responders and nonresponders from different compared conditions are similar or different
clear condition
for cp = 1:length(condition2plot)
    % Determine which conditions to include
    if strcmp(condition2plot{cp}, 'AllConditions')
        condition = fieldnames(mydata);
    else
        temp_condition = fieldnames(mydata);
        if ~contains(condition2plot{cp}, '_')
            idx = cell2mat(cellfun(@(x) contains(x, condition2plot{cp}), temp_condition, 'UniformOutput', 0));
            idx = find(idx);
        else %to pull out conditions with different strings (for example: All RuBiGABA conditions and NoRuBi_LaserUncagingControl would be used with 'RuBiGABA_NoRuBi'
            idx_1 = cell2mat(cellfun(@(x) contains(x, extractBefore(condition2plot{cp}, '_')), temp_condition, 'UniformOutput', 0));
            idx_2 = cell2mat(cellfun(@(x) contains(x, extractAfter(condition2plot{cp}, '_')), temp_condition, 'UniformOutput', 0));
            idx = find(idx_1 | idx_2);
        end
        if ~isempty(EXCLUDE{cp})
            exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
            idx = setdiff(idx, exclude_idx);
        end
        condition = temp_condition(idx, 1);
        clear temp_condition idx exclude_idx idx_1 idx_2
    end
    
    % Determine which color to use for each condition
    if ~isempty(COLOR_Scheme{cp})
        cmap = cbrewer(COLOR_Scheme{cp}, COLORS{cp}, 8,'PCHIP');
        if contains(condition2plot{cp}, '_')
            cmap(1,:) = [0.5 0.5 0.5];
        end
        if contains(condition2plot{cp}, '_') || strcmp(COLOR_Scheme{cp}, 'div')
            cmap_num = flip(floor(linspace(1, 8, length(condition))));
        else
            cmap_num = flip(floor(linspace(3, 8, length(condition))));
        end
    else
        cmap = [1 0 1; 0 1 1]; %magneta and cyan
        cmap_num = [1, 2];
    end

    % Calculate the total number of post-stim frames above threshold & the
    % number of recordings in each condition with post-stim frames above
    % threshold
    struct2plot = struct();
    struct2plot_Responders = struct();
    ProportionResponders = struct();
    NUM_Responders = struct();
    Ns = struct();
    for cond = 1:length(condition)
        temp = TotalFramesAboveThresh.(WE).(condition{cond});
        struct2plot.(condition{cond}) = temp;
        struct2plot_Responders.(condition{cond}) = temp(temp>0);
        NumResp = sum(temp>0);
        TotalTS = sum(~isnan(temp));
        ProportionResponders.(condition{cond}) = NumResp/TotalTS;
        NUM_Responders.(condition{cond}) = [NumResp , TotalTS - NumResp ]; %column1: # of responders, column2: # of non-responders
        Ns.(condition{cond}) = [NumResp, TotalTS]; %column1: # of responders, column2: total number of TS with detected events following exclusion of ramping cells
        clear temp NumResp TotalTS
    end

    % Plot stacked bar graph for each condition for % of t-series with post-stim
    % frames over threshold (Fig. 3c, to the right of binarized raster
    % plots. In Fig. 3c, the height of each stacked bar graph is fit to the
    % height of the corresponding binarized raster, excluding an recordings
    % with out detected events [bottom-most rows greyed out])
    for cond = 1:length(condition)
        if exist('pulse_num', 'var')
            TITLE = {sprintf('%s %s %s %d-100ms pulses', indicator, experiment_type, condition{cond}, pulse_num),...
                'proportion of tseries with post-stim frames above threshold(dF/F)',...
                sprintf('%s %s: n=%d responders, %d total TS', WE, WTS, Ns.(condition{cond})(1), Ns.(condition{cond})(2)),...
                sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
        else
            TITLE = {sprintf('%s %s %s-100ms pulses', indicator, experiment_type, condition{cond}),...
                'proportion of tseries with post-stim frames above threshold(dF/F)',...
               sprintf('%s %s: n=%d responders, %d total TS', WE, WTS, Ns.(condition{cond})(1), Ns.(condition{cond})(2)),...
               sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
        end
        RESPONDERS = ProportionResponders.(condition{cond});
        NONRESPONDERS = 1-RESPONDERS;
        
        Y = [NONRESPONDERS, RESPONDERS; 0,0];
        figure()
        b = bar(Y, "stacked"); %cmap, cmap_num
        b(2).FaceColor = cmap(cmap_num(cond), :); %responders will be the color specified by cmap
        b(1).FaceColor = [0.5 0.5 0.5]; %non-responders will be grey
        legend({'non-responders', 'responders'}, 'Location', 'best')
        title(TITLE)
        xlim([.3 1.7])

        cd(save_dir)
        saveas(gcf, sprintf('%s_ProportionTSWithPostFramesAbvThresh_StackedBar_dFF_ExRamp%dExBl%dS_%s%s.tif', condition{cond}, ExcludeRampingCells, ExcludeBaseline, WE, WTS))
        saveas(gcf, sprintf('%s_ProportionTSWithPostFramesAbvThresh_StackedBar_dFF_ExRamp%dExBl%dS_%s%s.svg', condition{cond}, ExcludeRampingCells, ExcludeBaseline, WE, WTS))
        close
        cd(start_dir)
    
        clear NONRESPONDERS RESPONDERS b TITLE Y
    end
    
    %Run a Fisher's exact test to see if there's a difference between the # of responders from different conditions
    TwoByTwo = NaN(2,2);
    if length(condition) == 2 
        TwoByTwo = [NUM_Responders.(condition{1}); NUM_Responders.(condition{2})];
        Args_FishersExact = table(TwoByTwo(:,1),TwoByTwo(:,2), 'VariableNames', {'Responders', 'Nonresponders'}, 'RowNames', {condition{1}, condition{2}});
        [~, FishersExact_pval.(condition2plot{cp})] = fishertest(Args_FishersExact);
        clear TwoByTwo Args_FishersExact
    end
        
    clear ExcludeBaselineBins struct2plot PlotTS cond TITLE YLABEL YLIMS XTICKS...
        XLIMS XTICKLABELS DataLabels cmap cmap_num ProportionTS_aboveThresh STAT cond...
        Mod_TITLE Additional_TITLE struct2plot_Responders ProportionResponders
end

clear save_dir start_dir Xsec FramesAboveThresh_XSecBins condition2plot cp...
    EXCLUDE COLORS COLOR_Scheme condition NumBins TotalFramesAboveThresh PAIRED FishersExact_pval...
    BinSize EBL ERRORBARS ExcludeBaseline ExcludeBaselineFrames FramesAboveThresh Ns NUM_Responders...
    PostStimTime ResNonRes TestDifference UseAreaThresh WE WTS XsecBin_dffMat TSERIES
%% Section D. Fig. 3d&e: Latency to post-stim frames reaching a threshold for dF/F values
% Modified from FluorescenceTraceFromAQuAEvtsNoRegion_MC20211117.m
% Run section A before running this section
% 1. Record time to first post-stim frame above threshold for each t-series
%    of each condition for the chosen groups of cells
%    For Fig. 3: EVENTS = {'NeighboringEvents'};
% 2. Distribution of time to first frame above threshold
%    a. Fig. 3d overlay of response onset (cumulative proportion of
%       t-series with post-stim frames above threshold, out of all t-series)
%    b. Fig. 3e violin plots for onset latency
%       Xsec = 1; condition2plot = {'RuBiGABA', 'RuBiGlu'}; PAIRED = [0,0];
%       EXCLUDE = {'', ''};
% calls function violinplot_STAT.m and cbrewer.m

ExcludeBaseline = 60; %number of seconds at the beginning of baseline to exclude, 0 to include all baseline period
EBL = sprintf('Exclude%dSBL', ExcludeBaseline);
UseAreaThresh = 0; %1 to plot Res or NonRes, 0 to plot all TS in a condition
ERRORBARS = 'sem'; %'sem' or 'std'
pulse_num = 10;

condition = fieldnames(mydata);
TSERIES = struct();
for cond = 1:length(condition)
    if exist('pulse_num', 'var')
        pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num));
    else
        pulse_cond = extractAfter(condition{cond}, '_');
    end
    if UseAreaThresh
        if isempty(PostStimTime)
            RESPONDER_TS = dFFResponders.poststim;
        else
            temptime = sprintf('poststim%ds', PostStimTime);
            RESPONDER_TS = dFFResponders.(temptime);
            clear temptime
        end
        
        if ResNonRes == 1
            TSERIES.(condition{cond}) = intersect(RESPONDER_TS.UncagingEvents.(condition{cond}), pulse_x_idx.(pulse_cond)); %return the t-series of responders within that pulse condition
        elseif ResNonRes == 0
            TSERIES.(condition{cond}) = setdiff(pulse_x_idx.(pulse_cond), RESPONDER_TS.UncagingEvents.(condition{cond})); %return the t-series of non-responders within that pulse condition
        end
    else
        TSERIES.(condition{cond}) = intersect(pulse_x_idx.(pulse_cond), 1:length(mydata));
    end
    clear pulse_cond RESPONDER_TS
end
clear cond temp_idx

% 1. Record time to first post-stim frame above threshold for each t-series
%   of each condition for the chosen groups of cells
dFF_mat_props = struct();
EVENTS = {'NeighboringEvents'}; %{'UncagingEvents', 'NeighboringEvents'};
condition = fieldnames(mydata);
for e = 1:length(EVENTS)
    for cond = 1:length(condition)
        tempTimeFirstFrameAbvThresh = NaN(length(mydata), 1);
        for tseries = TSERIES.(condition{cond})
            if ~isempty(mydata(tseries).(condition{cond}))
                SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
                FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe;
                ExcludeBaselineFrames = round(ExcludeBaseline / SecPerFrame);
                
                temp_trace = mean(mydata(tseries).(condition{cond}).(EBL).(EVENTS{e}).dffMatFilter2, 1, 'omitnan');%(:, ExcludeBaselineFrames+1:end);
                if ~isempty(temp_trace)
                    tempTotalFramesAboveThresh(tseries, 1) = sum(mydata(tseries).(condition{cond}).(EBL).(EVENTS{e}).dffMatFilter2_Thresh3std(FrameOfIncrease:end));
                    %Find time to the first post-stim frame >= threshold
                    %(uncaging frame would be considered 0s)
                    if tempTotalFramesAboveThresh(tseries, 1) > 0
                        PostStimFramesAbvThresh = mydata(tseries).(condition{cond}).(EBL).(EVENTS{e}).dffMatFilter2_Thresh3std(FrameOfIncrease:end);
                        FirstPostStimFrameAbvThresh = find(PostStimFramesAbvThresh, 1);
                        tempTimeFirstFrameAbvThresh(tseries, 1) = (FirstPostStimFrameAbvThresh - 1) * SecPerFrame;
                        clear PostStimFramesAbvThresh FirstPostStimFrameAbvThresh
                    end
                end
                clear SecPerFrame FrameOfIncrease ExcludeBaselineFrames temp_trace...
                    First_Peak p_num p t TB BLBinsExclude Nonzero_cells ramping_cell_idx
            end
        end
        dFF_mat_props.(EVENTS{e}).(condition{cond}).TimeFirstAbvThreshFrame = tempTimeFirstFrameAbvThresh;
        clear tempTotalFramesAboveThresh tempAvgdFF_Pre tempTimeFirstPeak tempdFF_FirstPeak tseries...
            tempAQuAEvt_params p_num p t tempTimeFirstFrameAbvThresh
    end
end
clear cond e EVENTS PreOrPost params

% 2. Distribution of time to first frame above threshold
%    Choose which comparisons to make (using 'condition2plot')
%   a. Fig. 3d overlay: CDF plots showing cumulative proportion of t-series
%      with 1st frame above threshold: out of all t-series 
%   b. Fig. 3e: violin plots for onset latency
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig3d_e';
Xsec = 1;
condition2plot = {'RuBiGABA', 'RuBiGlu'};% Which conditions should be plotted and compared? Each cell/string in condition2plot will be run separately and will include all conditions containing the string.  'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' Which conditions should be plotted
PAIRED = [0,0]; %1 if the data within that set is paired or 0 if the compare conditions are unpaired
EXCLUDE = {'', ''}; %Within that dataset, should any condition be left out?
COLORS = {'RdPu', 'Blues'}; %Color corresponding to cbrewer; 'MgCy' is not a cbrewer color scheme and can be used for datasets with 2 conditions, coloring one magenta and one cyan
COLOR_Scheme = {'seq', 'seq'}; %Color schemes corresponding to cbrewer
EVENTS = {'NeighboringEvents'}; %'UncagingEvents',

clear condition
for cp = 1:length(condition2plot)
    if strcmp(condition2plot{cp}, 'AllConditions')
        condition = fieldnames(mydata);
    else
        temp_condition = fieldnames(mydata);
        idx = find(cell2mat(cellfun(@(x) contains(x, condition2plot{cp}), temp_condition, 'UniformOutput', 0)));
        if ~isempty(EXCLUDE{cp})
            exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
            idx = setdiff(idx, exclude_idx);
        end
        condition = temp_condition(idx, 1);
        clear temp_condition idx exclude_idx
    end
    
    if ~isempty(COLOR_Scheme{cp})
        cmap = cbrewer(COLOR_Scheme{cp}, COLORS{cp}, 8,'PCHIP');
        if contains(condition2plot{cp}, '_')
            cmap(1,:) = [0.5 0.5 0.5];
        end
        if contains(condition2plot{cp}, '_') || strcmp(COLOR_Scheme{cp}, 'div')
            cmap_num = flip(floor(linspace(1, 8, length(condition))));
        else
            cmap_num = flip(floor(linspace(3, 8, length(condition))));
        end
    else
        cmap = [1 0 1; 0 1 1]; %magneta and cyan
        cmap_num = [1, 2];
    end
    
    for e=1:length(EVENTS)
        struct2plot = struct();
        for cond = 1:length(condition)
            temp = dFF_mat_props.(EVENTS{e}).(condition{cond}).TimeFirstAbvThreshFrame;
            struct2plot.(condition{cond}) = temp;
            clear temp
        end
        
        if exist('pulse_num', 'var')
            TITLE = {sprintf('%s %s %s %d-100ms pulses', indicator, experiment_type, condition2plot{cp}, pulse_num),...
                'latency to first frame above thresh',...
                sprintf('%s Exclude Ramping Cells: %d, %s', EVENTS{e}, ExcludeRampingCells, EBL)};
        else
            TITLE = {sprintf('%s %s %s', indicator, experiment_type, condition2plot{cp}),...
                'latency to first frame above thresh',...
                sprintf('%s Exclude Ramping Cells: %d, %s', EVENTS{e}, ExcludeRampingCells, EBL)};
        end

        YLABEL = 'latency to first frame above thresh';
        XLABEL = 'probability';
        DATA_LABELS = cellfun(@(x) strrep(x, '_', ' '), condition, 'UniformOutput', 0);
         
        EDGES = 0:Xsec:150;
        hist2plot = cell(1,length(condition));
        cdf2plot = cell(1,length(condition));
        cdf2plot_resp = cell(1,length(condition));
        PlotTS = cell(1,length(condition));
        for cond = 1:length(condition)
            TotalTS = numel(TSERIES.(condition{cond}));
            TotalResp = sum(~isnan(struct2plot.(condition{cond})));
            temp_count = histcounts(struct2plot.(condition{cond}), EDGES);
            temp_propTotal = temp_count ./ TotalTS;
            temp_propResp = temp_count ./ TotalResp;
            hist2plot{1,cond} = temp_propTotal;
            cdf2plot{1,cond} = cumsum(temp_propTotal, 2);
            cdf2plot_resp{1,cond} = cumsum(temp_propResp, 2);
            PlotTS{1,cond} = {1};
            clear TotalTS temp_count temp_propTotal
        end

        XTICKS = 1:numel(EDGES)-1;
        XTICKS = XTICKS * Xsec - (Xsec/2);
        XLIMS = [EDGES(1) EDGES(end)+1];
        XTICKLABELS = {};
        DataLabels = cellfun(@(x) strrep(x, '_', ' '), condition, 'UniformOutput', 0);
        
        %plot the cumulative proportion of tseries with onset times in the given time
        %bin out of all t-series
        if length(cdf2plot) == 2
            [h, p]= kstest2(cdf2plot{1,1}, cdf2plot{1,2});
        else
            h = NaN;
            p = NaN;
        end
        
        TITLE_pvals = [TITLE, sprintf('ks test p=%.15f', p)];
        PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt( cdf2plot, PlotTS, 1, TITLE_pvals,...
            'cumulative proportion ts with response onset', [0 1], sprintf('%dS time bins', Xsec), XLIMS, XTICKS,XTICKLABELS,...
            DataLabels, {}, ERRORBARS, cmap, cmap_num )
%         if Xsec == 1
            xticks([0:10:150])
%         end
        
        cd(save_dir)
        saveas(gcf, sprintf('%sEx%s_LatencyToFirstFrameAbvThreshCumSum_%dSBins_dFF_ExRamp%dExBl%dS_%s.tif', condition2plot{cp}, EXCLUDE{cp}, Xsec, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}))
        saveas(gcf, sprintf('%sEx%s_LatencyToFirstFrameAbvThreshCumSum_%dSBins_dFF_ExRamp%dExBl%dS_%s.svg', condition2plot{cp}, EXCLUDE{cp}, Xsec, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}))
        close
        cd(start_dir)
        clear h p TITLE_pvals
        

        clear EDGES hist2plot cdf2plot cdf2plot_resp PlotTS cond XTICKS XLIMS XTICKLABELS DataLabels h p TITLE_pvals
        
        if PAIRED(cp) == 0
            if length(condition) == 2
                STAT = 'RankSum';
            else
                STAT = 'ANOVA';
            end
            
            violinplot_STATS(struct2plot, [], TITLE, YLABEL, cmap, cmap_num, STAT)
            cd(save_dir)
            saveas(gcf, sprintf('%sEx%s_LatencyToFirstFrameAbvThresh_Violin_dFF_ExRamp%dExBl%dS_%s.tif', condition2plot{cp}, EXCLUDE{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}))
            saveas(gcf, sprintf('%sEx%s_LatencyToFirstFrameAbvThresh_Violin_dFF_ExRamp%dExBl%dS_%s.svg', condition2plot{cp}, EXCLUDE{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}))
            close
            cd(start_dir)
            
        elseif PAIRED(cp) == 1
            Param_PreVPost = [struct2plot.(condition{1}),struct2plot.(condition{2})];
            NUM_TS = TSERIES.(condition{1});
            
            PlotTSMean_PreVPost(Param_PreVPost, NUM_TS, TITLE, YLABEL, [], ERRORBARS)
            xticklabels({strrep(condition{1}, '_', ' '), strrep(condition{2}, '_', ' ')})
            legend('off')
            cd(save_dir)
            saveas(gcf, sprintf('%sEx%s_LatencyToFirstFrameAbvThresh_dFF_ExRamp%dExBl%dS_%s.tif', condition2plot{cp}, EXCLUDE{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}))
            saveas(gcf, sprintf('%sEx%s_LatencyToFirstFrameAbvThresh_dFF_ExRamp%dExBl%dS_%s.svg', condition2plot{cp}, EXCLUDE{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}))
            close
            cd(start_dir)
            clear Param_PreVPost NUM_TS
        end
        clear struct2plot PlotTS cond TITLE YLABEL STAT cond
    end
    clear struct2plot PlotTS cond TITLE YLABEL...
        cmap cmap_num STAT cond e condition
end

clear cond EVENTS MEASUREMENT start_dir save_dir dFF_mat_props params PreOrPost...
    condition2plot PAIRED EXCLUDE COLORS COLOR_Scheme cp TSERIES ERRORBARS...
    DATA_LABELS temp_propResp TotalResp UseAreaThresh XLABEL Xsec pulse_num 
%% Section E. Fig. 3f: Response persistence. Proportion of post-stim frames
%  above threshold following first post-stim frame above threshold. 
%  Modified from FluorescenceTraceFromAQuAEvtsNoRegion_MC20211117.m
%  Run section A before running this section
% 1. Record proportion of post-stim frames above threshold for each t-series
%    of each condition for the chosen groups of cells
% 2. Fig. 3f: Violin plots of proportion of post-stim frames above threshold
%    Choose which comparisons to make (using 'condition2plot')
%    Calls function violinplot_STATS.m and cbrewer.m
%
%    UseAreaThresh = 0; EVENTS = {'NeighboringEvents'}; 
%    condition2plot = {'RuBiGABA', 'RuBiGlu'}; PAIRED = [0,0]; EXCLUDE = {'', ''}; 

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig3f';

UseAreaThresh = 0; %1 to plot Res or NonRes, 0 to plot all TS in a condition
ResNonRes = 1; %If UseAreaThresh == 1, 1: plot responders, 0: plot non-responders
PostStimTime = []; %time period (s) post-stim for determining responders if UseAreaThresh == 1; [] for all time post-stim 
pulse_num = 10;

if UseAreaThresh == 1
    if ResNonRes == 1
        WTS = 'Res';
    else
        WTS = 'NRes';
    end
else
    WTS = 'AllTS';
end

% 1. Record proportion of post-stim frames above threshold for each t-series
%   of each condition for the chosen groups of cells
condition = fieldnames(mydata);
PropFramesAboveThresh = struct();
EVENTS = {'NeighboringEvents'}; %{'UncagingEvents', 'NeighboringEvents'}; 
for e = 1:length(EVENTS)
    TSERIES = struct();
    for cond = 1:length(condition)
        pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num));
        %         pulse_cond = extractAfter(condition{cond}, '_');
        if UseAreaThresh
            if isempty(PostStimTime)
                RESPONDER_TS = dFFResponders.poststim;
            else
                temptime = sprintf('poststim%ds', PostStimTime);
                RESPONDER_TS = dFFResponders.(temptime);
                clear temptime
            end
            
            if ResNonRes == 1
                TSERIES.(condition{cond}) = intersect(RESPONDER_TS.(EVENTS{e}).(condition{cond}), pulse_x_idx.(pulse_cond)); %return the t-series of responders within that pulse condition
            elseif ResNonRes == 0
                TSERIES.(condition{cond}) = setdiff(pulse_x_idx.(pulse_cond), RESPONDER_TS.(EVENTS{e}).(condition{cond})); %return the t-series of non-responders within that pulse condition
            end
        else
            TSERIES.(condition{cond}) = intersect(pulse_x_idx.(pulse_cond), 1:length(mydata));
        end
        clear pulse_cond RESPONDER_TS
    end
    clear cond temp_idx
    for cond = 1:length(condition)
        FramesAboveThresh_temp = NaN(length(mydata),1);
        for tseries = TSERIES.(condition{cond})
            if ~isempty(mydata(tseries).(condition{cond}))
                FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe;
                tempFT = mydata(tseries).(condition{cond}).(EBL).(EVENTS{e}).dffMatFilter2_Thresh3std(FrameOfIncrease:end);
                if ~isempty(tempFT)
                    TotalPostStimAboveThresh = sum(tempFT); %Count total post-stim frames above threshold
                    if TotalPostStimAboveThresh > 0
                        FirstPostStimFrame = find(tempFT == 1);
                        FirstPostStimFrame = FirstPostStimFrame(1); %Find the first post-stim frame above threshold
                        TotalPostStimFrames = numel(tempFT);
                        
                        FramesAboveThresh_temp(tseries,1) = TotalPostStimAboveThresh / (TotalPostStimFrames - FirstPostStimFrame + 1);
                    end
                end
                clear FrameOfIncrease tempFT TotalPostStimAboveThresh FirstPostStimFrame TotalPostStimFrames
            end
        end
        PropFramesAboveThresh.(EVENTS{e}).(condition{cond}) = FramesAboveThresh_temp;
        clear FramesAboveThresh_temp tseries
    end
end

% 2. Fig. 3f: Violin plots of proportion of post-stim frames above threshold
%    Choose which comparisons to make (using 'condition2plot')
%    Calls function violinplot_STATS.m and cbrewer.m
condition2plot = {'RuBiGABA', 'RuBiGlu'}; % 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' Which conditions should be plotted
PAIRED = [0,0]; %1 if the data within that set is paired or unpaired
EXCLUDE = {'', ''}; %Within that dataset, should any condition be left out?
COLORS = {'RdPu', 'Blues'}; %'RdPu', 'Blues', 'YlGnBu' %'RdPu', 'Blues', 'MgCy', 'MgCy', 'MgCy'
COLOR_Scheme = {'seq', 'seq'};

clear condition
for cp = 1:length(condition2plot)
    if strcmp(condition2plot{cp}, 'AllConditions')
        condition = fieldnames(mydata);
    else
        temp_condition = fieldnames(mydata);
        if ~contains(condition2plot{cp}, '_')
            idx = cell2mat(cellfun(@(x) contains(x, condition2plot{cp}), temp_condition, 'UniformOutput', 0));
        else %to pull out conditions with different strings (for example: All RuBiGABA conditions and NoRuBi_LaserUncagingControl would be used with 'RuBiGABA_NoRuBi'
            idx_1 = cell2mat(cellfun(@(x) contains(x, extractBefore(condition2plot{cp}, '_')), temp_condition, 'UniformOutput', 0));
            idx_2 = cell2mat(cellfun(@(x) contains(x, extractAfter(condition2plot{cp}, '_')), temp_condition, 'UniformOutput', 0));
            idx = find(idx_1 | idx_2);
        end
        if ~isempty(EXCLUDE{cp})
            exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
            idx = setdiff(idx, exclude_idx);
        end
        condition = temp_condition(idx, 1);
        clear temp_condition idx exclude_idx idx_1 idx_2
    end
    
    if ~isempty(COLOR_Scheme{cp})
        cmap = cbrewer(COLOR_Scheme{cp}, COLORS{cp}, 8,'PCHIP');
        if contains(condition2plot{cp}, '_')
            cmap(1,:) = [0.5 0.5 0.5];
        end
        if contains(condition2plot{cp}, '_') || strcmp(COLOR_Scheme{cp}, 'div')
            cmap_num = flip(floor(linspace(1, 8, length(condition))));
        else
            cmap_num = flip(floor(linspace(3, 8, length(condition))));
        end
    else
        cmap = [1 0 1; 0 1 1]; %magneta and cyan
        cmap_num = [1, 2];
    end
    
    for e=1:length(EVENTS)
        struct2plot = struct();
        for cond = 1:length(condition)
            temp = PropFramesAboveThresh.(EVENTS{e}).(condition{cond});
            struct2plot.(condition{cond}) = temp;
            clear temp
        end
        
        if exist('pulse_num', 'var')
            TITLE = {sprintf('%s %s %s %d-100ms pulses', indicator, experiment_type, condition2plot{cp}, pulse_num),...
                'proportion of post-stim frames above threshold (dF/F) after first frame above thresh',...
                sprintf('%s Exclude Ramping Cells: %d, %s', EVENTS{e}, ExcludeRampingCells, EBL),...
                sprintf('%s (%dS PostStim)', WTS, PostStimTime)};
        else
            TITLE = {sprintf('%s %s %s', indicator, experiment_type, condition2plot{cp}),...
                'proportion of post-stim frames above threshold (dF/F) after first frame above thresh',...
                sprintf('%s Exclude Ramping Cells: %d, %s', EVENTS{e}, ExcludeRampingCells, EBL),...
                sprintf('%s (%dS PostStim)', WTS, PostStimTime)};
        end

        YLABEL = {'proportion of post-stim frames above threshold', 'after first post-stim frame above thresh'};
        
        if PAIRED(cp) == 0
            if length(condition) == 2
                STAT = 'RankSum';
            else
                STAT = 'ANOVA';
            end
            
            violinplot_STATS(struct2plot, [], TITLE, YLABEL, cmap, cmap_num, STAT)
            cd(save_dir)
            saveas(gcf, sprintf('%sEx%s_PropFramAbvThreshAfterInitRes_Vio_dFF_ExRamp%dExBl%dS_%s_%s_%dSpost.tif', condition2plot{cp}, EXCLUDE{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime))
            saveas(gcf, sprintf('%sEx%s_PropFramAbvThreshAfterInitRes_Vio_dFF_ExRamp%dExBl%dS_%s_%s_%dSpost.svg', condition2plot{cp}, EXCLUDE{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime))
            close
            cd(start_dir)

        elseif PAIRED(cp) == 1
            Param_PreVPost = [struct2plot.(condition{1}),struct2plot.(condition{2})];
            NUM_TS = TSERIES.(condition{cond});
            
            PlotTSMean_PreVPost(Param_PreVPost, NUM_TS, TITLE, YLABEL, [], ERRORBARS)
            xticklabels({strrep(condition{1}, '_', ' '), strrep(condition{2}, '_', ' ')})
            legend('off')
            cd(save_dir)
            saveas(gcf, sprintf('%sEx%s_PropFramAbvThreshAfterInitRes_dFF_ExRamp%dExBl%dS_%s_%s_%dSpost.tif', condition2plot{cp}, EXCLUDE{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime))
            saveas(gcf, sprintf('%sEx%s_PropFramAbvThreshAfterInitRes_dFF_ExRamp%dExBl%dS_%s_%s_%dSpost.svg', condition2plot{cp}, EXCLUDE{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime))
            close
            cd(start_dir)
            clear Param_PreVPost NUM_TS
        end
    clear struct2plot PlotTS cond TITLE YLABEL STAT cond    
    end
    clear struct2plot PlotTS cond TITLE YLABEL...
        cmap cmap_num STAT cond e condition
end

clear PropFramesAboveThresh EVENTS condition2plot PAIRED EXCLUDE COLORS COLOR_Scheme cp...
    PostStimTime pulse_num TSERIES UseAreaThresh WTS
%%
clear EBL ExcludeBaseline ExcludeRampingCells RampingCellThresh ResNonRes...
    save_dir start_dir WE