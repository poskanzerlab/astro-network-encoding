%% Mean dF/F per 30s time bins from active astrocytes in each t-series
% Extended Data Fig. 1d: time-series traces of average dF/F in 30s windows
% from active astrocytes in each slice. Permutation testing used to
% determine time windows with dF/F values significantly greater than chance. 
% dF/F values calculated using moving 10s baseline window, averaging the 
% lower 50% of values in the window  
%
% Extended Data Fig. 1e: Mean dF/F before and after bath application of
% agonist. Average dF/F pre-stim (300-0s pre-stim?) vs. 30s time window 
% with highest dF/F value for each slice (post-stim peak). Change in
% average dF/F values, comparing the response of each slice to both
% agonists
% Modified from:
%     'FittingSigmoids_ImageJTracesByRegion_MultCond_MC20211202.m
%     'Average_dFF_PreVPostByRegion_MultCond_MC20220121.m (section
%     3a,c,d&e)

%  1. Load 'AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat' workspace
%  2. Run'Fig1_PreppingDataStruct.m':
%     Section A: experiment_type = 'BathApp'; method_AgonistEntry = 'MaxCurve'
%  3. Run the following code:

% A. Calculate dF/F values for each region of each t-series using raw 
%    fluorescence traces from FIJI/ImageJ.
%    method_AgonistEntry = 'MaxCurve'; indicator = 'CytoGCaMP';
%    MovingWindow = 10;
%    Calls the following functions: 
%       'dFF_lower_20211202.m'
%
% B. Extended Data Fig. 1d: plotting dF/F in 30sec bins averaging across
%    t-series/recordings for the chosen agonist and treatment value (concentration). 
%    Comparing experimental mean to permuted means to determine 30s bins with
%    dF/F values significantly greater than chance.
%    Trace2Plot = 'ImageJTraces'
%    AVERAGE_ACROSS = 'TSERIES'
%    XSec = 30; ERRORBARS = 'sem'; NormBy = ''; dFF_T = 'dFF_lower10s'
%    ExcludeBaseline = 30; ExcludeEnd = 30;
%    PerformPerm = 1; perm_num = 10000;
%    Calls the following functions: 
%      'MeanXSecBins_ActiveRegions_Perm.m'
%      'fdr_bh.m', 'MultiLineLabels.m' and 'vline.m'
%
% C. Extended Data Fig. 1e: Mean dF/F before and after bath application of
%     agonist. Run sections A & B before running this section 
%     Calls the following functions:
%       'PlotTSMean_PreVPost.m'
%       'cbrewer.m' 

% Michelle Cahill 20240110
%% A. Calculate dF/F values for each region of each t-series using raw fluorescence traces from FIJI/ImageJ
% Raw fluorescence traces are loaded in data struct in the field
% ImageJByRegion.RawTrace
% Calls the function dFF_lower_20211202.m (alternatively could call
% dFF_lowerXPercent_20220127.m to specify the percent of lowest values from
% moving window to calculate baseline, to change from using lower 50%)
start_dir = pwd;
LoadRawData = 0; %1 if raw fluorescence traces per region from FIJI need to be loaded into exisiting data struct, 0 if traces are already loaded 
if LoadRawData
    data_dir = ''; %Choose the file containing the raw fluorescence .csv files if loading raw fluorescence traces per region from FIJI
end
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1d';
method_AgonistEntry = 'MaxCurve'; %'Threshold', 'MaxCurve'
indicator = 'CytoGCaMP'; %'CytoGCaMP', 'pinkFlamindo'
condition = fieldnames(mydata);

% 1. Load the fluorescence traces by regions (.csv files) for each t-series
%    in the struct if traces have not yet been loaded. Columns are regions and rows are frames
if LoadRawData
    for cond = 1:length(condition)
        for tseries = 1:length(mydata)
            curr_file = mydata(tseries).(condition{cond}).file;
            curr_file = extractAfter(curr_file, 'AQuA\');
            trace_file = strcat(curr_file, '_ZAxisProfile.csv');
            cd(data_dir)
            temp_trace = table2array(readtable(trace_file));
            temp_trace = temp_trace(:,2:end);
            mydata(tseries).(condition{cond}).ImageJByRegion.RawTrace = temp_trace;
            
            figure()
            plot(mydata(tseries).(condition{cond}).ImageJByRegion.RawTrace)
            hold on
            vline(mydata(tseries).(condition{cond}).uncagingframe, 'r')
            title(sprintf('ImageJ Traces By Region: %s TS %d', condition{cond}, tseries))
            xlabel('time (frames)')
            ylabel('average fluorescence')
            
            cd(save_dir)
            saveas(gcf, sprintf('ImageJTracesByRegion_%s_%s_TS%d.tif', indicator, condition{cond}, tseries))
            saveas(gcf, sprintf('ImageJTracesByRegion_%s_%s_TS%d.svg', indicator, condition{cond}, tseries))
            close
            cd(start_dir)
            
            clear curr_file trace_file temp_trace
        end
    end
    clear cond tseries data_dir
end

% 2. Calculate dF/F per region from the raw trace by using a moving 10s 
%    baseline window, averaging the lower 50% of values in the window.
%    Calls the function dFF_lower_20211202.m
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        RawTrace = mydata(tseries).(condition{cond}).ImageJByRegion.RawTrace;
        SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
        MovingWindow = 10;
        dFF_T = sprintf('dFF_lower%ds', MovingWindow);
        mydata(tseries).(condition{cond}).ImageJByRegion.(dFF_T) = dFF_lower_20211202(RawTrace, SecPerFrame, MovingWindow);
        
        figure()
        plot(mydata(tseries).(condition{cond}).ImageJByRegion.(dFF_T))
        hold on
        vline(mydata(tseries).(condition{cond}).uncagingframe, 'r')
        title({sprintf('ImageJ Traces By Region: %s TS %d', condition{cond}, tseries),...
            sprintf('dF/F lower: %dS moving window', MovingWindow)})
        xlabel('time (frames)')
        ylabel('dF/F lower')
        
        cd(save_dir)
        saveas(gcf, sprintf('ImageJTracesByRegion_dFF%dSWind_%s_%s_TS%d.tif', MovingWindow, indicator, condition{cond}, tseries))
        saveas(gcf, sprintf('ImageJTracesByRegion_dFF%dSWind_%s_%s_TS%d.svg', MovingWindow, indicator, condition{cond}, tseries))
        close
        cd(start_dir)
        clear RawTrace SecPerFrame MovingWindow
    end
end

clear cond tseries LoadRawData
%% B. Extended Data Fig. 1d: plotting dF/F in 30sec bins averaging across t-series/recordings for chosen agonist and treatment value (concentration) 
% Comparing experimental mean to permuted mean to determine 30s bins with
% dF/F values significantly greater than chance (one-tailed since test was
% run on absolute values and not change from average baseline). Exclude the
% first and final 30s bins to avoid edge arifacts
% Note: run time for 10,000 permutations for a single condition including 4-6 recordings is ~1 min 
% Calls the following functions: 
%      'MeanXSecBins_ActiveRegions_Perm.m'
%      'fdr_bh.m', 'MultiLineLabels.m' and 'vline.m'

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1d';

Trace2Plot = 'ImageJTraces'; %'dff', 'ActivePx', 'EventCountPerFrame', 'EventCountPerFrame_NormAllEvts', 'ImageJTraces'
AVERAGE_ACROSS = 'TSERIES'; % 'TSERIES', 'REGIONS'
Xsec = 30; %size of desired bins (in seconds)
BinSize = sprintf('TimeBins%dS', Xsec);
ERRORBARS = 'sem';

TreatmentVals = [5 25 50 100]; %what treatment groups do you want to separate by (if empty, include all t-series in each condition) [5 25 50 100]

PerformPerm = 1; %0: don't perform any permutation test 1: perform a permutation test
perm_num = 10000; %How many permutations should be carried out for each condition?

if strcmp(Trace2Plot, 'ImageJTraces')
    NormBy = ''; % '', '_NormAllActPx', '_NormAllPxFOV'
    Trace2Plot = dFF_T; %which dF/F normalization should be used? dFF_T should be a string that matches a field within the field ImageJByRegion
    ExcludeBaseline = 30; %number of seconds at the beginning of baseline to exclude, 0 to include all baseline period
    ExcludeEnd = 30; %number of seconds to exclude at the end 
    YL = 'dF/F'; %'dF/F lower'; %for y-label
elseif strcmp(Trace2Plot, 'ActivePx')
    NormBy = ''; % '', '_NormAllActPx', '_NormAllPxFOV'
    Trace2Plot = sprintf('ActivePxPerFrame%s', NormBy);
    if ~isempty(NormBy)
        YL = 'per field active'; %for y-label
    else
        YL = 'num px active'; %for y-label 
    end
elseif strcmp(Trace2Plot, 'dff')
    NormBy = ''; % '', '_NormAllActPx', '_NormAllPxFOV'
    Trace2Plot = 'dffMatFilter2';
    ExcludeBaseline = 0; %number of seconds at the beginning of baseline to exclude, 0 to include all baseline period
    EBL = sprintf('Exclude%dSBL', ExcludeBaseline);
    YL = 'dF/F'; %for y-label
else
    NormBy = ''; % '', '_NormAllActPx', '_NormAllPxFOV'
    Trace2Plot = Trace2Plot;
    YL = Trace2Plot; %for y-label
end

if ~exist('XsecBin_perm')
    XsecBin_perm = struct();
end

NormName = strrep(NormBy, '_', '');
if isempty(NormBy)
    NormName = 'NoNorm';
end

% 1. Calculate average dF/F value/30s for each t-series/recording using active 
%    astrocytes in each FOV.
%    a. Combine together t-series (in XsecBin_AllActPx) and perform a permutation test if desired
%    b. Calculate average baseline dF/F across pre-stim 30s time windows and
%    post-stim peak for each t-series/recording (this will be used in
%    Section C)

XsecBin_AllActPx = struct();
AvgPre_PeakPost = struct();
for cond = 1:length(condition)
    for TV = TreatmentVals
        concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
        temp_TS_idx = pulse_x_idx.(concentration_TS_idx);
        %         temp_TS_idx = 1:length(mydata); %pulse_x_idx.WT_100;
        
        if strcmp(AVERAGE_ACROSS, 'TSERIES')
            XsecBin_AllActPx.(concentration_TS_idx) = NaN(numel(temp_TS_idx),2000);
            AvgPre_PeakPost.(concentration_TS_idx) = NaN(numel(temp_TS_idx),2); %column 1: average across baseline bins, column 2: value of peak post-stim bin
        elseif strcmp(AVERAGE_ACROSS, 'REGIONS')
            XsecBin_AllActPx.(concentration_TS_idx) = NaN(NUM_REGIONS.(concentration_TS_idx),2000);
            AvgPre_PeakPost.(concentration_TS_idx) = NaN(NUM_REGIONS.(concentration_TS_idx),2);
            
            %             XsecBin_AllActPx.(concentration_TS_idx) = NaN(sum(ActiveCells.(condition{cond})),2000);
            %             AvgPre_PeakPost.(concentration_TS_idx) = NaN(sum(ActiveCells.(condition{cond})),2);
            
            %             XsecBin_AllActPx.(concentration_TS_idx) = NaN(147,2000);
            %             AvgPre_PeakPost.(concentration_TS_idx) = NaN(147,2);
        end
        Region_count = 0;
        for tseries = temp_TS_idx
            if ~isempty(mydata(tseries).(condition{cond}))  && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
                %             if ~isempty(mydata(tseries).(condition{cond}))
                if contains(Trace2Plot, 'dff')
                    PerFieldActive = mydata(tseries).(condition{cond}).(EBL).(Trace2Plot);
                elseif contains(Trace2Plot, dFF_T)
                    if strcmp(indicator, 'CytoGCaMP')
                        % Finding cells with at lest 1 AQuA event within across all
                        % conditions in that FOV
                        CellsWithActivityCond1 = sum(mydata(tseries).(condition{1}).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                        CellsWithActivityCond2 = sum(mydata(tseries).(condition{2}).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                        CellsWithActivity = sum([CellsWithActivityCond1; CellsWithActivityCond2], 1);
                        % Finding cells with at lest 1 AQuA event within this
                        % condition
                        %                         CellsWithActivity = sum(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                        ACTIVE_REGIONS = find(CellsWithActivity);
                        
                        if strcmp(AVERAGE_ACROSS, 'TSERIES')
                            PerFieldActive = mean(mydata(tseries).(condition{cond}).ImageJByRegion.(Trace2Plot)(:,ACTIVE_REGIONS),2, 'omitnan')';
                        elseif strcmp(AVERAGE_ACROSS, 'REGIONS')
                            PerFieldActive = mydata(tseries).(condition{cond}).ImageJByRegion.(Trace2Plot)(:,ACTIVE_REGIONS)';
                        end
                        clear CellsWithActivity ACTIVE_REGIONS CellsWithActivityCond1 CellsWithActivityCond2
                        
                    elseif strcmp(indicator, 'pinkFlamindo')
                        if strcmp(AVERAGE_ACROSS, 'TSERIES')
                            PerFieldActive = mean(mydata(tseries).(condition{cond}).ImageJByRegion.(Trace2Plot),2, 'omitnan')';
                        elseif strcmp(AVERAGE_ACROSS, 'REGIONS')
                            PerFieldActive = mydata(tseries).(condition{cond}).ImageJByRegion.(Trace2Plot)';
                        end
                    end
 
                else
                    PerFieldActive = mydata(tseries).(condition{cond}).(Trace2Plot);
                end
                FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe;
                SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
                NumFrames_Xsec = round(Xsec/SecPerFrame); %The number of frames in x seconds, listed above as 5sec (Xsec = 5)
                EXCLUDEBINS_START = ceil(ExcludeBaseline/Xsec); %number of Xsec bins to exclude from the start of all recordings (used at this point for calculating baseline average)             
                EXCLUDEBINS_END = ceil(ExcludeEnd/Xsec); %number of Xsec bins to exclude from the end of all recordings (used at this point for calculating post-stim peak frame)
                
                for reg = 1:size(PerFieldActive, 1)
                    Region_count = Region_count + 1;
                    FramesCroppedBeginning = mod(FrameOfIncrease-1, NumFrames_Xsec);
                    if FramesCroppedBeginning == 0
                        PerFieldActive_Pre = PerFieldActive(reg, 1 : FrameOfIncrease-1);
                    else
                        PerFieldActive_Pre = PerFieldActive(reg, FramesCroppedBeginning+1 : FrameOfIncrease-1);
                    end
                    
                    FramesCroppedEnd = mod((length(PerFieldActive)- FrameOfIncrease), NumFrames_Xsec);
                    PerFieldActive_Post = PerFieldActive(reg, FrameOfIncrease : end-(FramesCroppedEnd+1));
                    
                    if contains(Trace2Plot, 'ActivePx')
                        PerFieldActivePerXSec_Pre = sum(reshape(PerFieldActive_Pre, NumFrames_Xsec, []),1);
                        PerFieldActivePerXSec_Post = sum(reshape(PerFieldActive_Post, NumFrames_Xsec, []),1);
                    else
                        PerFieldActivePerXSec_Pre = mean(reshape(PerFieldActive_Pre, NumFrames_Xsec, []),1);
                        PerFieldActivePerXSec_Post = mean(reshape(PerFieldActive_Post, NumFrames_Xsec, []),1);
                    end
%                     PerFieldActivePerXSec_Pre = PerFieldActivePerXSec_Pre(1+EXCLUDEBINS_START:end);
%                     PerFieldActivePerXSec_Post = PerFieldActivePerXSec_Post(1:end-EXCLUDEBINS_END);
                    
                    XsecBin_AllActPx.(concentration_TS_idx)(Region_count, (1000-length(PerFieldActivePerXSec_Pre)+1):1000) = PerFieldActivePerXSec_Pre;
                    XsecBin_AllActPx.(concentration_TS_idx)(Region_count, (1000+1):(1000+length(PerFieldActivePerXSec_Post))) = PerFieldActivePerXSec_Post;
                    
                    temp_pre = mean(PerFieldActivePerXSec_Pre(1+EXCLUDEBINS_START:end));
                    temp_post = max(PerFieldActivePerXSec_Post(1:end-EXCLUDEBINS_END));
                    AvgPre_PeakPost.(concentration_TS_idx)(Region_count, 1) = temp_pre;
                    AvgPre_PeakPost.(concentration_TS_idx)(Region_count, 2) = temp_post;
                    
                    clear FramesCroppedBeginning PerFieldActive_Pre...
                        FramesCroppedEnd PerFieldActive_Post PerFieldActivePerXSec_Pre...
                        PerFieldActivePerXSec_Post temp_pre temp_post
                end
                %             NumXsecBins.(condition{cond}).Pre(tseries) = numel(PerFieldActivePerXSec_Pre);
                %             NumXsecBins.(condition{cond}).Post(tseries) = numel(PerFieldActivePerXSec_Post);
                
                clear PerFieldActive FrameOfIncrease SecPerFrame NumFrames_Xsec...
                    ax1 ax2 EXCLUDEBINS_START EXCLUDEBINS_END reg
            end
        end
    end
    
    %If desired, perform a permutation test perm_num times, shuffling the
    %frames/Xsec bins to create a null_distribution for each treatment
    %value
    if PerformPerm == 1
        if ~isempty(TreatmentVals) %For permuation tests per treatment value (concentration or pulse number etc.)
            for TV = TreatmentVals
                concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
                if ~isfield(XsecBin_perm, BinSize) || ~isfield(XsecBin_perm.(BinSize), NormName) ||...
                        ~isfield(XsecBin_perm.(BinSize).(NormName), concentration_TS_idx) %checking to make sure this specific permutation test hasn't already been run
                    sprintf('START: %s', datetime)
                    %                     [XsecBin_perm.(BinSize).(NormName).(concentration_TS_idx)] = MeanXSecBins_ActiveRegions_Perm(mydata,condition, condition{cond}, NUM_REGIONS.(concentration_TS_idx), pulse_x_idx.(concentration_TS_idx),...
                    %                         Xsec, NormBy, perm_num, 2000, Trace2Plot, AVERAGE_ACROSS, ExcludeBaseline, ExcludeEnd);
                    %                     sprintf('END: %s', datetime)
                    [XsecBin_perm.(BinSize).(NormName).(concentration_TS_idx)] = MeanXSecBins_ActiveRegions_Perm(mydata,condition, condition{cond}, [], pulse_x_idx.(concentration_TS_idx),...
                        Xsec, NormBy, perm_num, 2000, Trace2Plot, AVERAGE_ACROSS, ExcludeBaseline, ExcludeEnd);
                    sprintf('END: %s', datetime)
                end
            end
        else %To include all tseries of a condition in the permuation test
%             if ~isfield(XsecBin_perm, BinSize) || ~isfield(XsecBin_perm.(BinSize), NormName) ||...
%                     ~isfield(XsecBin_perm.(BinSize).(NormName), concentration_TS_idx) %checking to make sure this specific permutation test hasn't already been run
%                 sprintf('START: %s', datetime)
%                 [XsecBin_perm.(BinSize).(NormName).(condition{cond})] = MeanXSecBins_ActiveRegions_Perm(mydata,condition{cond}, 1:length(mydata),...
%                     Xsec, NormBy, perm_num, 2000);
%                 sprintf('END: %s', datetime)
%             end
            
        end
    end
    clear TV concentration_TS_idx
end
clear cond tseries

% 2. Plot mean+-error bars for each condition separated by
% treatment value (concentration, pulse number, etc.) indicating which bins
% are significantly increased above chance if a permutation test was run
CMAP.Baclofen = [1 0 1]; %magenta
CMAP.tACPD = [0 1 1]; %cyan
if PerformPerm == 1
    if ~isempty(TreatmentVals)
        pval_PerFrame = struct();
        for cond = 1:length(condition)
            for TV = TreatmentVals
                concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
                temp_TS_idx = pulse_x_idx.(concentration_TS_idx);
                mat2plot = XsecBin_AllActPx.(concentration_TS_idx);
                ValBin = find(sum(~isnan(mat2plot),1)>0); %Find only the bins with values, not nan
                EXCLUDEBINS_START = ceil(ExcludeBaseline/Xsec); %number of Xsec bins to exclude from the start of all recordings
                EXCLUDEBINS_END = ceil(ExcludeEnd/Xsec); %number of Xsec bins to exclude from the end of all recordings
                first_bin = ValBin(1) + EXCLUDEBINS_START;
                last_bin = ValBin(end) - EXCLUDEBINS_END;
                SecPerFrame = mydata(temp_TS_idx(1)).(condition{cond}).SecPerFrame;
                
                mat2plot_cropped = mat2plot(:, first_bin:last_bin);
                
                AlignedOnsetFrame = (size(XsecBin_AllActPx.(concentration_TS_idx), 2)/2) + 1;
                if Xsec == 1
                    XVALS = 1:numel(first_bin:last_bin);
                    XVALS = (XVALS - numel(first_bin:AlignedOnsetFrame)) * mydata(1).(condition{cond}).SecPerFrame;
                else
                    XVALS = 1:numel(first_bin:last_bin);
                    XVALS = (XVALS - numel(first_bin:AlignedOnsetFrame)) * Xsec;
                end
                
                exp_mean = nanmean(mat2plot_cropped, 1);
                
                perm_mean = XsecBin_perm.(BinSize).(NormName).(concentration_TS_idx)(:, first_bin:last_bin);
                
                PermGreaterExp = perm_mean >= exp_mean;
                pval_PerFrame.(concentration_TS_idx) = (sum(PermGreaterExp, 1)+1) ./ (perm_num+1);
                % correct for multiple comparisons using Benjamini Yekutieli FDR
                % method with FDR of 5%
                [h.(concentration_TS_idx), ~, ~, adj_pval_PerFrame.(concentration_TS_idx)] = fdr_bh(pval_PerFrame.(concentration_TS_idx), 0.05, 'dep');
                sig_frames.(concentration_TS_idx) = find(h.(concentration_TS_idx));
                
                XTICKLABELS_PVAL = cell(3, length(XVALS));
                for BINS = 1:length(XVALS)
                    XTICKLABELS_PVAL{1, BINS} = sprintf('%d ', XVALS(BINS));
                    XTICKLABELS_PVAL{2, BINS} = sprintf('p=%.4f ', pval_PerFrame.(concentration_TS_idx)(BINS));
                    XTICKLABELS_PVAL{3, BINS} = sprintf('adj p=%.4f ', adj_pval_PerFrame.(concentration_TS_idx)(BINS));
                    clear idx
                end
                clear BINS
                tickLabels = MultiLineLabels(XTICKLABELS_PVAL); %allows for multiline xlabels
                
                figure() %For each condition plot the average dF/F across regions or t-series per x-sec bins
                if strcmp(ERRORBARS, 'std') %plot permutation mean in light grey
                    mean_plot = nanmean(perm_mean, 1);
                    errorbar(XVALS, mean_plot, nanstd(perm_mean,1), '-o',...
                        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none',...
                        'Color', [0.5 0.5 0.5])
                elseif strcmp(ERRORBARS, 'sem')
                    mean_plot = nanmean(perm_mean, 1);
                    SEM = nanstd(perm_mean,1) / (sqrt(size(perm_mean, 1)-1));
                    errorbar(XVALS, mean_plot, SEM, '-o',...
                        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none',...
                        'Color', [0.5 0.5 0.5])
                end
                hold on
                clear mean_plot SEM
                
                if strcmp(ERRORBARS, 'std') %plot experimental mean in magenta or cyan
                    mean_plot = nanmean(mat2plot_cropped, 1);
                    errorbar(XVALS, mean_plot, nanstd(mat2plot_cropped,1), '-o',...
                        'MarkerFaceColor', CMAP.(condition{cond}), 'MarkerEdgeColor', 'none',...
                        'Color', CMAP.(condition{cond}))
                elseif strcmp(ERRORBARS, 'sem')
                    mean_plot = nanmean(mat2plot_cropped, 1);
                    SEM = nanstd(mat2plot_cropped,1) / (sqrt(size(mat2plot_cropped, 1)-1));
                    errorbar(XVALS, mean_plot, SEM, '-o',...
                        'MarkerFaceColor', CMAP.(condition{cond}), 'MarkerEdgeColor', 'none',...
                        'Color', CMAP.(condition{cond}))
                end
                clear mean_plot SEM

                %plot the frames that are significantly above shuffled data
                MAX_VAL = max(max(nanmean(mat2plot_cropped, 1)), max(nanmean(perm_mean, 1)));
                BUFFER_OrderMag = floor(log(abs(MAX_VAL))./log(10));
                BUFFER = 10^BUFFER_OrderMag;
                
                Pval_placement = MAX_VAL + BUFFER;
                clear MAX_VAL BUFFER_OrderMag BUFFER
                if Xsec == 1
                    SF = (sig_frames.(concentration_TS_idx) - numel(first_bin:AlignedOnsetFrame)) * SecPerFrame;
                else
                    SF = (sig_frames.(concentration_TS_idx) - numel(first_bin:AlignedOnsetFrame)) * Xsec;
                end
                scatter(SF, ones(1, length(sig_frames.(concentration_TS_idx)))* Pval_placement,  'k', 'fill')
                
                vline(-0.5*Xsec, '--k')
                if Xsec == 1
                    xlabel('time (s)')
                else
                    xlabel(sprintf('time (%dS bins)', Xsec))
                end
                
                ylabel(sprintf('mean %s +- %s', YL, ERRORBARS))
%                 if contains(Trace2Plot, 'ActivePx')
%                     if ~isempty(NormBy)
%                         ylabel(sprintf('mean per field active +- %s', ERRORBARS))
%                     else
%                         ylabel(sprintf('mean num px active +- %s', ERRORBARS))
%                     end
%                 elseif contains(Trace2Plot, 'dff')
%                     ylabel(sprintf('mean dF/F +- %s', ERRORBARS))
%                 end
                title({sprintf('%s %s %s (%duM):', indicator, experiment_type, strrep(condition{cond}, '_', ' '), TV), sprintf('%s: %dS bins', strrep(Trace2Plot, '_', ' '), Xsec),...
                    sprintf('%d permutations, n=%d %s', perm_num, size(mat2plot,1), lower(AVERAGE_ACROSS))})
                xticks(XVALS)
                xlim([-300 270])
                xticklabels(tickLabels)
                set(gcf,'position', [100, 200, 2000, 400])
                
                cd(save_dir)
                saveas(gcf, sprintf('%s%duM_Perm%s_%dS_mean%s_%s.tif', condition{cond}, TV, Trace2Plot, Xsec, ERRORBARS, AVERAGE_ACROSS))
                saveas(gcf, sprintf('%s%duM_Perm%s_%dS_mean%s_%s.svg', condition{cond}, TV, Trace2Plot, Xsec, ERRORBARS, AVERAGE_ACROSS))
                close
                cd(start_dir)
                clear mat2plot ValBin first_bin last_bin mat2plot_cropped exp_mean perm_mean...
                    PermGreaterExp exclude_bins XLIMS SecPerFrame concentration_TS_idx temp_TS_idx...
                    XVALS AlignedOnsetFrame SF adj_pval_PerFrame h pval_PerFrame Pval_placement...
                    sig_frames EXCLUDEBINS_START EXCLUDEBINS_END tickLabels XTICKLABELS_PVAL
            end
            clear TV
        end
        clear cond
    end
end

clear BinSize CMAP ExcludeBaseline ExcludeEnd NormBy NormName...
    PerformPerm perm_num Region_count save_dir start_dir Trace2Plot XsecBin_AllActPx YL...
    ans TreatmentVals
%% C. Extended Data Fig. 1e: Mean dF/F before and after bath application of
%     agonist. Run sections A & B before running this section 
% 1. Average dF/F pre-stim (300-0s pre-stim) vs. 30s time window 
%    with highest dF/F value for each slice (post-stim peak), calculated above in Section B. 
% 2. Change in average dF/F values (post peak - avg baseline), comparing the
%    response of each slice to both agonists
% Calls the following functions:
%     'PlotTSMean_PreVPost.m'
%     'cbrewer.m' 

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1e';

TreatmentVals = [5 25 50 100];

% 1. Plot baseline mean vs post-stim peak for each condition separated by
% treatment value, running a paired t-test (Extended Data Fig. 1e left and
% center)
for cond = 1:length(condition)
    for TV = TreatmentVals
        concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
        Param_PreVPost = AvgPre_PeakPost.(concentration_TS_idx);
        NUM_TS = 1:size(Param_PreVPost, 1);
        YLABEL = 'mean dF/F';
        
        %plotting each t-series mean pre and mean post and the group mean
        PlotTSMean_PreVPost( Param_PreVPost, NUM_TS,...
            {sprintf('Average dF/F: %s %duM', strrep(condition{cond}, '_', ' '), TV),...
            sprintf('avg %dS bins pre vs peak post bin, n=%d %s', Xsec, size(Param_PreVPost, 1), lower(AVERAGE_ACROSS))},...
            YLABEL, [], ERRORBARS )
        
        cd(save_dir)
        saveas(gcf, sprintf('%s_%s%duMAvgdFFByTseries_ImageJTraces_Avg%dSPrePeakPost.tif', indicator, condition{cond}, TV, Xsec))
        saveas(gcf, sprintf('%s_%s%duMAvgdFFByTseries_ImageJTraces_Avg%dSPrePeakPost.svg', indicator, condition{cond}, TV, Xsec))
        close
        cd(start_dir)
        
        MEAN_VALS_dFF_BLvPeak.(AVERAGE_ACROSS).(concentration_TS_idx).mean = mean(Param_PreVPost,1,'omitnan');
        if strcmp(ERRORBARS, 'std')
            MEAN_VALS_dFF_BLvPeak.(AVERAGE_ACROSS).(concentration_TS_idx).(ERRORBARS) = std(Param_PreVPost,1,'omitnan');
        elseif strcmp(ERRORBARS, 'sem')
            MEAN_VALS_dFF_BLvPeak.(AVERAGE_ACROSS).(concentration_TS_idx).(ERRORBARS) = std(Param_PreVPost,1,'omitnan')/(sqrt(size(Param_PreVPost, 1)-1));
        end
        
        clear Param_PreVPost YLABEL NUM_TS concentration_TS_idx
    end
    clear TV
end
clear cond

% 2. Plot change from baseline (post-stim peak bin - mean baseline) 
% between conditions run on the same FOV/slice, running a paired t-test
% (Extended Data Fig. 1e, right)
condition2plot = {'AllConditions'}; %'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions', 'RuBiGABA_NoRuBi', 'RuBiGlu_NoRuBi', 'WT' Which conditions should be plotted
COLORS = {'MgCy'}; %'RdPu', 'Blues', 'YlGnBu' %'RdPu', 'Blues', 'MgCy', 'MgCy', 'MgCy'
COLOR_Scheme = {''}; %'seq', 'seq', '', '', ''

for TV = TreatmentVals
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
                idx = idx_1 | idx_2;
            end
            condition = temp_condition(idx, 1);
            clear temp_condition idx idx_1 idx_2
        end
        
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
        
        change_Cond1VCond2 = [];
        for cond = 1:length(condition)
            concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
            Temp_PreVPost = AvgPre_PeakPost.(concentration_TS_idx);
            Temp_Change = Temp_PreVPost(:,2) -  Temp_PreVPost(:,1);
            change_Cond1VCond2 = [change_Cond1VCond2, Temp_Change];
            clear concentration_TS_idx Temp_PreVPost Temp_Change
        end
        
        Param_PreVPost = change_Cond1VCond2;
        NUM_TS = 1:size(Param_PreVPost, 1);
        YLABEL = 'change dF/F (avg pre - peak post)';
        
        %plotting each t-series mean pre and mean post and the group mean
        PlotTSMean_PreVPost( Param_PreVPost, NUM_TS,...
            {sprintf('Change in dF/F: %s & %s %duM', strrep(condition{1}, '_', ' '), strrep(condition{2}, '_', ' '), TV),...
            sprintf('peak post bin - avg %dS bins pre, n=%d %s', Xsec, size(Param_PreVPost, 1), lower(AVERAGE_ACROSS))},...
            YLABEL, [], ERRORBARS )
        xticklabels(condition)
        
        cd(save_dir)
        saveas(gcf, sprintf('%s_%duMAvgdFFByTseries_ImageJTraces_ChangePeakPostAvg%dSPre.tif', indicator, TV, Xsec))
        saveas(gcf, sprintf('%s_%duMAvgdFFByTseries_ImageJTraces_ChangePeakPostAvg%dSPre.svg', indicator, TV, Xsec))
        close
        cd(start_dir)
        
        temp_TV = sprintf('TV%d', TV);
        MEAN_VALS_ChangedFF_BetweenCond.(AVERAGE_ACROSS).(temp_TV).mean = mean(Param_PreVPost,1,'omitnan');
        if strcmp(ERRORBARS, 'std')
            MEAN_VALS_ChangedFF_BetweenCond.(AVERAGE_ACROSS).(temp_TV).(ERRORBARS) = std(Param_PreVPost,1,'omitnan');
        elseif strcmp(ERRORBARS, 'sem')
            MEAN_VALS_ChangedFF_BetweenCond.(AVERAGE_ACROSS).(temp_TV).(ERRORBARS) = std(Param_PreVPost,1,'omitnan')/(sqrt(size(Param_PreVPost, 1)-1));
        end
        
        clear Param_PreVPost YLABEL NUM_TS concentration_TS_idx change_Cond1VCond2 cmap cmap_num temp_TV
    end
    clear cp
end

clear cond condition2plot COLORS COLOR_Scheme TV Xsec AVERAGE_ACROSS ERRORBARS...
    AvgPre_PeakPost dFF_T XsecBin_perm save_dir start_dir MEAN_VALS_ChangedFF_BetweenCond...
    MEAN_VALS_dFF_BLvPeak