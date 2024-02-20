%% Fig. 1e: dose-response curves for percent of field active (avg. change in percent active post-pre)
% Modified from 'PercentFieldAtivePerFrame_MultCond_MC20210817.m'

%  1. Load 'AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat' workspace
%  2. Run'Fig1_PreppingDataStruct.m':
%     Section A: experiment_type = 'BathApp'; method_AgonistEntry = 'MaxCurve'
%     Section B: Calculating parameters/features for population-wide
%     activity
%  3. Run the following code

% Section A. Calculate percent of the field active for a flexible number of conditions
%            for population-wide activity with normalization of choice (not separated by regions)
% Section B. For each condition, combine together tseries for chosen
%            normalization:
%            Trace2Plot: = 'ActivePx'; 
%            Xsec = 1; %Active pixels are calculated for each 1 second
%            period (ends up being/frame with rounding of frame rate)
%            NormBy = '_NormAllActPx'; %Active pixels per frame will be
%            normalized by the number of pixels in the FOV active at some
%            point during the recording
% Section C. Dose-response curves for percent of field active per frame (out of all active px in t-series)
%            For a given time window pre- and post-stim (avg post- avg
%            pre). Calls the function DoseResponseCurve.m: 
%            time_window = 120;

%Michelle Cahill 20240109
%% Section A. Calculate percent of the field active for each t-series of each condition for population-wide activity with normalization of choice (not separated by regions)
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond})) && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
            clear num_evts ActivePxPerFrame ActivePxPerFrameSum xSize ySize ActivePx...
                numActivePx numTotalPx ActivePxPerFrame_NormAllActPx ActivePxPerFrame_NormAllPxFOV
            num_evts = length(mydata(tseries).(condition{cond}).tLocs); % # of evts in that tseries
            ActivePxPerFrame = zeros(num_evts, mydata(tseries).(condition{cond}).totalframes); %matrix of 0s rows:evts columns:frames
            xSize = mydata(tseries).(condition{cond}).res.opts.sz(1)+1; %pixels in x (minus what AQuA takes off of borders) +1 to account for x-vals of 0
            ySize = mydata(tseries).(condition{cond}).res.opts.sz(2)+1; %pixels in y (minus what AQuA takes off of borders) +1 to account for y-vals of 0
            ActivePx = zeros(xSize,ySize); %Keep track of which px in the field of view are active in the tseries
            for evt = 1:num_evts
                if ~isempty(mydata(tseries).(condition{cond}).tLocs{1,evt})
                    clear uv idx num_px ln_idx
                    [uv, ~, idx] = unique(mydata(tseries).(condition{cond}).tLocs{1,evt});
                    num_px = accumarray(idx(:),1)';
                    ActivePxPerFrame(evt, uv) = num_px;
                    %to deal with cases in which either yLocs or xLocs has a index
                    %with 0, add 1 to all x and y coordinates
                    ln_idx = sub2ind([ySize, xSize], mydata(tseries).(condition{cond}).yLocs{1,evt}+1, mydata(tseries).(condition{cond}).xLocs{1,evt}+1);
                    ActivePx(ln_idx) = 1; %active px will have a one to keep track of which px in the field of view are active in the tseries
                end
            end
            ActivePxPerFrameSum = sum(ActivePxPerFrame, 1); %sum across evts to get # of active px per frame
            numActivePx = sum(sum(ActivePx)); %# of total active px in tseries
            numTotalPx = xSize * ySize; %total Px in FOV
            ActivePxPerFrame_NormAllActPx = ActivePxPerFrameSum / numActivePx; %normalizing the #px/frame by the # of active px throughout the tseries
            ActivePxPerFrame_NormAllPxFOV = ActivePxPerFrameSum /  numTotalPx; %normalizing the #px/frame by the # of px in FOV
            
            mydata(tseries).(condition{cond}).ActivePxPerFrame = ActivePxPerFrameSum; %no normalization
            mydata(tseries).(condition{cond}).ActivePxPerFrame_NormAllActPx = ActivePxPerFrame_NormAllActPx; %normalizing the #px/frame by the # of px active at some point during the recording
            mydata(tseries).(condition{cond}).ActivePxPerFrame_NormAllPxFOV = ActivePxPerFrame_NormAllPxFOV; %normalizing the #px/frame by the # of px in FOV
            clear ActivePxPerFrameSum numActivePx numTotalPx ActivePxPerFrame_NormAllActPx ActivePxPerFrame_NormAllPxFOV...
                evt uv idx num_px ln_idx ActivePx ySize xSize ActivePxPerFrame num_evts
        end
    end
end

clear cond tseries
%% Section B. For each condition, combine together tseries for chosen normalization 
% Active pixels with normalization of choice.
Trace2Plot = 'ActivePx'; %'dff', 'ActivePx', 'EventCountPerFrame', 'EventCountPerFrame_NormAllEvts', 'ImageJTraces'
Xsec = 1; %size of desired bins (in seconds)

TreatmentVals = [5 25 50 100]; %what treatment groups do you want to separate by (if empty, include all t-series in each condition) [5 25 50 100]

if strcmp(Trace2Plot, 'ActivePx')
    NormBy = '_NormAllActPx'; % '', '_NormAllActPx', '_NormAllPxFOV' (corresponding to the normalizations from Section A
    Trace2Plot = sprintf('ActivePxPerFrame%s', NormBy);
    if ~isempty(NormBy)
        YL = 'per field active'; %for y-label
    else
        YL = 'num px active'; %for y-label 
    end
end

NormName = strrep(NormBy, '_', '');
if isempty(NormBy)
    NormName = 'NoNorm';
end

% combine together t-series
XsecBin_AllActPx = struct();
for cond = 1:length(condition)
    XsecBin_AllActPx.(condition{cond}) = NaN(length(mydata),2000);
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))  && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
            PerFieldActive = mydata(tseries).(condition{cond}).(Trace2Plot);
            FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe;
            SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
            NumFrames_Xsec = round(Xsec/SecPerFrame); %The number of frames in x seconds, listed above as 5sec (Xsec = 5)
            
            FramesCroppedBeginning = mod(FrameOfIncrease-1, NumFrames_Xsec);
            if FramesCroppedBeginning == 0
                PerFieldActive_Pre = PerFieldActive(1 : FrameOfIncrease-1);
            else
                PerFieldActive_Pre = PerFieldActive(FramesCroppedBeginning+1 : FrameOfIncrease-1);
            end
            
            FramesCroppedEnd = mod((length(PerFieldActive)- FrameOfIncrease), NumFrames_Xsec);
            PerFieldActive_Post = PerFieldActive(FrameOfIncrease : end-(FramesCroppedEnd+1));
            
            if contains(Trace2Plot, 'ActivePx')
                PerFieldActivePerXSec_Pre = sum(reshape(PerFieldActive_Pre, NumFrames_Xsec, []),1);
                PerFieldActivePerXSec_Post = sum(reshape(PerFieldActive_Post, NumFrames_Xsec, []),1);
            else
                PerFieldActivePerXSec_Pre = mean(reshape(PerFieldActive_Pre, NumFrames_Xsec, []),1);
                PerFieldActivePerXSec_Post = mean(reshape(PerFieldActive_Post, NumFrames_Xsec, []),1);
            end
            
            XsecBin_AllActPx.(condition{cond})(tseries, (1000-length(PerFieldActivePerXSec_Pre)+1):1000) = PerFieldActivePerXSec_Pre;
            XsecBin_AllActPx.(condition{cond})(tseries, (1000+1):(1000+length(PerFieldActivePerXSec_Post))) = PerFieldActivePerXSec_Post;
            
            NumXsecBins.(condition{cond}).Pre(tseries) = numel(PerFieldActivePerXSec_Pre);
            NumXsecBins.(condition{cond}).Post(tseries) = numel(PerFieldActivePerXSec_Post);
            
            clear PerFieldActive FrameOfIncrease SecPerFrame NumFrames_Xsec...
                FramesCroppedBeginning PerFieldActive_Pre PerFieldActivePerXSec_Pre...
                FramesCroppedEnd PerFieldActive_Post PerFieldActivePerXSec_Post ax1 ax2
        end
    end
    clear TV concentration_TS_idx
end

clear cond tseries
%% Section C. Dose-response curves for percent of field active per frame (out of all active px in t-series)
%   For a given time window pre- and post-stim (avg post- avg pre)
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig1E';

time_window = 120; %time window (in seconds) taken before the stim and after the stim, [] for all time pre and post
ERRORBARS = 'sem'; %'sem' or 'std'

clear condition
condition = fieldnames(mydata);

Change_PreVPost = struct(); %record the difference post-pre for the average of the feature in the time window given
for cond = 1:length(condition)
    mat2plot = XsecBin_AllActPx.(condition{cond});
    SecPerFrame = mydata(1).(condition{cond}).SecPerFrame;
    if ~isempty(time_window)
        if Xsec == 1
            num_bins = round(time_window / SecPerFrame);
        else
            num_bins = round(time_window / Xsec);
        end
        baseline_temp = mat2plot(:, 1000-num_bins :1000); %baseline bins
        poststim_temp = mat2plot(:, 1001:1001+num_bins); %post-stim bins
    else
        baseline_temp = mat2plot(:, 1:1000); %baseline bins
        poststim_temp = mat2plot(:, 1001:end); %post-stim bins
    end
    
    avg_baseline = nanmean(baseline_temp, 2); %get the average baseline value for each row/t-series/region
    avg_poststim = nanmean(poststim_temp, 2); %get the average post-stim value for each row/t-series/region
    avg_diff = (avg_poststim - avg_baseline) * 100;
    Change_PreVPost.(condition{cond}) = avg_diff;
    clear mat2plot SecPerFrame num_bins baseline_temp poststim_temp avg_baseline avg_poststim avg_diff
end
clear cond

condition2plot = {{'Baclofen', 'tACPD'}}; %choose which conditions to plot together e.g.{{X, XX}, {Y,YY}} will plot 2 graphs with the datasets X&XX and Y&YY
TTEST_Type = {'paired'}; % 'paired' or 'unpaired'
for cp = 1:length(condition2plot)
    CONDITION = condition2plot{cp}; %{'Baclofen', 'tACPD'}; %{'Baclofen_CBX', 'tACPD_CBX'}
    TITLE = {'Cyto-GCaMP 16X Bath App:',sprintf('Percent Field Active %s', NormName),...
        sprintf('%dS pre & post', time_window)};
    YLABEL = '\Delta % field active (post-pre)';
    YLIMS = [-5 30];
    XLABEL = 'concentration (uM)';
    XLIMS = [0 max(TreatmentVals)+10];
    DATACOLORS = []; %{'-mo', '-co'}; %{'-mo', '-co'}
    for cond = 1:length(CONDITION) %determining the line style and color for each condition to be plotted
        c = CONDITION{cond};
        if contains(c, 'CBX')
            tempLINE = '--';
        else
            tempLINE = '-';
        end
        if contains(c, 'Baclofen')
            tempCOLOR = 'm';
        elseif contains(c, 'tACPD')
            tempCOLOR = 'c';
        end
        tempDC = strcat(tempLINE, tempCOLOR, 'o');
        DATACOLORS = [DATACOLORS, {tempDC}];
        clear c tempLINE tempCOLOR tempDC
    end
    DoseResponseCurves(Change_PreVPost, pulse_x_idx, TreatmentVals, CONDITION,...
        TITLE, YLABEL, YLIMS, XLABEL, XLIMS,...
        DATACOLORS, ERRORBARS, 2, 1, TTEST_Type{cp})
    
    WhichDataSets = string(strrep(join(CONDITION), ' ', ''));
    cd(save_dir)
    saveas(gcf, sprintf('DoseResponsePerPxAct%s_%dSPrePost_%s.tif', NormBy, time_window, WhichDataSets))
    saveas(gcf, sprintf('DoseResponsePerPxAct%s_%dSPrePost_%s.svg', NormBy, time_window, WhichDataSets))
    close
    cd(start_dir)
    clear CONDITION TITLE YLABEL YLIMS XLABEL XLIMS DATACOLORS WhichDataSets
end

clear ans BinSize ERRORBARS NormBy NormName NumXsecBins save_dir start_dir time_window...
    Trace2Plot Xsec XsecBin_AllActPx XsecBin_AllActPx_perm YL cp condition2plot...
    Change_PreVPost TTEST_Type cond