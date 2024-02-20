%% Extended Data Fig. 2b, Ext. Data Fig. 3g-h & Ext. Data Fig. 7d: 
%  Change in event number/30s in the chosen cells/regions following uncaging. 
%  Comparing the event frequency change to shuffled data in each post-stim 
%  time bin to determine which time bins are siginifcantly different from chance.
%
%  Modified from PlotEvts_XSecStatBinsPostMinusPreAvgBL_MC20220317.m
% 
%  For Extended Data Fig. 3g:
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%     Section E: time_window = 30;
%  3. Run the following code:
%     A. Aggregate together the event number/30s for neighboring cells 
%     for all t-series/recordings in each condition. For each t-series/recording,
%     calculate the average change from baseline in each post-stim time bin 
%     Perform a permutation test to determine post-stim bins significantly
%     different from chance:
%          time_window = 30; param_name = 'Number'; WhichCells = 2;
%          UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh =
%          'pval_1'; ExcludeBaseline = 60; PerformPerm = 1; perm_num =
%          10000; pulse_num = 10; PlotIndTS = 0;
%    B. Plot change in event frequency from baseline for each condition vs 
%       permuted data in 30s post-stim time bins:
%           ERRORBARS = 'sem'; plotchange = [1]; AltXTicks = 0; TimeS_BeforeAfterUncaging isn't relevant
%
%  For Extended Data Fig. 2b & Extended Data Fig. 3h:
%  1. Load 'AQuA_CytoGCaMP_2PUncaging_WTRecAntaLaserCtrl_ExtDataFig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'CGP', 'LY', 'LaserUncagingControl'};
%     Section C & D as is
%     Section E: time_window = 30;
%  3. Run the following sections of code:
%
%     For Ext. Data Fig. 2b:
%         A. Aggregate together the event number/30s for uncaginging cells 
%            for all t-series/recordings in each condition. For each t-series/recording,
%            calculate the average change from baseline in each post-stim time bin 
%            Perform a permutation test to determine post-stim bins significantly
%            different from chance:
%                time_window = 30; param_name = 'Number'; WhichCells = 1;
%                UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh =
%                'pval_1'; ExcludeBaseline = 60; PerformPerm = 1; perm_num =
%                10000; pulse_num = 10; PlotIndTS = 0;
%          B. Plot change in event frequency from baseline for each condition vs 
%             permuted data in 30s post-stim time bins:
%             ERRORBARS = 'sem'; plotchange = [1]; AltXTicks = 0; TimeS_BeforeAfterUncaging isn't relevant
%
%     For. Ext Data Fig. 3h:
%         A. Aggregate together the event number/30s for neighboring cells 
%            for all t-series/recordings in each condition. For each t-series/recording,
%            calculate the average change from baseline in each post-stim time bin 
%            Perform a permutation test to determine post-stim bins significantly
%            different from chance:
%                time_window = 30; param_name = 'Number'; WhichCells = 2;
%                UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh =
%                'pval_1'; ExcludeBaseline = 60; PerformPerm = 1; perm_num =
%                10000; pulse_num = 10; PlotIndTS = 0;
%          B. Plot change in event frequency from baseline for each condition vs 
%             permuted data in 30s post-stim time bins:
%             ERRORBARS = 'sem'; plotchange = [1]; AltXTicks = 0; TimeS_BeforeAfterUncaging isn't relevant
%
%  For Extended Data Fig. 7d (top):
%  1. Load 'AQuA_CytoGCaMP_2PMultiRoundGluUncaging_70AU_ExtDataFig7.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B, C & D: as is
%  3. Run ExtendedDataFig7_PreppingDataStruct.m: time_window = 30;
%  4. Run the following sections of code:
%         A. Aggregate together the event number/30s for neighboring cells 
%            for all t-series/recordings in each condition. For each t-series/recording,
%            calculate the average change from baseline in each post-stim time bin 
%            Perform a permutation test to determine post-stim bins significantly
%            different from chance:
%                time_window = 30; param_name = 'Number'; WhichCells = 2;
%                UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1';
%                ExcludeBaseline = 30; (in this dataset ExcludeBaseline = 30 corresponds to excluding the first 60s of baseline, since
%                the first 30s time bin has been cropped in ExtendedDataFig7_PreppingDataStruct.m)
%                PerformPerm = 1; perm_num = 10000; pulse_num = 10; PlotIndTS = 0;
%          B. Plot change in event frequency from baseline for each condition vs 
%             permuted data in 30s post-stim time bins:
%             ERRORBARS = 'sem'; plotchange = [1]; AltXTicks = 1;  TimeS_BeforeAfterUncaging = [-120 570];
%
%  For Extended Data Fig. 7d (bottom):
%  1. Load 'AQuA_CytoGCaMP_2PMultiRoundGluUncagingCtrl_25AU_ExtDataFig7.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B & D: as is
%  3. Run ExtendedDataFig7_PreppingDataStruct.m: time_window = 30;
%  4. Run the following sections of code:
%         A. Aggregate together the event number/30s for neighboring cells 
%            for all t-series/recordings in each condition. For each t-series/recording,
%            calculate the average change from baseline in each post-stim time bin 
%            Perform a permutation test to determine post-stim bins significantly
%            different from chance:
%                time_window = 30; param_name = 'Number'; WhichCells = 2;
%                UseAreaThresh = 0; ExcludeRampingCells = 0; (as some recordings in this dataset don't have ramping cell info) 
%                RampingCellThresh = 'pval_1'; (not relevant)
%                ExcludeBaseline = 30; (in this dataset ExcludeBaseline = 30 corresponds to excluding the first 60s of baseline, since
%                the first 30s time bin has been cropped in ExtendedDataFig7_PreppingDataStruct.m)
%                PerformPerm = 1; perm_num = 10000; pulse_num = 10; PlotIndTS = 0;
%          B. Plot change in event frequency from baseline for each condition vs 
%             permuted data in 30s post-stim time bins:
%             ERRORBARS = 'sem'; plotchange = [1]; AltXTicks = 1;  TimeS_BeforeAfterUncaging = [-120 570];
%
% Michelle Cahill 20240116
%% A. Aggregate together the chosen event feature for the chosen regions 
%     for chosen t-series/recordings in each condition. For each t-series/recording,
%     calculate the average change from baseline in each post-stim time bin 
%     Perform a permutation test to determine post-stim bins significantly
%     different from chance
%  Calls the function 'MeanFt_XSecStatBins_PermEventOnsetByRegion.m'
%
% For Ext. Data Fig. 2b: time_window = 30; param_name = 'Number'; WhichCells = 1;
%                UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh =
%                'pval_1'; ExcludeBaseline = 60; PerformPerm = 1; perm_num =
%                10000; pulse_num = 10; PlotIndTS = 0;
%
% For Ext. Data Fig. 3g-h: time_window = 30; param_name = 'Number'; WhichCells = 2;
%          UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh =
%          'pval_1'; ExcludeBaseline = 60; PerformPerm = 1; perm_num =
%          10000; pulse_num = 10; PlotIndTS = 0;
%
% For Ext. Data Fig. 7d (top): time_window = 30; param_name = 'Number'; WhichCells = 2;
%          UseAreaThresh = 0; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1';
%          ExcludeBaseline = 30; (in this dataset ExcludeBaseline = 30 corresponds to excluding the first 60s of baseline, since
%          the first 30s time bin has been cropped in ExtendedDataFig7_PreppingDataStruct.m)
%          PerformPerm = 1; perm_num = 10000; pulse_num = 10; PlotIndTS = 0;
%
% For Ext. Data Fig. 7d (bottom): time_window = 30; param_name = 'Number'; WhichCells = 2;
%          UseAreaThresh = 0; ExcludeRampingCells = 0; (as some recordings in this dataset don't have ramping cell info) 
%          RampingCellThresh = 'pval_1'; (not relevant)
%          ExcludeBaseline = 30; (in this dataset ExcludeBaseline = 30 corresponds to excluding the first 60s of baseline, since
%          the first 30s time bin has been cropped in ExtendedDataFig7_PreppingDataStruct.m)
%          PerformPerm = 1; perm_num = 10000; pulse_num = 10; PlotIndTS = 0;

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig7d\bottom_RampingCellsNotExcluded';

time_window = 30; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

param_name = 'Number'; % What parameter to plot? 'Number', 'Area', 'Duration' or 'Amplitude'

WhichCells = 2; %1=uncaging cell, 2=all neighbors combined
PlotIndTS = 0; %1 = plot raw counts for each cell by tseries, 0 = skip this step

UseAreaThresh = 0; %1 to plot only responding or non-responding cells/t-series, 0 to plot all tseries/recordings in a condition
ResNonRes = 1; %0 to plot only non-responders, 1 to plot responders (if UseAreaThresh == 1)
ResBasedOn = 'UncagingEvents' ; % If UseAreaThresh == 1, should responder TS be based on response of uncaging cell or neighboring network? 'NeighboringEvents' or 'UncagingEvents'
PostStimTime = []; %% If UseAreaThresh == 1, time period post-stim for determining responders (should match time_window used in 'Fig2_dFFResponders.m')

ExcludeRampingCells = 0; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 30; %60; %number of seconds at the beginning of baseline to exclude (must be a multiple of time_window)

PerformPerm = 1; %0: don't perform any permutation test 1: perform a permutation test
perm_num = 10000; %How many permutations should be carried out for each condition?
if PerformPerm == 1
    XSecBinsPrePost_PermMean = struct(); %Keep track of the permuted means for each X-sec bin pre- and post-stim
    XSecBinsChangeAvgBaseline_PermMean = struct(); %Keep track of the permuted means for the avgerage change from baseline each X-sec post-stim bin
end

pulse_num = 10;

if strcmp(param_name, 'Number')
    param_pre = 'numevents_Pre'; %Which parameter would you like to aggregate?
    param_post = 'numevents_Post';
    per_change = 'EvtCt_PerChange'; %What are the units of this parameter? put in (), e.g. '(s)'
    param_unit = '';
    if WhichCells == 1
        YLIMS = [-1 5]; %[-10 10];  [-1 7];
    else
        if exist('data_set', 'var')
            YLIMS = [-1 2]; %[-1 4]
        else
            YLIMS = [-1 1];
        end
    end
elseif strcmp(param_name, 'Area')
    param_pre = 'avgarea_Pre'; %Which parameter would you like to aggregate?
    param_post = 'avgarea_Post';
    per_change = 'area_PerChange';
    param_unit = '(um^2)'; %What are the units of this parameter? put in (), e.g. '(s)'
    if WhichCells == 1
        YLIMS = [-50 200]; %[-10 10];
    else
        YLIMS = [-50 50];
    end
elseif strcmp(param_name, 'Duration')
    param_pre = 'avgduration_Pre'; %Which parameter would you like to aggregate?
    param_post = 'avgduration_Post';
    per_change = 'duration_PerChange';
    param_unit = '(s)'; %What are the units of this parameter? put in (), e.g. '(s)'
    if WhichCells == 1
        YLIMS = [-10 15]; %[-20 40] %[-10 10];
    else
        YLIMS = [-5 5];
    end
elseif strcmp(param_name, 'Amplitude')
    param_pre = 'avgdffMax_Pre'; %Which parameter would you like to aggregate?
    param_post = 'avgdffMax_Post';
    per_change = 'dffMax_PerChange';
    param_unit = '(dF/F max)'; %What are the units of this parameter? put in (), e.g. '(s)'
    if WhichCells == 1
        YLIMS = [-5 5]; %[-20 20];
    else
        YLIMS = [-2 2];
    end
end

condition = fieldnames(mydata);
TSERIES = struct();
for cond = 1:length(condition)
    if exist('pulse_num', 'var')
        if exist('data_sets', 'var') && contains(data_sets{1,2}, 'Cx43')
            pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num)); %For WT, CX43 and CBX dataset
        else
            pulse_cond = strcat(condition{cond}, sprintf('_%d',pulse_num)); %For laser uncaging control and receptor antagonist dataset
        end
    else
        pulse_cond = extractAfter(condition{cond}, '_');
    end
    if UseAreaThresh
        temptime = sprintf('poststim%ds', PostStimTime);
        RESPONDER_TS = dFFResponders.(temptime);
        clear temptime
        
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

for cond = 1:length(condition)
    TS_combined.(condition{cond}) = []; %rows will be tseries, columns will be 30s time bins (column 1 is the baseline 30S, column 2 is the first 30S following stimulation, etc.)
    TS_combined_change.(condition{cond}) = []; %rows will be tseries, columns will be 30s time bins
    TS_combined_PrePost.(condition{cond}) = []; %rows will be tseries, columns will be 30s time bins of all pre-stim time bins followed by all post-stim time bins
    TS_combined_PerChange.(condition{cond}) = []; %rows will be tseries, columns will be 30s time bins: percent change from baseline average for that t-series
    TS_combined_PerChange_PrePost.(condition{cond}) = []; %rows will be tseries, columns will be 30s time bins of all pre-stim time bins followed by all post-stim time bins: percent change from baseline average for that t-series
    TS_combined_ChangeNormTotal.(condition{cond}) = []; %rows will be tseries, columns will be 30s time bins for post-stim activity: (post - avg(pre)) / total activity * 100
    TS_combined_ChangeNormTotal_PrePost.(condition{cond}) = []; %rows will be tseries, columns will be 30s time bins of all pre-stim time bins followed by all post-stim time bins: (time_bin - avg(pre)) / total activity * 100
%     Responder.(condition{cond}) = []; %keep track of which tseries had at least one post-event above the baseline area threshold
end

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
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
            
            BLBinsExclude = ExcludeBaseline / time_window; %how many baseline bins should be excluded?
            
            if WhichCells == 1 %To include only the uncaging cells
                Cell_idx = find(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk == 1);
                if ~ismember(Cell_idx, ramping_cell_idx) %if the uncaging cell is not a ramping cell
                    TS_val = [nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end)), mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :)];
                    TS_combined.(condition{cond}) = [TS_combined.(condition{cond}); TS_val];
                    TS_combined_PrePost.(condition{cond}) = [TS_combined_PrePost.(condition{cond});...
                        mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end), mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :)];
                    
                    PerChangePost = ((mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :) - nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end)))./...
                        nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end))) * 100;
                    PerChangePre = ((mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end) - nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end)))./...
                        nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end))) * 100;
                    TS_combined_PerChange.(condition{cond}) = [TS_combined_PerChange.(condition{cond}); PerChangePost];
                    TS_combined_PerChange_PrePost.(condition{cond}) = [TS_combined_PerChange_PrePost.(condition{cond}); PerChangePre, PerChangePost];
                    
                    TotalPrePost = sum([mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end), mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :)]); %total # of events in the cell
                    ChangeNormTotalPost = ((mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :) - nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end)))./...
                        TotalPrePost);
                    ChangeNormTotalPre = ((mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end) - nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end)))./...
                        TotalPrePost);
                    TS_combined_ChangeNormTotal.(condition{cond}) = [TS_combined_ChangeNormTotal.(condition{cond}); ChangeNormTotalPost];
                    TS_combined_ChangeNormTotal_PrePost.(condition{cond}) = [TS_combined_ChangeNormTotal_PrePost.(condition{cond}); ChangeNormTotalPre, ChangeNormTotalPost];
                    
                    PostMinusPre = mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :) - nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end));
                    TS_combined_change.(condition{cond}) = [TS_combined_change.(condition{cond});PostMinusPre];
                    
                    if PlotIndTS
                        %plotting individual cells for each t-series
                        plot([1:size(TS_combined_PrePost.(condition{cond}),2)], [mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, :), mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :)])
                        xlim([0 (size(TS_combined_PrePost.(condition{cond}),2)+1)])
                        ylabel(sprintf('event %s', lower(param_name)))
                        title(sprintf('%s: t-series %d uncaging cells', strrep(condition{cond}, '_', ' '), tseries))
                        xtickangle(45)
                        xticks([1:size(TS_combined_PrePost.(condition{cond}),2)])
                        xticklabels({'-150-120', '-120-90', '-90-60', '-60-30', '-30-0', '0-30', '30-60','60-90','90-120','120-150'});
                        cd(save_dir)
                        saveas(gcf, sprintf('Cyto_2P%sUncagingCells_Evt%s%dSbins_TS%d.tif', condition{cond}, param_name, time_window, tseries))
                        cd(start_dir)
                        close
                    end
                else %if the uncaging cell is a ramping cell, place nans in the row for that t-series
                    NumBinsPre = size(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end),2);
                    NumBinsPost = size(mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :),2);
                    
                    TS_val = [nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, BLBinsExclude+1:end)), mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :)];
                    TS_combined.(condition{cond}) = [TS_combined.(condition{cond}); nan(1, size(TS_val,2))];
                    TS_combined_PrePost.(condition{cond}) = [TS_combined_PrePost.(condition{cond}); nan(1, NumBinsPre + NumBinsPost)];
                    
                    TS_combined_PerChange.(condition{cond}) = [TS_combined_PerChange.(condition{cond}); nan(1, NumBinsPost)];
                    TS_combined_PerChange_PrePost.(condition{cond}) = [TS_combined_PerChange_PrePost.(condition{cond}); nan(1, NumBinsPre + NumBinsPost)];
                    
                    TS_combined_ChangeNormTotal.(condition{cond}) = [TS_combined_ChangeNormTotal.(condition{cond}); nan(1, NumBinsPost)];
                    TS_combined_ChangeNormTotal_PrePost.(condition{cond}) = [TS_combined_ChangeNormTotal_PrePost.(condition{cond}); nan(1, NumBinsPre + NumBinsPost)];
                    
                    TS_combined_change.(condition{cond}) = [TS_combined_change.(condition{cond}); nan(1, NumBinsPost)];
                    clear NumBinsPre NumBinsPost
                    
                end
                clear Cell_idx TS_val PostMinusPre PerChangePost PerChangePre TotalPrePost ChangeNormTotalPost ChangeNormTotalPre ramping_cell_idx
                
            elseif WhichCells == 2 %To include neighboring cells with activity
                Cell_idx = find(isnan(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk)); %Identifying the neighboring cells
                Cell_idx = Cell_idx(~ismember(Cell_idx, ramping_cell_idx)); %exclude any ramping cells from the group (if ExcludeRampingCells = 1)
                sum_evts = sum(mydata(tseries).(condition{cond}).numevents_PreVPost,2); %Sum the events pre & post in each cell
                change_inf = mean(mydata(tseries).(condition{cond}).(TB).(per_change),2);
                Nonzero_cells = Cell_idx(sum_evts(Cell_idx) > 0); %Modified 3/21/2022 Find the indices of neighboring cells with at least one event
%                 Nonzero_cells = Cell_idx(sum_evts(Cell_idx) > 0 & isfinite(change_inf(Cell_idx))); %Find the indices of neighboring cells with at least one event and without infinity as a value
                TS_val = [nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Nonzero_cells, BLBinsExclude+1:end),'all'), nanmean(mydata(tseries).(condition{cond}).(TB).(param_post)(Nonzero_cells, :),1)];
                TS_combined.(condition{cond}) = [TS_combined.(condition{cond}); TS_val];
                TS_combined_PrePost.(condition{cond}) = [TS_combined_PrePost.(condition{cond});...
                    nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Nonzero_cells, BLBinsExclude+1:end),1), nanmean(mydata(tseries).(condition{cond}).(TB).(param_post)(Nonzero_cells, :),1)];
                
                AvgPre = nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Nonzero_cells, BLBinsExclude+1:end),2); %average # events/30S for each cell at baseline
                ChangePost = mydata(tseries).(condition{cond}).(TB).(param_post)(Nonzero_cells, :) - AvgPre; %each post-stim bin - avgPre for each cell
                ChangePre = mydata(tseries).(condition{cond}).(TB).(param_pre)(Nonzero_cells, BLBinsExclude+1:end) - AvgPre; %each pre-stim bin - avgPre for each cell
                PerChangePost = nanmean((ChangePost ./ AvgPre) * 100, 1); %average the percent change at each time point across cells/regions
                PerChangePre = nanmean((ChangePre ./ AvgPre) * 100, 1); %average the percent change at each time point across cells/regions
                TS_combined_PerChange.(condition{cond}) = [TS_combined_PerChange.(condition{cond}); PerChangePost];
                TS_combined_PerChange_PrePost.(condition{cond}) = [TS_combined_PerChange_PrePost.(condition{cond}); PerChangePre, PerChangePost];
                
                TotalPrePost = sum([mydata(tseries).(condition{cond}).(TB).(param_pre)(Nonzero_cells, BLBinsExclude+1:end), mydata(tseries).(condition{cond}).(TB).(param_post)(Nonzero_cells, :)], 2); %total # of events in the cell
                ChangeNormTotalPost = nanmean((ChangePost ./ TotalPrePost), 1); %change from avgerage baseline normalized by total activity in that region at each time point averaged across cells/regions
                ChangeNormTotalPre = nanmean((ChangePre ./ TotalPrePost), 1); %change from average baseline normalized by total activity in that region at each time point averaged across cells/regions
                TS_combined_ChangeNormTotal.(condition{cond}) = [TS_combined_ChangeNormTotal.(condition{cond}); ChangeNormTotalPost];
                TS_combined_ChangeNormTotal_PrePost.(condition{cond}) = [TS_combined_ChangeNormTotal_PrePost.(condition{cond}); ChangeNormTotalPre, ChangeNormTotalPost];
                
                PostMinusPre = nanmean(ChangePost,1);
                TS_combined_change.(condition{cond}) = [TS_combined_change.(condition{cond});PostMinusPre];
                
                if PlotIndTS
                    %plotting individual cells for each t-series
                    plot([1:size(TS_combined_PrePost.(condition{cond}),2)], [mydata(tseries).(condition{cond}).(TB).(param_pre)(Cell_idx, :), mydata(tseries).(condition{cond}).(TB).(param_post)(Cell_idx, :)])
                    xlim([0 (size(TS_combined_PrePost.(condition{cond}),2)+1)])
                    ylabel(sprintf('event %s', lower(param_name)))
                    title(sprintf('%s: t-series %d neighboring cells', strrep(condition{cond}, '_', ' '), tseries))
                    xtickangle(45)
                    xticks([1:size(TS_combined_PrePost.(condition{cond}),2)])
                    xticklabels({'-150-120', '-120-90', '-90-60', '-60-30', '-30-0', '0-30', '30-60','60-90','90-120','120-150'});
                    cd(save_dir)
                    saveas(gcf, sprintf('Cyto_2P%sNeighboringCells_Evt%s%dSbins_TS%d.tif', condition{cond}, param_name, time_window, tseries))
                    cd(start_dir)
                    close
                end
                
                clear Cell_idx TS_val sum_evts Nonzero_cells change_inf PostMinusPre PerChangePost PerChangePre...
                    AvgPre ChangePost ChangePre TotalPrePost ChangeNormTotalPost ChangeNormTotalPre ramping_cell_idx
            end
        else
            NumBinsPre = size(TS_combined_PrePost.(condition{cond}),2) - size(TS_combined_change.(condition{cond}),2);
            NumBinsPost = size(TS_combined_change.(condition{cond}),2);
            
            TS_combined.(condition{cond}) = [TS_combined.(condition{cond}); nan(1, size(TS_combined.(condition{cond}),2))];
            TS_combined_PrePost.(condition{cond}) = [TS_combined_PrePost.(condition{cond}); nan(1, NumBinsPre + NumBinsPost)];
            
            TS_combined_PerChange.(condition{cond}) = [TS_combined_PerChange.(condition{cond}); nan(1, NumBinsPost)];
            TS_combined_PerChange_PrePost.(condition{cond}) = [TS_combined_PerChange_PrePost.(condition{cond}); nan(1, NumBinsPre + NumBinsPost)];
            
            TS_combined_ChangeNormTotal.(condition{cond}) = [TS_combined_ChangeNormTotal.(condition{cond}); nan(1, NumBinsPost)];
            TS_combined_ChangeNormTotal_PrePost.(condition{cond}) = [TS_combined_ChangeNormTotal_PrePost.(condition{cond}); nan(1, NumBinsPre + NumBinsPost)];
            
            TS_combined_change.(condition{cond}) = [TS_combined_change.(condition{cond}); nan(1, NumBinsPost)];
            clear NumBinsPre NumBinsPost
        end
    end
    
    if PerformPerm == 1 %only include t-series with the chosen pulse number 
        sprintf('START: %s', datetime)
        
        NBinsPost = size(TS_combined_change.(condition{cond}),2);
        NBinsPre = size(TS_combined_PrePost.(condition{cond}), 2)-NBinsPost;
        [XSecBinsPrePost_PermMean.(condition{cond}),XSecBinsChangeAvgBaseline_PermMean.(condition{cond})] = MeanFt_XSecStatBins_PermEventOnsetByRegion(mydata, condition{cond}, TSERIES.(condition{cond}),...
            time_window, ExcludeBaseline, perm_num, param_name, WhichCells, ExcludeRampingCells, RampingCellThresh, NBinsPre, NBinsPost);
        clear NBinsPost NBinsPre
        
        sprintf('END: %s', datetime)
    end
end
%% B. Plot each condition vs permuted data individually in post-stim time bins of a chosen duration
%
% For Ext. Data Fig. 2b & Ext. Data Fig. 3g-h: ERRORBARS = 'sem';
% plotchange = [1]; AltXTicks = 0; TimeS_BeforeAfterUncaging isn't relevant
%
% For Ext. Data Fig. 7d: ERRORBARS = 'sem'; plotchange = [1]; AltXTicks =
% 1; TimeS_BeforeAfterUncaging = [-120 570]; (in this dataset, use -120s
% instead of -150s, since the first 30s timebin was cropped in
% ExtendedDataFig7_PreppingDataStruct.m'

ERRORBARS = 'sem';
plotchange = [1]; %0 = plot pre- and post-stim bins, 1 = plot post-stim bins change from average baseline
AltXTicks = 1; %0 to have x-tick labels manually defined, with p-vals from permutation test, 1 to use TimeS_BeforeAfterUncaging to define tick marks
TimeS_BeforeAfterUncaging = [-120 570];%[-150 600]; %if AltXTicks = 1: first element is time (in s) recorded before uncaging, second element is time (in s) recorded after the uncaging

if PerformPerm == 1
    for cond = 1:length(condition)
        for PT = plotchange
            if WhichCells == 1
                Cell2Plot = 'uncaging';
            elseif WhichCells == 2
                Cell2Plot = 'neighboring';
            end
            
            if PT == 0
                exp_data = TS_combined_PrePost.(condition{cond});
                exp_mean = mean(exp_data(TSERIES.(condition{cond}),:), 1, 'omitnan');
                perm_mean = XSecBinsPrePost_PermMean.(condition{cond});
                
                if AltXTicks == 0
                    XTICKLABELS = {'-150-120', '-120-90','-90-60', '-60-30', '-30-0', '0-30', '30-60','60-90','90-120','120-150', '150-180', '180-210', '210-240', '240-270', '270-300'};
                    XTICKLABELS = XTICKLABELS(BLBinsExclude+1:end);
                elseif AltXTicks == 1
                    XTICKLABELS = TimeS_BeforeAfterUncaging(1): time_window : TimeS_BeforeAfterUncaging(2);
                    XTICKLABELS(XTICKLABELS==0) = [];
                    XTICKLABELS = XTICKLABELS(BLBinsExclude+1:end);
                end
                
                YLABEL = sprintf('mean %s of events %s', lower(param_name), param_unit);
                YMIN = -1;
                pchange = '';
                if exist('pulse_num', 'var')
                    TITLE = {sprintf('%s %s %s cell(s) %d-100ms pulses', indicator, strrep(condition{cond}, '_', ' '), Cell2Plot, pulse_num),...
                        sprintf('event %s per %dS bins', lower(param_name), time_window)};
                else
                    TITLE = {sprintf('%s %s-100ms pulses %s cell(s)', indicator, strrep(condition{cond}, '_', ' '), Cell2Plot),...
                        sprintf('event %s per %dS bins', lower(param_name), time_window)};
                end
            elseif PT == 1
                exp_data = TS_combined_change.(condition{cond});
                exp_mean = mean(exp_data(TSERIES.(condition{cond}),:), 1, 'omitnan');
                perm_mean = XSecBinsChangeAvgBaseline_PermMean.(condition{cond});
                
                if AltXTicks == 0
                    XTICKLABELS = {'0-30', '30-60','60-90','90-120','120-150'}; %{'0-30', '30-60','60-90','90-120'};
                elseif AltXTicks == 1
                    XTICKLABELS = time_window : time_window: TimeS_BeforeAfterUncaging(2);
                end
                
                YLABEL = sprintf('change in mean %s of events %s', lower(param_name), param_unit);
                YMIN = -1; %-2;
                pchange = 'Change';
                if exist('pulse_num', 'var')
                    TITLE = {sprintf('%s %s %s cell(s) %d-100ms pulses', indicator, strrep(condition{cond}, '_', ' '), Cell2Plot, pulse_num),...
                        sprintf('event %s per %dS bins', lower(param_name), time_window)};
                else
                    TITLE = {sprintf('%s %s-100ms pulses %s cell(s)', indicator, strrep(condition{cond}, '_', ' '), Cell2Plot),...
                        sprintf('event %s per %dS bins', lower(param_name), time_window)};
                end
            end
            
            if PT == 0
                PermGreaterExp = perm_mean >= exp_mean; %if plotting actual value, test specifically for increase from permuted mean (1-tailed)
            elseif PT == 1
                PermGreaterExp = abs(perm_mean) >= abs(exp_mean); %if plotting change from baseline, test if changes in experimental mean are greater than expected by chance (tests for increase or decrease 2-tailed)
            end
            
            pval_PerFrame.(condition{cond}) = (sum(PermGreaterExp, 1)+1) ./ (perm_num+1);
            % correct for multiple comparisons using Benjamini Yekutieli FDR
            % method with FDR of 5%
            [h.(condition{cond}), ~, ~, adj_pval_PerFrame.(condition{cond})] = fdr_bh(pval_PerFrame.(condition{cond}), 0.05, 'dep');
            sig_frames.(condition{cond}) = find(h.(condition{cond}));
            
            if AltXTicks == 0
                XTICKLABELS_PVAL = cell(3, length(XTICKLABELS));
                for BINS = 1:length(XTICKLABELS)
                    XTICKLABELS_PVAL{1, BINS} = sprintf('%s ', XTICKLABELS{BINS});
                    XTICKLABELS_PVAL{2, BINS} = sprintf('p=%.4f ', pval_PerFrame.(condition{cond})(BINS));
                    XTICKLABELS_PVAL{3, BINS} = sprintf('adj p=%.4f ', adj_pval_PerFrame.(condition{cond})(BINS));
                    clear idx
                end
                clear BINS
                tickLabels = MultiLineLabels(XTICKLABELS_PVAL); %allows for multiline xlabels
            elseif AltXTicks == 1
                tickLabels = XTICKLABELS;
            end
            
            Structs2Plot = {perm_mean, exp_data};
            TS2Plot = [{{1:perm_num}}, {{TSERIES.(condition{cond})}}];

            MAX_VAL = max(max(nanmean(perm_mean, 1)), max(nanmean(exp_mean, 1)));
            BUFFER_OrderMag = floor(log(abs(MAX_VAL))./log(10));
            BUFFER = 10^BUFFER_OrderMag;
            Pval_placement = MAX_VAL + BUFFER;
            clear MAX_VAL BUFFER_OrderMag BUFFER
            
            %             YLIMS = [YMIN Pval_placement+2];
            %             YLIMS = [-1 2];%[-1 2]; %[-1 3]; [-1 5]; %
            
            DataLabels = {'permuted data', strrep(condition{cond}, '_', '+')};
            if ~isempty(sig_frames.(condition{cond}))
                DataLabels = [DataLabels, 'p<0.05'];
            end
            
            if strcmp(extractBefore(condition{cond}, '_'), 'RuBiGABA') %if it's a GABA condition, the color will be magenta
                DataColors = {'-ko', '-mo'};
            elseif strcmp(extractBefore(condition{cond}, '_'), 'RuBiGlu') %if it's a Glu condition, the color will be cyan
                DataColors = {'-ko', '-co'};
            else
                DataColors = {'-ko', '-go'};
            end
            
            %Plot permuted mean in black and experimental mean in magenta or
            %cyan
            PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt( Structs2Plot, TS2Plot, 1,...
                TITLE, YLABEL, YLIMS, sprintf('%dS time bins', time_window), [0 size(Structs2Plot{1,1},2)+1], [1:1:size(Structs2Plot{1,1},2)],...
                tickLabels, '', DataColors, ERRORBARS)
%             xtickangle(45)
            
            %plot the bins that are significantly above shuffled data
            SF = sig_frames.(condition{cond});
            scatter(SF, ones(1, length(sig_frames.(condition{cond})))* Pval_placement,  'k', 'fill')
            
            legend(DataLabels, 'Location', 'best')
            set(gcf,'position', [100, 200, 600, 400])
            
            cd(save_dir)
            if exist('pulse_num', 'var')
                saveas(gcf, sprintf('%s%s%sCell%dPul_%sEvt%s%dSbins_Mean%s_Perm%d.tif', indicator, condition{cond}, upper(Cell2Plot), pulse_num, pchange, param_name, time_window, upper(ERRORBARS), perm_num))
                saveas(gcf, sprintf('%s%s%sCell%dPul_%sEvt%s%dSbins_Mean%s_Perm%d.svg', indicator, condition{cond}, upper(Cell2Plot), pulse_num, pchange, param_name, time_window, upper(ERRORBARS), perm_num))
            else
                saveas(gcf, sprintf('%s%sPul%sCell_%sEvt%s%dSbins_Mean%s_Perm%d.tif', indicator, condition{cond}, upper(Cell2Plot), pchange, param_name, time_window, upper(ERRORBARS), perm_num))
                saveas(gcf, sprintf('%s%sPul%sCell_%sEvt%s%dSbins_Mean%s_Perm%d.svg', indicator, condition{cond}, upper(Cell2Plot), pchange, param_name, time_window, upper(ERRORBARS), perm_num))
            end
            close
            cd(start_dir)
            
            clear Cell2Plot DataColors DataLabels pchange YLABEL XTICKLABELS TS2Plots Structs2Plot...
                SF Pval_placement exp_mean perm_mean PermGreaterExp pval_PerFrame h adj_pval_PerFrame...
                sig_frames YMIN TITLE TS2Plot XTICKLABELS_PVAL exp_data tickLabels
        end
    end
end

clear ERRORBARS cond plotchange PT TimeS_BeforeAfterUncaging AltXTicks
%%
clear WhichCells tseries Time param_pre param_post per_change TS_combined TS_combined_change TS_combined_PrePost...
    save_dir start_dir TB time_window pulse_num pulse_idx GABA_multcomp Glu_multcomp GABA_multcomp_array Glu_multcomp_array...
    UseAreaThresh ResNonRes PostStimTime ResBasedOn...
    YLIMS param_name param_unit condition cond DataColors data_set plotchange PlotIndTS TS_combined_PerChange TS_combined_PerChange_PrePost...
    plot_absChange TS_combined_ChangeNormTotal TS_combined_ChangeNormTotal_PrePost ERRORBARS ExcludeRampingCells BLBinsExclude ExcludeBaseline...
    PerformPerm perm_num RampingCellThresh XSecBinsPrePost_PermMean XSecBinsChangeAvgBaseline_PermMean TSERIES ans  YLIMS