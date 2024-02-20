%% Fig. 1C bottom: plot the change in event count/region for defined time bins
% Test the difference from permuted data
%
%  1. Load 'AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat' workspace
%  2. Run'Fig1_PreppingDataStruct.m':
%     Section A: experiment_type = 'BathApp'; method_AgonistEntry = 'MaxCurve'
%     Section B: Run as is
%     Section C: time_window = 60
%  3. Run the following code
%
% For Fig. 1C bottom (change in event #/60s compared to baseline):
%  Section A.
%       time_window = 60; param_name = 'Number'; PlotIndTS = 0;
%       UseAreaThresh = 0; ExcludeRampingCells = 0; ExcludeBaseline = 0;
%       PerformPerm = 1; perm_num = 10000; pulse_num = [50];
%  Section B.
%       ERRORBARS = 'sem'; plotchange = [1]; IncludePreBins = [1];
%
% Calls the following functions:
% 'MeanFt_XSecStatBins_PermEventOnsetByRegion_NoLmk.m' for permuation testing
% 'PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt.m' for creating
% the plots
% 'fdr_bh.m' for correcting for multiple comparisons
%
% Michelle Cahill 20231213
%% Section 1: Calculate the change in the number of events/region/60s for each t-series in each condition
% Run a permuation test to to test the difference of experimental data from
% permuted/shuffled data. Note: Running 10,000 permutations for all
% treatment values (5, 25, 50 & 100 uM) takes about 2hrs, with each
% condition with 20-24 recordings taking ~15min
time_window = 60; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig1C\Bottom';

param_name = 'Number'; % What parameter to plot? 'Number', 'Area', 'Duration', 'Amplitude' or 'PropFrequency'

PlotIndTS = 0; %1 = plot raw counts for each cell by tseries, 0 = skip this step

UseAreaThresh = 0; %0= take all tseries from that pulse number, 1= take only tseries exceeding a specific threshold
TimeThreshold = 60; %If UseAreaThresh = 1, define the time period (in sec) in which the large event must occur to be considered

ResNonRes = 1; %0 to plot only non-responders, 1 to plot responders (if UseAreaThresh == 1)

ExcludeRampingCells = 0; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If ExcludeRampingCells = 1; If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 0; %number of seconds at the beginning of baseline to exclude (must be a multiple of time_window)

PerformPerm = 1; %0: don't perform any permutation test 1: perform a permutation test
perm_num = 10000; %How many permutations should be carried out for each condition?
if PerformPerm == 1
    XSecBinsPrePost_PermMean = struct(); %Keep track of the permuted means for each X-sec bin pre- and post-stim
    XSecBinsChangeAvgBaseline_PermMean = struct(); %Keep track of the permuted means for the avgerage change from baseline each X-sec post-stim bin
    XSecBinsChangeAvgBLPrePost_PermMean = struct(); %Keep track of the permuted means for the avgerage change from baseline each X-sec pre- & post-stim bin
end

pulse_num = [50]; %[5 25 50 100]; % pulse_num and TreatmentVals refer to the same thing (agonist concentration for each recording)

if strcmp(param_name, 'Number')
    param_pre = 'numevents_Pre'; %Which parameter would you like to aggregate?
    param_post = 'numevents_Post';
    per_change = 'EvtCt_PerChange'; %What are the units of this parameter? put in (), e.g. '(s)'
    param_unit = '';
    YLIMS = [-1 10]; %[-1 4]
elseif strcmp(param_name, 'PropFrequency')
    param_pre = 'FreqGrowingProp_Pre'; %Which parameter would you like to aggregate?
    param_post = 'FreqGrowingProp_Post';
    per_change = '';
    param_unit = ''; %What are the units of this parameter? put in (), e.g. '(s)'
        YLIMS = [-1 4];
elseif strcmp(param_name, 'Area')
    param_pre = 'avgarea_Pre'; %Which parameter would you like to aggregate?
    param_post = 'avgarea_Post';
    per_change = 'area_PerChange';
    param_unit = '(um^2)'; %What are the units of this parameter? put in (), e.g. '(s)'
        YLIMS = [-50 50];
elseif strcmp(param_name, 'Duration')
    param_pre = 'avgduration_Pre'; %Which parameter would you like to aggregate?
    param_post = 'avgduration_Post';
    per_change = 'duration_PerChange';
    param_unit = '(s)'; %What are the units of this parameter? put in (), e.g. '(s)'
        YLIMS = [-5 5];
elseif strcmp(param_name, 'Amplitude')
    param_pre = 'avgdffMax_Pre'; %Which parameter would you like to aggregate?
    param_post = 'avgdffMax_Post';
    per_change = 'dffMax_PerChange';
    param_unit = '(dF/F max)'; %What are the units of this parameter? put in (), e.g. '(s)'
        YLIMS = [-2 2];
end

condition = fieldnames(mydata);

for cond = 1:length(condition)
    for TV = pulse_num
        concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
        TS_combined.(concentration_TS_idx) = []; %rows will be tseries, columns will be 30s time bins (column 1 is the baseline 30S, column 2 is the first 30S following stimulation, etc.)
        TS_combined_change.(concentration_TS_idx) = []; %rows will be tseries, columns will be 30s time bins
        TS_combined_changePrePost.(concentration_TS_idx) = []; %rows will be tseries, columns will be 30s time bins showing change from average baseline: all pre-stim time bins followed by all post-stim time bins
        TS_combined_PrePost.(concentration_TS_idx) = []; %rows will be tseries, columns will be 30s time bins of all pre-stim time bins followed by all post-stim time bins
        TS_combined_PerChange.(concentration_TS_idx) = []; %rows will be tseries, columns will be 30s time bins: percent change from baseline average for that t-series
        TS_combined_PerChange_PrePost.(concentration_TS_idx) = []; %rows will be tseries, columns will be 30s time bins of all pre-stim time bins followed by all post-stim time bins: percent change from baseline average for that t-series
        TS_combined_ChangeNormTotal.(concentration_TS_idx) = []; %rows will be tseries, columns will be 30s time bins for post-stim activity: (post - avg(pre)) / total activity * 100
        TS_combined_ChangeNormTotal_PrePost.(concentration_TS_idx) = []; %rows will be tseries, columns will be 30s time bins of all pre-stim time bins followed by all post-stim time bins: (time_bin - avg(pre)) / total activity * 100
        Responder.(concentration_TS_idx) = []; %keep track of which tseries had at least one post-event above the baseline area threshold
        clear concentration_TS_idx
    end
    clear TV
end

for cond = 1:length(condition)
    for TV = pulse_num
        concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
        for tseries = pulse_x_idx.(concentration_TS_idx)
            if ~isempty(mydata(tseries).(condition{cond}))
                if UseAreaThresh
                    if AboveThresh.(condition{cond})(tseries) == 1
                        Responder.(condition{cond}) = [Responder.(condition{cond}), tseries];
                    end
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
                
                BLBinsExclude = ExcludeBaseline / time_window; %how many baseline bins should be excluded?
                
                Cell_idx = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2); %Identifying the neighboring cells
                Cell_idx = Cell_idx(~ismember(Cell_idx, ramping_cell_idx)); %exclude any ramping cells from the group (if ExcludeRampingCells = 1)
                EvtsPre = cell2mat(cellfun(@(x) ~isempty(x), mydata(tseries).(condition{cond}).(TB).event_idx_pre, 'UniformOutput', 0));
                EvtsPost = cell2mat(cellfun(@(x) ~isempty(x), mydata(tseries).(condition{cond}).(TB).event_idx_post, 'UniformOutput', 0));
                sum_evts = sum([EvtsPre,EvtsPost],2); %Sum the events pre & post in each cell
                Nonzero_cells = Cell_idx(sum_evts(Cell_idx) > 0); %Modified 3/21/2022 Find the indices of neighboring cells with at least one event and without infinity as a value
                
                if ~isempty(TS_combined.(concentration_TS_idx))
                    CurrSizePost = size(TS_combined.(concentration_TS_idx),2)-1;
                    CurrSizePre = size(TS_combined_PerChange_PrePost.(concentration_TS_idx),2) - CurrSizePost;
                    Post_Size_Diff = CurrSizePost - size(EvtsPost,2);
                    Pre_Size_Diff = CurrSizePre - size(EvtsPre,2);
                else
                    Post_Size_Diff = 0;
                    Pre_Size_Diff = 0;
                end
                
                PRE_VALS = mydata(tseries).(condition{cond}).(TB).(param_pre)(Nonzero_cells, BLBinsExclude+1:end); 
                if Pre_Size_Diff > 0
                    PRE_VALS = [NaN(size(PRE_VALS,1), Pre_Size_Diff), PRE_VALS];
                else
                    PRE_VALS = PRE_VALS(:,1-Pre_Size_Diff:end);
                end
                
                POST_VALS = mydata(tseries).(condition{cond}).(TB).(param_post)(Nonzero_cells, :);
                if Post_Size_Diff > 0
                    POST_VALS = [POST_VALS, NaN(size(POST_VALS,1), Post_Size_Diff)];
                else
                    POST_VALS = POST_VALS(:,1:end+Post_Size_Diff);
                end
                
                TS_val = [nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Nonzero_cells, BLBinsExclude+1:end),'all'), nanmean(POST_VALS,1)];

                TS_combined.(concentration_TS_idx) = [TS_combined.(concentration_TS_idx); TS_val];
                TS_combined_PrePost.(concentration_TS_idx) = [TS_combined_PrePost.(concentration_TS_idx);...
                    nanmean(PRE_VALS,1), nanmean(POST_VALS,1)];
                
                AvgPre = nanmean(mydata(tseries).(condition{cond}).(TB).(param_pre)(Nonzero_cells, BLBinsExclude+1:end),2); %average # events/30S for each cell at baseline
                ChangePost = POST_VALS - AvgPre; %each post-stim bin - avgPre for each cell
                ChangePre = PRE_VALS - AvgPre; %each pre-stim bin - avgPre for each cell
                PerChangePost = nanmean((ChangePost ./ AvgPre) * 100, 1); %average the percent change at each time point across cells/regions
                PerChangePre = nanmean((ChangePre ./ AvgPre) * 100, 1); %average the percent change at each time point across cells/regions
                TS_combined_PerChange.(concentration_TS_idx) = [TS_combined_PerChange.(concentration_TS_idx); PerChangePost];
                TS_combined_PerChange_PrePost.(concentration_TS_idx) = [TS_combined_PerChange_PrePost.(concentration_TS_idx); PerChangePre, PerChangePost];
                
                TotalPrePost = sum([PRE_VALS, POST_VALS], 2); %total # of events in the cell
                ChangeNormTotalPost = nanmean((ChangePost ./ TotalPrePost), 1); %change from avgerage baseline normalized by total activity in that region at each time point averaged across cells/regions
                ChangeNormTotalPre = nanmean((ChangePre ./ TotalPrePost), 1); %change from average baseline normalized by total activity in that region at each time point averaged across cells/regions
                TS_combined_ChangeNormTotal.(concentration_TS_idx) = [TS_combined_ChangeNormTotal.(concentration_TS_idx); ChangeNormTotalPost];
                TS_combined_ChangeNormTotal_PrePost.(concentration_TS_idx) = [TS_combined_ChangeNormTotal_PrePost.(concentration_TS_idx); ChangeNormTotalPre, ChangeNormTotalPost];
                
                AvgChangePre = nanmean(ChangePre,1);
                AvgChangePost = nanmean(ChangePost,1);
                TS_combined_change.(concentration_TS_idx) = [TS_combined_change.(concentration_TS_idx);AvgChangePost];
                TS_combined_changePrePost.(concentration_TS_idx) = [TS_combined_changePrePost.(concentration_TS_idx); AvgChangePre, AvgChangePost];
                
                if PlotIndTS
                    %plotting individual cells for each t-series
                    plot([1:size(TS_combined_PrePost.(concentration_TS_idx),2)], [PRE_VALS, POST_VALS])
                    xlim([0 (size(TS_combined_PrePost.(concentration_TS_idx),2)+1)])
                    ylabel(sprintf('event %s', lower(param_name)))
                    title(sprintf('%s (%duM): t-series %d neighboring cells', strrep(condition{cond}, '_', ' '), TV, tseries))
                    xtickangle(45)
                    xticks([1:size(TS_combined_PrePost.(concentration_TS_idx),2)])
                    %                     xticklabels({'-150-120', '-120-90', '-90-60', '-60-30', '-30-0', '0-30', '30-60','60-90','90-120','120-150'});
                    cd(save_dir)
                    saveas(gcf, sprintf('%s_%s_%sAllCells_Evt%s%dSbins_TS%d.tif', indicator, experiment_type, concentration_TS_idx, param_name, time_window, tseries))
                    cd(start_dir)
                    close
                end
                
                clear Cell_idx TS_val sum_evts Nonzero_cells change_inf PostMinusPre PerChangePost PerChangePre...
                    AvgPre ChangePost ChangePre TotalPrePost ChangeNormTotalPost ChangeNormTotalPre ramping_cell_idx...
                    EvtsPre EvtsPost CurrSizePost CurrSizePre Post_Size_Diff Pre_Size_Diff PRE_VALS POST_VALS...
                    AvgChangePost AvgChangePre
            end
        end
        
        if PerformPerm == 1 %only include t-series with the chosen pulse number
            sprintf('START: %s', datetime)
            
            NBinsPost = size(TS_combined_change.(concentration_TS_idx),2);
            NBinsPre = size(TS_combined_PrePost.(concentration_TS_idx), 2)-NBinsPost;
            [XSecBinsPrePost_PermMean.(concentration_TS_idx),XSecBinsChangeAvgBaseline_PermMean.(concentration_TS_idx), XSecBinsChangeAvgBLPrePost_PermMean.(concentration_TS_idx)] = MeanFt_XSecStatBins_PermEventOnsetByRegion_NoLmk(mydata,...
                condition{cond}, pulse_x_idx.(concentration_TS_idx),time_window, ExcludeBaseline, perm_num, param_name, ExcludeRampingCells, RampingCellThresh, NBinsPre, NBinsPost);
            clear NBinsPost NBinsPre
            
            sprintf('END: %s', datetime)
        end
        clear tseries concentration_TS_idx
    end
    clear TV
end
clear cond
%% Section 2: Plot each condition individually mean+-sem experimental & permutation (plus individual data points from experimental data): 
% 1. average value per xsec bin pre and post-stim (plotchange = 0,
% IncludePreBins = 1)
% 2. change from average baseline post-stim (plotchange = 1, IncludePreBins = 0)
% 3. change from average baseline pre- and post-stim (plotchange = 1, IncludePreBins = 1)
ERRORBARS = 'sem';
plotchange = [1] ;%[0, 1, 1]; %0 = plot value of bins, 1 = plot change from average baseline
IncludePreBins = [1]; %[1, 1, 0] %1: plot pre- and post-stim bins, 0:plot only post-stim bins (length of plotchange and IncludePreBins  must match)

if PerformPerm == 1
    for cond = 1:length(condition)
        for TV = pulse_num
            concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
            for PLOT_TYPE = 1:length(plotchange)
                if plotchange(PLOT_TYPE) == 0 && IncludePreBins(PLOT_TYPE) == 1
                    exp_data = TS_combined_PrePost.(concentration_TS_idx);
                    exp_mean = mean(exp_data(:,:), 1, 'omitnan');
                    perm_mean = XSecBinsPrePost_PermMean.(concentration_TS_idx);
                    
                    NumBinsPost = size(TS_combined_change.(concentration_TS_idx),2);
                    NumBinsPre = size(exp_data,2)-NumBinsPost;
                    
                    XTICKLABELS = {};
                    START_TIME = -(time_window*NumBinsPre);
                    for BIN = 1:NumBinsPre+NumBinsPost
                        temp_bin_label = sprintf('%d-%d', START_TIME, START_TIME+time_window);
                        XTICKLABELS = [XTICKLABELS, temp_bin_label];
                        START_TIME = START_TIME + time_window;
                        clear temp_bin_label
                    end
                    clear BIN START_TIME NumBinsPre NumBinsPost
                    
                    YLABEL = sprintf('mean %s of events %s', lower(param_name), param_unit);
                    YMIN = -1;
                    pchange = '';
                    TITLE = {sprintf('%s %s (%duM)', indicator, strrep(condition{cond}, '_', ' '), TV),...
                        sprintf('event %s per %dS bins', lower(param_name), time_window)};
                elseif plotchange(PLOT_TYPE) == 1 && IncludePreBins(PLOT_TYPE) == 1
                    exp_data = TS_combined_changePrePost.(concentration_TS_idx);
                    exp_mean = mean(exp_data(:,:), 1, 'omitnan');
                    perm_mean = XSecBinsChangeAvgBLPrePost_PermMean.(concentration_TS_idx);
                    
                    NumBinsPost = size(TS_combined_change.(concentration_TS_idx),2);
                    NumBinsPre = size(exp_data,2)-NumBinsPost;
                    
                    XTICKLABELS = {};
                    START_TIME = -(time_window*NumBinsPre);
                    for BIN = 1:NumBinsPre+NumBinsPost
                        temp_bin_label = sprintf('%d-%d', START_TIME, START_TIME+time_window);
                        XTICKLABELS = [XTICKLABELS, temp_bin_label];
                        START_TIME = START_TIME + time_window;
                        clear temp_bin_label
                    end
                    clear BIN START_TIME NumBinsPre NumBinsPost
                    
                    YLABEL = sprintf('change in mean %s of events %s', lower(param_name), param_unit);
                    YMIN = -1;
                    pchange = 'ChangePrePost';
                    TITLE = {sprintf('%s %s (%duM)', indicator, strrep(condition{cond}, '_', ' '), TV),...
                        sprintf('event %s per %dS bins', lower(param_name), time_window)};
                elseif plotchange(PLOT_TYPE) == 1 && IncludePreBins(PLOT_TYPE) == 0
                    exp_data = TS_combined_change.(concentration_TS_idx);
                    exp_mean = mean(exp_data(:,:), 1, 'omitnan');
                    perm_mean = XSecBinsChangeAvgBaseline_PermMean.(concentration_TS_idx);
                    
                    XTICKLABELS = {};
                    START_TIME = 0;
                    for BIN = 1:size(exp_data,2)
                        temp_bin_label = sprintf('%d-%d', START_TIME, START_TIME+time_window);
                        XTICKLABELS = [XTICKLABELS, temp_bin_label];
                        START_TIME = START_TIME + time_window;
                        clear temp_bin_label
                    end
                    clear BIN START_TIME
                    YLABEL = sprintf('change in mean %s of events %s', lower(param_name), param_unit);
                    YMIN = -2;
                    pchange = 'Change';
                    TITLE = {sprintf('%s %s (%duM)', indicator, strrep(condition{cond}, '_', ' '), TV),...
                        sprintf('event %s per %dS bins', lower(param_name), time_window)};
                end
                
                if plotchange(PLOT_TYPE) == 0 %if plotting actual value, test specifically for increase from permuted mean (1-tailed)
                    PermGreaterExp = perm_mean >= exp_mean;
                elseif plotchange(PLOT_TYPE) == 1 %if plotting change from baseline, test if changes in experimental mean are greater than expected by chance (tests for increase or decrease 2-tailed)
                    PermGreaterExp = abs(perm_mean) >= abs(exp_mean);
                end
                
                pval_PerFrame.(concentration_TS_idx) = (sum(PermGreaterExp, 1)+1) ./ (perm_num+1);
                % correct for multiple comparisons using Benjamini Yekutieli FDR
                % method with FDR of 5%
                [h.(concentration_TS_idx), ~, ~, adj_pval_PerFrame.(concentration_TS_idx)] = fdr_bh(pval_PerFrame.(concentration_TS_idx), 0.05, 'dep');
                sig_frames.(concentration_TS_idx) = find(h.(concentration_TS_idx));
                
                XTICKLABELS_PVAL = cell(3, length(XTICKLABELS));
                for BINS = 1:length(XTICKLABELS)
                    XTICKLABELS_PVAL{1, BINS} = sprintf('%s ', XTICKLABELS{BINS});
                    XTICKLABELS_PVAL{2, BINS} = sprintf('p=%.4f ', pval_PerFrame.(concentration_TS_idx)(BINS));
                    XTICKLABELS_PVAL{3, BINS} = sprintf('adj p=%.4f ', adj_pval_PerFrame.(concentration_TS_idx)(BINS));
                    clear idx
                end
                clear BINS
                tickLabels = MultiLineLabels(XTICKLABELS_PVAL); %allows for multiline xlabels
                
                Structs2Plot = {perm_mean, exp_data};
                TS2Plot = [{{1:perm_num}}, {{1:size(exp_data,1)}}];
                
                MAX_VAL = max(max(exp_data));
                BUFFER_OrderMag = floor(log(abs(MAX_VAL))./log(10));
                BUFFER = 10^BUFFER_OrderMag;
                Pval_placement = MAX_VAL + BUFFER;
                
                MIN_VAL = min(min(exp_data));
                YMIN = MIN_VAL - BUFFER;
                clear MAX_VAL BUFFER_OrderMag BUFFER MIN_VAL
                
                %                 YLIMS = [YMIN Pval_placement+1]; %[YMIN Pval_placement+20]; [YMIN Pval_placement+2];
                YLIMS = [-1 7];

                DataLabels = {'permuted data', strrep(condition{cond}, '_', '+')};
                if ~isempty(sig_frames.(concentration_TS_idx))
                    DataLabels = [DataLabels, 'p<0.05'];
                end
                
                %                 if strcmp(extractBefore(condition{cond}, '_'), 'RuBiGABA') %if it's a GABA condition, the color will be magenta
                %                     DataColors = {'-ko', '-mo'};
                %                 elseif strcmp(extractBefore(condition{cond}, '_'), 'RuBiGlu') %if it's a Glu condition, the color will be cyan
                %                     DataColors = {'-ko', '-co'};
                %                 end
                if contains(condition{cond}, 'Baclofen') %if it's a GABA condition, the color will be magenta
                    DataColors = {'-ko', '-mo'};
                    temp_color = [1 0 1];
                elseif contains(condition{cond}, 'tACPD') %if it's a Glu condition, the color will be cyan
                    DataColors = {'-ko', '-co'};
                    temp_color = [0 1 1];
                end
                
                %Plot permuted mean in black and experimental mean in magenta or
                %cyan
                PlotMeanTSRepeatedMeasures_PercentInc_FlexibleCondInt( Structs2Plot, TS2Plot, TV,...
                    TITLE, YLABEL, YLIMS, sprintf('%dS time bins', time_window), [0 size(Structs2Plot{1,1},2)+1], [1:1:size(Structs2Plot{1,1},2)],...
                    tickLabels, '', DataColors, ERRORBARS)
                set(gcf,'position', [100, 200, 1200, 400])

                %Plot individual data points for experimental data
                hold on
                temp_xvals = 1:size(exp_data, 2);
                XVALS = repmat(temp_xvals, size(exp_data, 1),1);
                XVALS = reshape(XVALS, 1, []);
                YVALS = reshape(exp_data, 1, []);
                scatter(XVALS, YVALS, [],...
                    temp_color, 'fill', 'MarkerFaceAlpha', .4)
                clear temp_xvals XVALS YVALS temp_color
                
                %plot the bins that are significantly above shuffled data
                SF = sig_frames.(concentration_TS_idx);
                scatter(SF, ones(1, length(sig_frames.(concentration_TS_idx)))* Pval_placement,  'k', 'fill')
                
                legend(DataLabels, 'Location', 'best')
                
                cd(save_dir)
                saveas(gcf, sprintf('%s%s%duM_%sEvt%s%dSbins_Mean%s_Perm%d.tif', indicator, condition{cond}, TV, pchange, param_name, time_window, upper(ERRORBARS), perm_num))
                saveas(gcf, sprintf('%s%s%duM_%sEvt%s%dSbins_Mean%s_Perm%d.svg', indicator, condition{cond}, TV, pchange, param_name, time_window, upper(ERRORBARS), perm_num))
                close
                cd(start_dir)
                
                clear DataColors DataLabels pchange YLABEL XTICKLABELS TS2Plots Structs2Plot...
                    SF Pval_placement exp_mean perm_mean PermGreaterExp pval_PerFrame h adj_pval_PerFrame...
                    sig_frames YLIMS YMIN TITLE TS2Plot XTICKLABELS_PVAL exp_data
            end
            clear concentration_TS_idx PLOT_TYPE
        end
    end
end
clear ERRORBARS cond plotchange IncludePreBins TV
%%
clear BLBinsExclude condition ExcludeBaseline ExcludeRampingCells param_name...
    param_post param_pre param_unit per_change PerformPerm perm_num PlotIndTS...
    RampingCellThresh ResNonRes Responder save_dir start_dir TB tickLabels time_window...
    TimeThreshold pulse_num TS_combined TS_combined_change TS_combined_ChangeNormTotal...
    TS_combined_ChangeNormTotal_PrePost TS_combined_changePrePost TS_combined_PerChange...
    TS_combined_PerChange_PrePost TS_combined_PrePost UseAreaThresh XSecBinsChangeAvgBaseline_PermMean...
    XSecBinsChangeAvgBLPrePost_PermMean XSecBinsPrePost_PermMean ans