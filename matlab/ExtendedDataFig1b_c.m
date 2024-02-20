%% Event frequency per astrocyte (event number/60s) pre- and post-stim
% Extended Data Fig. 1b-c event frequency scatter plots
% Extended Data Fig. 1b: comparing baseline event frequency per astrocyte
% pre-Baclofen and pre-tACPD
% Extended Data Fig. 1b: baseline event frequency vs post-stim event
% frequency per astrocyte 60-120s following agonist entry into bath
% Modified from
% 'XSecStatBinsPreAvgBLvsChosenPostStim_MultCondAllReg_MC20230814.m'

%  1. Load 'AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat' workspace
%  2. Run'Fig1_PreppingDataStruct.m':
%     Section A: experiment_type = 'BathApp'; method_AgonistEntry = 'MaxCurve'
%     Section B: as is
%     Section C: time_window = 60;
%  3. Run the following code:

% A. Extended Data Fig. 1c: scatter plot of event frequency pre- and
%    post-agonist entry per region. Calls the function 'Scatter_Mean.m'
%    time_window = 60; param_name = 'Number'; PostStimBin = 2;
% B. Extended Data Fig. 1b: scatter plot comparing baseline event frequency
%    in the same astrocytes prior to agonist entry into bath 

% Michelle Cahill 20240109
%% A. Extended Data Fig. 1c: scatter plot of event frequency pre- and post-agonist entry per region 
% For each condition, record parameter values pre- and post-stim (values per active astrocyte and averaged by t-series)
% pre-stim value: averaged across all baseline time bins
% post-stim value: values from the post-stim time bin of choice (defined by
% PostStimBin)

% 1. Aggregate together values for each condition (per region and per
%    t-series)
% 2. For each condition: plot average baseline vs post-stim bin as a
%    scatter plot (calls 'Scatter_Mean.m')

time_window = 60; %time window (in seconds) for static time bins
TB = sprintf('StaticTimeBins%ds', time_window);

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1b_c';

param_name = 'Number'; % What parameter to plot? 'Number', 'Area', 'Duration', 'Amplitude' or 'PropFrequency'
PostStimBin = 2; %Which post-stim bin to plot? With time_window = 60, PostStimBin = 2 plots event frequency 60-120s following agonist entry into bath
ERRORBARS = 'sem'; %'sem' or 'std'

UseAreaThresh = 0; %0= take all tseries from that pulse number, 1= take only tseries with uncaging cells showing at least one large post-stim event

ExcludeRampingCells = 0; %0 to plot all cells of that type, 1 to exclude ramping cells
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 0; %number of seconds at the beginning of baseline to exclude (must be a multiple of time_window)

pulse_num = [50]; %[5 25 50 100]; % pulse_num and TreatmentVals refer to the same thing (agonist concentration for each recording)

if strcmp(param_name, 'Number')
    param_pre = 'numevents_Pre'; %Which parameter would you like to aggregate?
    param_post = 'numevents_Post';
    param_unit = '';
    YLIMS = [-1 10]; %[-1 4]
end

condition = fieldnames(mydata);

FTs = struct();

% 1. Aggregate together values for each condition (per region and per
% t-series)
for cond = 1:length(condition)
    for TV = pulse_num
        concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
        Agg_ByRegion = [];
        Agg_ByTseries = [];
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
                
                Cell_idx_full = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx,2); %Identifying the neighboring cells
                Cell_idx = Cell_idx_full(~ismember(Cell_idx_full, ramping_cell_idx)); %exclude any ramping cells from the group (if ExcludeRampingCells = 1)
                EvtsPre = cell2mat(cellfun(@(x) ~isempty(x), mydata(tseries).(condition{cond}).(TB).event_idx_pre, 'UniformOutput', 0));
                EvtsPost = cell2mat(cellfun(@(x) ~isempty(x), mydata(tseries).(condition{cond}).(TB).event_idx_post, 'UniformOutput', 0));
                sum_evts = sum([EvtsPre,EvtsPost],2); %Sum the events pre & post in each cell
                Nonzero_cells = Cell_idx(sum_evts(Cell_idx) > 0); %Modified 3/21/2022 Find the indices of neighboring cells with at least one event and without infinity as a value
                Zero_cells = setdiff(Cell_idx_full, Nonzero_cells);
                
                PRE_VALS = mydata(tseries).(condition{cond}).(TB).(param_pre)(:, BLBinsExclude+1:end); 
                POST_VALS = mydata(tseries).(condition{cond}).(TB).(param_post)(:, :);
                
                %Fill in regions with no events or that are ramping with
                %NaN so that region indices are maintained
                PRE_VALS(Zero_cells,:) = NaN;
                POST_VALS(Zero_cells,:) = NaN;
                
                % 1. Aggregate avg baseline per region and per t-series
                temp_AvgPreByReg = mean(PRE_VALS, 2);
                temp_AvgPreByTS = mean(temp_AvgPreByReg,'omitnan');
                
                % 2. Aggregate avg from chosen post-stim bin per region and per t-series
                temp_ChosenPostByReg = POST_VALS(:, PostStimBin);
                temp_AvgChosenPostByTS = mean(temp_ChosenPostByReg,'omitnan');
                
                temp_PrePostByReg = [temp_AvgPreByReg, temp_ChosenPostByReg];
                temp_PrePostByTS = [temp_AvgPreByTS, temp_AvgChosenPostByTS];
                
                Agg_ByRegion = [Agg_ByRegion; temp_PrePostByReg];
                Agg_ByTseries = [Agg_ByTseries; temp_PrePostByTS];
                
                clear Cell_idx Cell_idx_full EvtsPre EvtsPost sum_evts Nonzero_cells Zero_cells...
                    PRE_VALS POST_VALS temp_AvgPreByReg temp_AvgPreByTS temp_ChosenPostByReg temp_AvgChosenPostByTS...
                    temp_PrePostByReg temp_PrePostByTS ramping_cell_idx BLBinsExclude
            end
        end
        FTs.ByRegion.(param_name).(concentration_TS_idx) = Agg_ByRegion;
        FTs.ByTS.(param_name).(concentration_TS_idx) = Agg_ByTseries;
        clear Agg_ByRegion Agg_ByTseries
    end
    clear TV
end
clear cond

% 2. For each condition: plot average baseline vs post-stim period (scatter plot) 
for cond = 1:length(condition)
    for TV = pulse_num
        concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
        PLOT_BY = fieldnames(FTs);
        for p =1:length(PLOT_BY)
            struct2plot = FTs.(PLOT_BY{p}).(param_name).(concentration_TS_idx);
            TS = 1:size(struct2plot,1);
            
            YLABEL = sprintf('event %s /%dS', lower(param_name), time_window);
            TITLE = {sprintf('%s %s (%duM)', indicator, strrep(condition{cond}, '_', ' '), TV),...
            sprintf('event %s per %dS bins: %s', lower(param_name), time_window, PLOT_BY{p}),...
            sprintf('avg baseline v post-stim bin %d', PostStimBin)};
              
            % plot pre and post-stim values as scatter plots
            XLABEL = (sprintf('pre-stim %s', param_name));
            YLABEL = (sprintf('post-stim %s', param_name));
            if contains(PLOT_BY{p}, 'Region')
                AX_LIMS = [0 20];
            elseif contains(PLOT_BY{p}, 'TS')
                AX_LIMS = [0 6];
            end
            Scatter_Mean(struct2plot, TS, TITLE, XLABEL, YLABEL, ERRORBARS)
            xlim(AX_LIMS)
            ylim(AX_LIMS)
            
            cd(save_dir)
            saveas(gcf, sprintf('%s_%s%duM_Evt%s%dSBins_ScatterAvgPreVsPostStimBin%d_%s.tif', indicator, condition{cond}, TV, lower(param_name), time_window, PostStimBin, PLOT_BY{p}))
            saveas(gcf, sprintf('%s_%s%duM_Evt%s%dSBins_ScatterAvgPreVsPostStimBin%d_%s.svg', indicator, condition{cond}, TV, lower(param_name), time_window, PostStimBin, PLOT_BY{p}))
            close
            cd(start_dir)
            
            clear TITLE YLABEL struct2plot TS AX_LIMS
        end
        clear p PLOT_BY concentration_TS_idx
    end
    clear TV
end

clear cond
%% B. Extended Data Fig. 1b: compare baseline event frequency in the same astrocytes prior to agonist entry into bath (currently designed for 2 conditions)
% For each treatment value/level, combine pre-stim data for regions in all eligible recordings and plot
%    a. baseline values for each condition (with lines connecting paired
%    measurements). *This is where the paired t-test is carried out*
%    b. baseline values for each condition as a scatter plot
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1b_c';

condition = {'Baclofen', 'tACPD'}; %which two conditions to combine/compare? 
pulse_num = [50]; %[5 25 50 100]; % pulse_num and TreatmentVals refer to the same thing (agonist concentration for each recording)
EXCLUDE_NotMatching.ByRegion = {[16:22]}; %{[];[];[16:22];[249:257]}; %each pulse_num/TV is a cell, listing the indices of regions from the combined TS to exclude 
EXCLUDE_NotMatching.ByTS = {[]};%{[];[];[];[]};

CombineAcrossPulses.ByRegion = []; %Used for section b
CombineAcrossPulses.ByTS = [];
for TV = pulse_num
    PLOT_BY = fieldnames(FTs);
        for p =1:length(PLOT_BY)
            temp_comb = [];
            for cond = 1:length(condition)
                concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
                temp_comb = [temp_comb, FTs.(PLOT_BY{p}).(param_name).(concentration_TS_idx)(:,1)];
                clear concentration_TS_idx
            end
            temp_pulseIDX = find(pulse_num == TV);
            temp_exIDX = EXCLUDE_NotMatching.(PLOT_BY{p}){temp_pulseIDX, 1};
            temp_comb(temp_exIDX, :) = NaN; %For any regions or tseries that don't match between conditions, fill values with NaN
            clear temp_pulseIDX temp_exIDX 
            
            CombineAcrossPulses.(PLOT_BY{p}) = [CombineAcrossPulses.(PLOT_BY{p}); temp_comb];
            
            struct2plot = temp_comb;
            TS = 1:size(temp_comb, 1);
            
            YLABEL = sprintf('event %s /%dS pre-stim', lower(param_name), time_window);
            TITLE = {sprintf('%s (%duM)', indicator, TV),...
            sprintf('event %s per %dS bins: %s', lower(param_name), time_window, PLOT_BY{p}),...
            'avg baseline'};
        
            % a. plotting each condition's mean pre-stim value 
            PlotTSMean_PreVPost(struct2plot, TS, TITLE, YLABEL, [], ERRORBARS)
            DataLabels = cellfun(@(x) sprintf('%s Pre', x), condition, 'UniformOutput', 0);
            xticks([1 2])
            xticklabels(DataLabels)
            
            cd(save_dir)
            saveas(gcf, sprintf('%s_%duM_Evt%s%dSBins_CompAvgPre_%s.tif', indicator, TV, lower(param_name), time_window, PLOT_BY{p}))
            saveas(gcf, sprintf('%s_%duM_Evt%s%dSBins_CompAvgPre_%s.svg', indicator, TV, lower(param_name), time_window, PLOT_BY{p}))
            close
            cd(start_dir)
            
%             %b. plot pre and post-stim values as violin plots
%             figure()
%             violinplot(struct2plot)
%             title(TITLE)
%             xticklabels(DataLabels)
%             xtickangle(45)
%             ylabel(YLABEL)
%             
%             cd(save_dir)
%             saveas(gcf, sprintf('%s_%duM_Evt%s%dSBins_VioCompAvgPre_%s.tif', indicator, TV, lower(param_name), time_window, PLOT_BY{p}))
%             saveas(gcf, sprintf('%s_%duM_Evt%s%dSBins_VioCompAvgPre_%s.svg', indicator, TV, lower(param_name), time_window, PLOT_BY{p}))
%             close
%             cd(start_dir)
            
            %b. plot pre and post-stim values as scatter plots
            XLABEL = (sprintf('%s %s', DataLabels{1}, param_name));
            YLABEL = (sprintf('%s %s', DataLabels{2}, param_name));
            Scatter_Mean(struct2plot, TS, TITLE, XLABEL, YLABEL, ERRORBARS)
            xlim([0 3.5])
            ylim([0 3.5])
            
            cd(save_dir)
            saveas(gcf, sprintf('%s_%duM_Evt%s%dSBins_ScatterCompAvgPre_%s.tif', indicator, TV, lower(param_name), time_window, PLOT_BY{p}))
            saveas(gcf, sprintf('%s_%duM_Evt%s%dSBins_ScatterCompAvgPre_%s.svg', indicator, TV, lower(param_name), time_window, PLOT_BY{p}))
            close
            cd(start_dir)
            
            clear TITLE YLABEL  XLABEL struct2plot TS temp_comb cond
        end
        clear p PLOT_BY
end

clear TV
%%
clear CombineAcrossPulses FTs DataLabels DEFINED_YLIMS ERRORBARS EXCLUDE_NotMatching...
    ExcludeBaseline ExcludeRampingCells param_name param_pre param_post param_unit PostStimBin...
    pulse_num save_dir start_dir TB temp_conditions time_window tseries UseAreaThresh YLIMS RampingCellThresh