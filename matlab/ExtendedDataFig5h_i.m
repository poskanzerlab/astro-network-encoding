%% Extended Data Fig. 5h-i: Spearman correlation between features of baseline activity 
%  and relative post-stim propagative event frequency in neighboring astrocytes
%  Ext. Data Fig. 5h: the feature of baseline activity is baseline propagative
%                     event rate (# propagative events/30s)
%  Ext. Data Fig. 5i: the feature of baseline activity is fraction of 
%       propagative events at baseline (out of all baseline events within a
%       single cell, what fraction of them are propagative)
%
% Modified from PropEvtFreq_BLvsRelativeResponse_MC20231106.m
%
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%     Section E: time_window = 30;
%  3. Run the following code:
%
%     For Fig. 5h-i:
%       A. Aggregate together eligible cells for each condition and sort on
%          baseline propagation probability for each cell (out of all baseline events in a cell, what proportion of them are propagative): 
%          time_window = 30; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; 
%          ExcludeBaseline = 90; TotalBaseline = 150; PostStimPeriod = 120;
%          pulse_num = 10; UseAreaThresh = 0;
%
%     For Fig. 5h:
%       B. As is
%     For Fig. 5i:
%       C. As is
%     
%
% Michelle Cahill 20240115
%% A. Aggregate eligible cells (neighboring cells, non-ramping, with some propagative events in the baseline period (60-0s pre-stim))
% One struct with each eligible cell as its own row
%
% For Ext. Data Fig. 5h-i: time_window = 30; ExcludeRampingCells = 1;
%    RampingCellThresh = 'pval_1'; ExcludeBaseline = 90; TotalBaseline = 150;
%    PostStimPeriod = 120; pulse_num = 10; UseAreaThresh = 0;
%
time_window = 30; %time window (in seconds) taken before the stim and after the stim
TB = sprintf('StaticTimeBins%ds', time_window);

ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 90; %number of seconds at the beginning of baseline to exclude (must be a multiple of time_window)
TotalBaseline = 150; %number of seconds recorded for baseline
PostStimPeriod = 120; % number of seconds to include post-stim (must be a multiple of time_window)
UseAreaThresh = 0; %0= take all tseries from that pulse number, 1= take only tseries exceeding a specific threshold

pulse_num = 10;

condition = fieldnames(mydata);
condition = condition(1:2);

TSERIES = struct();
for cond = 1:length(condition)
    if exist('pulse_num', 'var')
        if contains(data_sets{1,2}, 'Cx43')
            pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num)); %For WT, CX43 and CBX dataset
        else
            pulse_cond = strcat(condition{cond}, sprintf('_%d',pulse_num)); %For laser uncaging control and receptor antagonist dataset
        end
    else
        pulse_cond = extractAfter(condition{cond}, '_');
    end
    if UseAreaThresh
%         if isempty(PostStimTime)
%             RESPONDER_TS = dFFResponders.poststim;
%         else
            temptime = sprintf('poststim%ds', PostStimTime);
            RESPONDER_TS = dFFResponders.(temptime);
            clear temptime
%         end
        
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
    PropProb_PropEvtRatePrePost.(condition{cond}) = []; %rows will be individual eligible cells, column 1: baseline propagation probability, 2: avg pre-stim (60-0s pre) propagative event rate, 3: post-stim (0-120s) propagative event rate, 4: relative post-stim response (post/pre)
end

%Aggregate eligible cells (neighboring cells, non-ramping, with some propagative events in the baseline period (60-0s pre-stim))
% One struct with each eligible cell as its own row
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
            PostStimBinsInclude = PostStimPeriod / time_window; % how many post-stim bins to include?
            
            temp_BLPropFreq = mydata(tseries).(condition{cond}).(TB).FreqGrowingProp_Pre; %rows are individual cells, columns are x second time bins pre-stim
            temp_SumBLPropEvt = sum(temp_BLPropFreq(:, BLBinsExclude+1:end), 2); %For each cell, sum the number of baseline propagative events in included baseline bins
            temp_PostStimPropFreq = mydata(tseries).(condition{cond}).(TB).FreqGrowingProp_Post; %rows are individual cells, columns are x second time bins post-stim
            
            temp_BLAllFreq = mydata(tseries).(condition{cond}).(TB).numevents_Pre;
            temp_SumBLAllEvt = sum(temp_BLAllFreq(:, BLBinsExclude+1:end), 2); %For each cell, sum the number of baseline events (all events) in included baseline bins
            
            Cell_idx = find(isnan(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk)); %Identifying the neighboring cells
            Cell_idx = Cell_idx(~ismember(Cell_idx, ramping_cell_idx)); %exclude any ramping cells from the group (if ExcludeRampingCells = 1)
            %             sum_evts = sum(mydata(tseries).(condition{cond}).numevents_PreVPost,2); %Sum the events pre & post in each cell
            Active_cells = Cell_idx(temp_SumBLPropEvt(Cell_idx) > 0); %Find the indices of neighboring cells with at least one baseline propagative event
            %             NoBLProp_cells = find(temp_SumBLPropEvt == 0);
            %             IncludeCellIdx = intersect(Active_cells,NoBLProp_cells);
            
            temp_PropagationProbability = (temp_SumBLPropEvt(Active_cells) ./ temp_SumBLAllEvt(Active_cells)); %For each cell calculate the proportion of baseline events that are propagative
            temp_PropEvtFreqPre = mean(temp_BLPropFreq(Active_cells, BLBinsExclude+1:end), 2, 'omitnan');
            temp_PropEvtFreqPost = mean(temp_PostStimPropFreq(Active_cells, 1:PostStimBinsInclude), 2, 'omitnan');
            temp_RelRespPost = temp_PropEvtFreqPost ./ temp_PropEvtFreqPre;
            
            PropProb_PropEvtRatePrePost.(condition{cond}) = [PropProb_PropEvtRatePrePost.(condition{cond});...
                temp_PropagationProbability, temp_PropEvtFreqPre, temp_PropEvtFreqPost, temp_RelRespPost]; %Include each eligible cell as it's own row
            
            clear Cell_idx sum_evts Active_cells ...
                ramping_cell_idx BLBinsExclude temp_BLPropFreq temp_SumBLPropEvt temp_PostStimPropFreq...
                PostStimBinsInclude temp_BLAllFreq temp_SumBLAllEvt temp_PropagationProbability...
                temp_PropEvtFreqPre temp_PropEvtFreqPost temp_RelRespPost
        end
    end
end
clear cond
%% B. Ext. Data Fig. 5h: Determine the correlation coefficient between baseline event frequency and relative post-stim response for neighboring cells with BL propagative activity
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig5h';

for cond = 1:length(condition)
    temp_data = PropProb_PropEvtRatePrePost.(condition{cond});
    tempData_PreFreqVPostRelResp = [temp_data(:,2), temp_data(:,4)];
    %scatter plot of pre- and relative post event rate including spearman and Pearson correlation
    [Spear_rho, Spear_corr_pval] = corr(tempData_PreFreqVPostRelResp, 'Type', 'Spearman'); %Spearman correlation
    [Pearson_rho, Pearson_corr_pval] = corr(tempData_PreFreqVPostRelResp); %Pearson correlation
    RegMod = fitlm(tempData_PreFreqVPostRelResp(:,1), tempData_PreFreqVPostRelResp(:,2));
    
    TITLE = {sprintf('%s propagative event number in neighbors with BL prop activity', strrep(condition{cond}, '_', ' ')),...
        sprintf('avg/%dS time bins %dS pre vs relative response %dS post', time_window, TotalBaseline-ExcludeBaseline, PostStimPeriod),...
        sprintf('spearman rho = %.15f, p-val = %.15f', Spear_rho(2,1), Spear_corr_pval(2,1)),...
        sprintf('Pearsons correlation coefficient = %.15f, p-val = %.15f', Pearson_rho(2,1), Pearson_corr_pval(2,1)),...
        sprintf('n = %d cells; exclude ramping cells: %d',size(temp_data, 1), ExcludeRampingCells)};
    
    XLABEL = sprintf('Baseline propagative event number/%dS', time_window);
    YLABEL = 'Propagative post-stim relative response (post/BL)';
    
    figure()
    hold on
    scatter(tempData_PreFreqVPostRelResp(:,1), tempData_PreFreqVPostRelResp(:,2), 'k', 'filled', 'MarkerFaceAlpha', 0.5)
    RM = plot(RegMod);
    delete(RM([1 3 4]))
%     hline(1.5, 'k')
    xlabel(XLABEL)
    ylabel(YLABEL)
    title(TITLE)
    ylim([0 6])
    xlim([0 7])
    axis square
    
    cd(save_dir)
    saveas(gcf, sprintf('%sAvgPropNum%dSNeighWBLProp_BLFreq%dSPre_RelResp%dSPost_CorrAndLinReg.tif', condition{cond}, time_window, TotalBaseline-ExcludeBaseline, PostStimPeriod))
    saveas(gcf, sprintf('%sAvgPropNum%dSNeighWBLProp_BLFreq%dSPre_RelResp%dSPost_CorrAndLinReg.svg', condition{cond}, time_window, TotalBaseline-ExcludeBaseline, PostStimPeriod))
    close
    cd(start_dir)
    clear TITLE RM temp_data tempData_PreFreqVPostRelResp Spear_rho Spear_corr_pval...
    Pearson_rho Pearson_corr_pval RegMod TITLE XLABEL YLABEL RM
end

clear cond start_dir save_dir
%% C. Ext. Data Fig. 5i: Determine the correlation coefficient between baseline propagation fraction and relative post-stim response for neighboring cells with BL propagative activity
%     Color-code cells by "low" and "high" baseline propagation fraction grouping 
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig5i';

for cond = 1:length(condition)
    temp_data = PropProb_PropEvtRatePrePost.(condition{cond});
    [~, sortingIDX] = sort(temp_data(:,1));
    temp_sorted_data = temp_data(sortingIDX, :);
    
    numTotalCells = size(temp_sorted_data, 1);
    Cutoff_IDX = ceil(numTotalCells/2);
    LOWER = 1:Cutoff_IDX;
    UPPER = Cutoff_IDX+1:numTotalCells;
    
    tempData_PreFreqVPostRelResp = [temp_sorted_data(:,1), temp_sorted_data(:,4)];
    %scatter plot of pre- and relative post event rate including spearman and Pearson correlation
    [Spear_rho, Spear_corr_pval] = corr(tempData_PreFreqVPostRelResp, 'Type', 'Spearman'); %Spearman correlation
    [Pearson_rho, Pearson_corr_pval] = corr(tempData_PreFreqVPostRelResp); %Pearson correlation
    RegMod = fitlm(tempData_PreFreqVPostRelResp(:,1), tempData_PreFreqVPostRelResp(:,2));
    
    TITLE = {sprintf('%s propagative event number in neighbors with BL prop activity', strrep(condition{cond}, '_', ' ')),...
        sprintf('baseline propagation probability %dS pre vs relative response %dS post', TotalBaseline-ExcludeBaseline, PostStimPeriod),...
        sprintf('spearman rho = %.15f, p-val = %.15f', Spear_rho(2,1), Spear_corr_pval(2,1)),...
        sprintf('Pearsons correlation coefficient = %.15f, p-val = %.15f', Pearson_rho(2,1), Pearson_corr_pval(2,1)),...
        sprintf('n = %d cells; exclude ramping cells: %d',size(temp_data, 1), ExcludeRampingCells)};
    
    XLABEL = 'Baseline propagative event probability';
    YLABEL = 'Propagative post-stim relative response (post/BL)';

    figure()
    hold on
    L = scatter(tempData_PreFreqVPostRelResp(LOWER,1), tempData_PreFreqVPostRelResp(LOWER,2), 'r', 'filled', 'MarkerFaceAlpha', 0.5);
    U = scatter(tempData_PreFreqVPostRelResp(UPPER,1), tempData_PreFreqVPostRelResp(UPPER,2), 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    RM = plot(RegMod);
    delete(RM([1 3 4]))
    hline(1.5, 'k')
    xlabel(XLABEL)
    ylabel(YLABEL)
    title(TITLE)
    legend([L, U], {'Low','High'})
    ylim([0 6])
    xlim([0 1])
    
    cd(save_dir)
    saveas(gcf, sprintf('%sNeighWBLProp_BLPropProb%dSPre_RelPropResp%dSPost_CorrAndLinReg.tif', condition{cond}, TotalBaseline-ExcludeBaseline, PostStimPeriod))
    saveas(gcf, sprintf('%sNeighWBLProp_BLPropProb%dSPre_RelPropResp%dSPost_CorrAndLinReg.svg', condition{cond}, TotalBaseline-ExcludeBaseline, PostStimPeriod))
    close
    cd(start_dir)
    clear TITLE RM temp_data tempData_PreFreqVPostRelResp Spear_rho Spear_corr_pval...
    Pearson_rho Pearson_corr_pval RegMod TITLE XLABEL YLABEL RM sortingIDX temp_sorted_data...
    numTotalCells Cutoff_IDX LOWER UPPER L U
end

clear cond start_dir save_dir
%%
clear PostStimPeriod ProProb_PropEvtRatePrePost TB time_window tseries TSERIES TotalBaseline...
    ExcludeBaseline ExcludeRampingCells PropProb_PropEvtRatePrePost pulse_num...
    RampingCellThresh UseAreaThresh condition