%% Fig. 2i, 2j, Extended Data Fig. 2a,c, Extended Data Fig. 3d: 
%   mean dF/F values per-and post-stim from identified AQuA events from 
%  chosen cells ('UncagingEvents' for events from the directly stimulated cell,
%  'NeighboringEvents' for events from cells not directly stimulated by uncaging)
% For a flexible number of conditions and the ability to choose which 
% cells to include for mean dF/F calculation for each recording.
% modified from 'AvgdFF_PrePost_ByTS_BySlice_MC20220531.m' (Section A & D)
%
% For Fig. 2i, 2j, Ext. Data Fig. 2c & Ext. Data Fig. 3d:
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%  3. Run section A of 'Fig2h_3c_f.m':
%     ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; ExcludeBaseline = 60; 
%     WE = 'UncagingEvents';
%  4. Run 'Fig2_dFFResponders.m': **only necessary for Fig. 2i**
%     WE = 'UncagingEvents'; ExcludeBaseline = 60; time_window = [];
%  5. Run the following sections of code:
%
%     For Fig. 2i: 
%     A. Calculate the average dF/F value pre- and post-stim from identified 
%        AQuA events from chosen regions with parameters as specificed for
%        Fig. 2i in comments of section
%     B. Plot average pre and post dF/F values for each t-series of each condition
%        for the chosen group of cells with parameters as specificed for
%        Fig. 2i in comments of section
%
%     For Fig. 2j:
%     A. Calculate the average dF/F value pre- and post-stim from identified 
%        AQuA events from chosen regions with parameters as specificed for
%        Fig. 2j in comments of section
%     C. Average change in dF/F values (post - pre) for chosen cells in
%        each recording for chosen conditions with parameters as specificed for
%        Fig. 2j in comments of section
%
%     For Extended Data Fig. 2c:
%     A. Calculate the average dF/F value pre- and post-stim from identified 
%        AQuA events from chosen regions with parameters as specificed for
%        Ext. Data Fig. 2c in comments of section
%     D. Scatter plots comparing mean baseline dF/F in uncaging cells 
%        prior to GABA and glutamate uncaging with parameters as specificed for
%        Ext. Data Fig. 2c in comments of section
%
%     For Extended Data Fig. 3d:
%     A. Calculate the average dF/F value pre- and post-stim from identified 
%        AQuA events from chosen regions with parameters as specified for
%        Ext. Data Fig. 3d in comments of section
%     C. Average change in dF/F values (post - pre) for chosen cells in
%        each recording for chosen conditions with parameters as specificed for
%        Ext. Data Fig. 3d in comments of section
%
% For Ext. Data Fig. 2a:
%  1. Load 'AQuA_CytoGCaMP_2PUncaging_WTRecAntaLaserCtrl_ExtDataFig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'CGP', 'LY', 'LaserUncagingControl'};
%     Section C & D as is
%  3. Run section A of 'Fig2h_3c_f.m':
%     ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; ExcludeBaseline = 60; 
%     WE = 'UncagingEvents';
%  4. Run the following sections of code:
%     A. Calculate the average dF/F value pre- and post-stim from identified 
%        AQuA events from chosen regions with parameters as specified for
%        Ext. Data Fig. 2a in comments of section
%     C. Average change in dF/F values (post - pre) for chosen cells in
%        each recording for chosen conditions with parameters as specificed for
%        Ext. Data Fig. 2a in comments of section
%
% Michelle Cahill 20240112
%% Section A. Calculate the average dF/F value pre- and post-stim from identified 
% AQuA events from chosen regions ('UncagingEvents' or
% 'NeighboringEvents') per tseries
% Combine together regions and tseries, calculating the average dF/F
% pre and post-stim per FOV for the chosen regions
% 
% For Fig. 2i (dF/F calculations for uncaging cells 'responding' to uncaging): 
%      ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; ExcludeBaseline = 60; 
%      time_window = []; UseAreaThresh = 1; ResNonRes = 1; ResBasedOn = 'UncagingEvents';
%      PostStimTime = time_window; pulse_num = 10; EVENTS = {'UncagingEvents', 'NeighboringEvents'};
%
% For Fig. 2j (mean dF/F calculations for the uncaging cell from all t-series),
% for Extended Data Fig. 2c (mean baseline dF/F in uncaging cells per t-series prior to uncaging) &
% for Extended Data Fig. 3d (mean dF/F calculations for neighboring cells per t-series): 
%      ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; ExcludeBaseline = 60; time_window = []; 
%      UseAreaThresh = 0; pulse_num = 10; EVENTS = {'UncagingEvents', 'NeighboringEvents'};
%      These aren't relevant: ResNonRes = 1; ResBasedOn = 'UncagingEvents'; PostStimTime = time_window;
%
% For Extended Data Fig. 2a (mean dF/F calculations for the uncaging cell
% from all t-series):
%      ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; ExcludeBaseline = 60; time_window = []; 
%      UseAreaThresh = 0; pulse_num = 10; EVENTS = {'UncagingEvents'};
%      These aren't relevant: ResNonRes = 1; ResBasedOn = 'UncagingEvents'; PostStimTime = time_window;

ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
ExcludeBaseline = 60; %number of seconds at the beginning of baseline to exclude, 0 to include all baseline period
EBL = sprintf('Exclude%dSBL', ExcludeBaseline);
time_window = []; %[] for all frames pre and post (except for excluded baseline frames) time in seconds to take pre- and post-stim

UseAreaThresh = 0; %1 to plot only responding or non-responding cells/t-series, 0 to plot all tseries/recordings in a condition
ResNonRes = 1;% If UseAreaThresh == 1, 1: plot responders, 0: plot non-responders
ResBasedOn = 'UncagingEvents' ; % If UseAreaThresh == 1, should responder TS be based on response of uncaging cell or neighboring network? 'NeighboringEvents' or 'UncagingEvents'
PostStimTime = time_window; %% If UseAreaThresh == 1, time period post-stim for determining responders (should match time_window used in 'Fig2_dFFResponders.m')

pulse_num = 10;

AVG_dFF_PrePost = struct();

if UseAreaThresh == 1
    if ResNonRes == 1
        WTS = sprintf('RespondersBy%s', extractBefore(ResBasedOn, 'Events'));
    else
        WTS = sprintf('NonRespondersBy%s', extractBefore(ResBasedOn, 'Events'));
    end
else
    WTS = 'AllTS';
end

condition = fieldnames(mydata);
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

% Combine together regions and tseries, calculating the average dF/F
% pre and post-stim per FOV for the chosen regions
condition = fieldnames(mydata);
EVENTS = {'UncagingEvents'}; %{'UncagingEvents', 'NeighboringEvents'};
for e = 1:length(EVENTS)
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
        AVG_dFF_PrePost.(EVENTS{e}).(condition{cond}) = NaN(max(pulse_x_idx.(pulse_cond)),2); %column 1: pre-stim avg column 2: post-stim average
%         AVG_dFF_PrePost.(EVENTS{e}).(condition{cond}) = NaN(max(TSERIES.(condition{cond})),2); %column 1: pre-stim avg column 2: post-stim average
        for tseries = TSERIES.(condition{cond})
            if ~isempty(mydata(tseries).(condition{cond}))
                SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
                ExcludeBaselineFrames = round(ExcludeBaseline / SecPerFrame);
                FrameOfIncrease = mydata(tseries).(condition{cond}).uncagingframe;
                NumFrames = round(time_window / SecPerFrame);
                StartFrame = FrameOfIncrease - 1 - NumFrames;
                tempFT = mydata(tseries).(condition{cond}).(EBL).(EVENTS{e}).dffMatFilter2;
                FT = mean(tempFT, 1, 'omitnan'); %average fluorescence per frame across AQuA events
                if ~isempty(FT)
                    if ~isempty(NumFrames)
                        MaxFrames = max([ExcludeBaselineFrames, StartFrame]);
                        AVG_PreStim = mean(FT(MaxFrames+1:FrameOfIncrease-1), 'omitnan');
                        AVG_PostStim = mean(FT(FrameOfIncrease:FrameOfIncrease+NumFrames), 'omitnan');
                        clear MaxFrames
                    else
                        AVG_PreStim = mean(FT(ExcludeBaselineFrames+1:FrameOfIncrease-1), 'omitnan');
                        AVG_PostStim = mean(FT(FrameOfIncrease:end), 'omitnan');
                    end
                    AVG_dFF_PrePost.(EVENTS{e}).(condition{cond})(tseries, 1) = AVG_PreStim;
                    AVG_dFF_PrePost.(EVENTS{e}).(condition{cond})(tseries, 2) = AVG_PostStim;
                end
                clear SecPerFrame ExcludeBaselineFrames FrameOfIncrease tempFT FT AVG_PreStim AVG_PostStim NumFrames StartFrame
            end
        end
    end
end

clear tseries cond e WE
%% Section B. Fig. 2i plot average pre and post dF/F values for each condition
%  for the chosen group of cells
%  Run section A prior to running this section
%  Calls the function 'PlotTSMean_PreVPost_ChooseStat.m'

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig2i';

STAT = 'nonparametric'; %'parametric' to run a t-test to compare pre- and post-stim or 'nonparametric' to run a Wilcoxon signed rank test
ERRORBARS = 'sem'; %'std' or 'sem'

GroupBy = {'tseries'}; %'tseries'
EVENTS = {'UncagingEvents'}; % {'UncagingEvents', 'NeighboringEvents'};
for g = 1:length(GroupBy)
    for e = 1:length(EVENTS)
        for cond = 1:length(condition)
            if strcmp(GroupBy{g}, 'tseries')
                AVG_Pre = AVG_dFF_PrePost.(EVENTS{e}).(condition{cond})(:,1);
                AVG_Post = AVG_dFF_PrePost.(EVENTS{e}).(condition{cond})(:,2);
                PLOT_TS = TSERIES.(condition{cond});
            end
            
            if exist('pulse_num', 'var')
                TITLE = {sprintf('%s %s %s %s %d-100ms pulses', indicator, experiment_type, strrep(condition{cond}, '_', ' '), EVENTS{e}, pulse_num),...
                    sprintf('Avg %s by %s', 'dF/F of Identified Events', GroupBy{g}),...
                    sprintf('%s (%dS post) Exclude Ramping Cells: %d, %s', WTS, PostStimTime, ExcludeRampingCells, EBL)};
            else
                TITLE = {sprintf('%s %s %s-100ms pulses %s', indicator, experiment_type, strrep(condition{cond}, '_', ' '), EVENTS{e}),...
                    sprintf('Avg %s by %s', 'dF/F of Identified Events', GroupBy{g}),...
                    sprintf('%s (%dS post) Exclude Ramping Cells: %d, %s', WTS, PostStimTime, ExcludeRampingCells, EBL)};
            end

            YLIMS = [0 1.5];
            YLABEL = 'mean fluorescence';
            
            PlotTSMean_PreVPost_ChooseStat([AVG_Pre,AVG_Post], PLOT_TS, TITLE, YLABEL, YLIMS, ERRORBARS, STAT)
            
            if exist('pulse_num', 'var')
                cd(save_dir)
                saveas(gcf, sprintf('%s%dPul_AvgFluor%dSPreVPostBy%s_ExRamp%dExBl%dS_%s%s_%dSpost_%s.tif', condition{cond}, pulse_num, time_window, GroupBy{g}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, STAT))
                saveas(gcf, sprintf('%s%dPul_AvgFluor%dSPreVPostBy%s_ExRamp%dExBl%dS_%s%s_%dSpost_%s.svg', condition{cond}, pulse_num, time_window, GroupBy{g}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, STAT))
                close
                cd(start_dir)
            else
                cd(save_dir)
                saveas(gcf, sprintf('%s_AvgFluor%dSPreVPostBy%s_ExRamp%dExBl%dS_%s%s_%dSpost_%s.tif', condition{cond}, time_window, GroupBy{g}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, STAT))
                saveas(gcf, sprintf('%s_AvgFluor%dSPreVPostBy%s_ExRamp%dExBl%dS_%s%s_%dSpost_%s.svg', condition{cond}, time_window, GroupBy{g}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, STAT))
                close
                cd(start_dir)
            end

            
            clear AVG_Pre AVG_Post PLOT_TS TITLE YLIMS YLABEL
        end
    end
end

clear e g GroupBy cond AVG_dFF_PrePost EBL ERRORBARS EVENTS ExcludeBaseline ExcludeRampingCells...
    PostStimTime pulse_cond pulse_num RampingCellThresh ResBasedOn ResNonRes...
    save_dir start_dir STAT time_window TSERIES UseAreaThresh WTS
%% Section C. Fig. 2j, Extended Data Fig. 2a & Extended Data Fig. 3d: 
%  average change in dF/F values (post - pre) for chosen cells in each recording for chosen conditions.
%  Fig. 2j: scatter plot of change in dF/F values comparing the uncaging cell 
%           response to GABA and glutamate, including a linear regression line and pearson correlation
%  Extended Data Fig. 2a & 3d: violin plots showing change in dF/F values,
%           including the chosen conditions on the same axes and testing
%           each condition's change from baseline
%  Run section A prior to running this section.
%  Calls the functions: fdr_bh.m, fwer_holmbonf.m, MultiLineLabels.m,
%  violinplot.m, cbrewer.m
%
% For Fig. 2j: condition2plot = {'WT'}; PAIRED = [1]; EXCLUDE = {''};
%              GroupBy = {'tseries'}; EVENTS = {'UncagingEvents'}; 
%              STAT = 'nonparametric'; PlotViolin = 0; PlotScatter = 1;
%
% For Extended Data Fig. 2a: condition2plot = {'AllConditions'}; PAIRED = [0]; EXCLUDE = {'WT'};
%              GroupBy = {'tseries'}; EVENTS = {'UncagingEvents'}; 
%              STAT = 'nonparametric'; PlotViolin = 1; PlotScatter = 0;
%              Correction for multiple comparisons using Bonferroni-Holm
%              correction returned in variable 'adj_pval_HolmBonferroni' (line 242),
%              is not listed on the violin plots
%
% For Extended Data Fig. 3d: condition2plot = {'WT'}; PAIRED = [1]; EXCLUDE = {''};
%              GroupBy = {'tseries'}; EVENTS = {'NeighboringEvents'}; 
%              STAT = 'nonparametric'; PlotViolin = 1; PlotScatter = 0;

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig2a';

% Choose which conditions to include and which graphs to return
condition2plot = {'AllConditions'}; %{'WT'}; % 'R1', 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' Which conditions should be plotted
PAIRED = [0]; %1 if the data within that set is paired or 0 if it's unpaired 
EXCLUDE = {'WT'};%{'WT'}; %Within that dataset, should any condition be left out? {''} to not exclude any datasets
COLORS = {'Blues'}; %'RdPu', 'Blues', 'YlGnBu' %'RdPu', 'Blues', 'MgCy', 'MgCy', 'MgCy'
COLOR_Scheme = {'seq'}; %'seq'
GroupBy = {'tseries'}; %'tseries', 'slice', 'mouse'
EVENTS = {'UncagingEvents'}; % {'UncagingEvents', 'NeighboringEvents'};
STAT = 'nonparametric'; %'parametric' to run a t-test to compare pre- and post-stim or 'nonparametric' to run a Wilcoxon signed rank test
PlotViolin = 1; %1 to plot violin plots showing change in dF/F from baseline; 0 to not return these plots
PlotScatter = 0; %1 to plot scatter plots showing change in dF/F from baseline comparing the chosen conditions (PAIRED == 1); 0 to not return these plots

clear condition
for cp = 1:length(condition2plot)
    if strcmp(condition2plot{cp}, 'AllConditions')
        temp_condition = fieldnames(mydata);
        if ~isempty(EXCLUDE{cp})
            idx = 1:length(temp_condition);
            exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
            idx = setdiff(idx, exclude_idx);
            condition = temp_condition(idx, 1);
        else
            condition = temp_condition;
        end
        clear temp_condition idx exclude_idx
        %         condition = fieldnames(mydata);
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
    
    for e = 1:length(EVENTS)
        % 1. Average change in dF/F post-pre
        struct2plot = struct();
        PlotTS = cell(1,length(condition));
        PVALS = NaN(1,length(condition));
        for cond = 1:length(condition)
            if strcmp(GroupBy{cp}, 'tseries')
                struct2plot.(condition{cond}) = AVG_dFF_PrePost.(EVENTS{e}).(condition{cond})(:,2)- AVG_dFF_PrePost.(EVENTS{e}).(condition{cond})(:,1);
                PlotTS{1,cond} = TSERIES.(condition{cond});
                if strcmp(STAT, 'parametric')
                    [~, PVALS(1,cond)] = ttest(struct2plot.(condition{cond})); %paired t-test
                elseif strcmp(STAT, 'nonparametric')
                    PVALS(1,cond) = signrank(struct2plot.(condition{cond})); %Wilcoxon signed rank test
                end
            end
        end
        
        YLABEL = 'mean change in dF/F (post-pre)';
        if strcmp(EVENTS{e}, 'UncagingEvents')
            %             YLIMS = [-0.5 0.5];
            %             YLIMS = [-2 2];
            %             YLIMS = [-1.5 1.5];
            YLIMS = [-4 4];
%             YLIMS = [-1 2];
        elseif strcmp(EVENTS{e},'NeighboringEvents')
            if strcmp(GroupBy{cp}, 'tseries')
                YLIMS = [-0.5 0.5];
            end
        end
        
        if PlotScatter
            if PAIRED(cp) == 1     
                %scatter plot of paired data including linear regression line and pearson correlation
                temp_PairedData = [struct2plot.(condition{1}), struct2plot.(condition{2})];
                temp_PairedData_selectTS = temp_PairedData(PlotTS{1,1}, :);
                temp_PairedData_NaNremoved = temp_PairedData_selectTS(all(~isnan(temp_PairedData_selectTS), 2), :);
                
                Scat_XLabel = {sprintf('%s', strrep(condition{1}, '_', ' ')), YLABEL};
                Scat_YLabel = {sprintf('%s', strrep(condition{2}, '_', ' ')), YLABEL};

                RegMod = fitlm(temp_PairedData_NaNremoved(:,1), temp_PairedData_NaNremoved(:,2));
                [rho, corr_pval] = corr(temp_PairedData_NaNremoved); %Pearson correlation
                TITLE = {sprintf('%s %s %s %d-100ms pulses', indicator, experiment_type, condition2plot{cp}, pulse_num),...
                    sprintf('Avg change in %s by %s with linear regression line','dF/F', GroupBy{cp}),...
                    sprintf('Pearsons correlation coefficient = %.4f, p-val = %.4f', rho(2,1), corr_pval(2,1)),...
                    sprintf('%s %s(%dS post): %dS pre and post', EVENTS{e}, WTS, PostStimTime, time_window), sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
                
                figure()
                hold on
                scatter(temp_PairedData_NaNremoved(:,1), temp_PairedData_NaNremoved(:,2), 'k', 'filled', 'MarkerFaceAlpha', 0.5)
                RM = plot(RegMod);
                delete(RM([1 3 4]))
                xlabel(Scat_XLabel)
                ylabel(Scat_YLabel)
                title(TITLE)
                axis square
                if strcmp(EVENTS{e}, 'UncagingEvents')
                    axis([-0.4 1.2 -0.4 1.2])
                elseif strcmp(EVENTS{e}, 'NeighboringEvents')
                    axis([-0.4 0.4 -0.4 0.4])
                end
                
                cd(save_dir)
                saveas(gcf, sprintf('%s%dPul_AvgChScatPearLinReg_dFF%dS_By%s_ExRamp%dExBl%dS_%s%s_%dSpost_%s.tif', condition2plot{cp}, pulse_num, time_window, GroupBy{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, STAT))
                saveas(gcf, sprintf('%s%dPul_AvgChScatPearLinReg_dFF%dS_By%s_ExRamp%dExBl%dS_%s%s_%dSpost_%s.svg', condition2plot{cp}, pulse_num, time_window, GroupBy{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, STAT))
                close
                cd(start_dir)
                clear TITLE rho corr_pval
                
                clear Scat_XLabel Scat_YLabel rho corr_pval temp_PairedData_NaNremoved temp_PairedData_selectTS temp_PairedData RM RegMod
            end
            clear TITLE DataLabels DataLabels_PVAL adj_pval
        end
        
        if PlotViolin
            % violin plots comparing change in each condition (post-pre) against no change
            [~, ~, ~, adj_pval] = fdr_bh(PVALS, 0.05, 'dep'); %Currently Benjamini-Yekutieli correction is reported as bottom line of data label
            adj_pval_HolmBonferroni = fwer_holmbonf(PVALS,0.05); %Correcting p-values for multiple comparison using Bonferroni-Holm correction (currently not shown in figure, corresponds to adjusted p-vals in Fig. 2A legend)
            DataLabels = cellfun(@(x) strrep(x, '_', ' '), condition, 'UniformOutput', 0);
            DataLabels_PVAL = cell(3, length(DataLabels));
            for cond = 1:length(DataLabels)
                DataLabels_PVAL{1, cond} = sprintf('%s ',DataLabels{cond});
                DataLabels_PVAL{2, cond} = sprintf('p=%.4f ',PVALS(cond));
                DataLabels_PVAL{3, cond} = sprintf('adj p=%.4f ',adj_pval(cond));
                clear idx
            end
            tickLabels = MultiLineLabels(DataLabels_PVAL); %allows for multiline xlabels
            clear cond
            
            if exist('pulse_num', 'var')
                TITLE = {sprintf('%s %s %s %d-100ms pulses', indicator, experiment_type, condition2plot{cp}, pulse_num),...
                    sprintf('Avg change in %s by %s','dF/F', GroupBy{cp}),...
                    sprintf('%s %s (%dS post): %dS pre and post', EVENTS{e}, WTS, PostStimTime, time_window), sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
            else
                TITLE = {sprintf('%s %s %s', indicator, experiment_type, condition2plot{cp}),...
                    sprintf('Avg change in %s by %s','dF/F', GroupBy{cp}),...
                    sprintf('%s %s (%dS post): %dS pre and post', EVENTS{e}, WTS, PostStimTime, time_window), sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
            end
            
            figure();
            v = violinplot(struct2plot);
            for cond = 1:length(condition)
                v(cond).ViolinColor = cmap(cmap_num(cond),:);
                v(cond).ScatterPlot.MarkerEdgeColor = cmap(cmap_num(cond),:);
                v(cond).ViolinAlpha = 0.25;
                v(cond).ShowMean = 1;
            end
            title(TITLE)
            xticklabels(tickLabels)
            ylabel(YLABEL)
            hline(0, '--k')
            ylim(YLIMS)
            
            cd(save_dir)
            if exist('pulse_num', 'var')
                saveas(gcf, sprintf('%s%dPul_Ex%s_AvgChVio_dFFPreVPost_By%s_ExRamp%dExBl%dS_%s%s%dSpost_%dSPrePost.tif', condition2plot{cp}, pulse_num, EXCLUDE{cp}, GroupBy{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, time_window))
                saveas(gcf, sprintf('%s%dPul_Ex%s_AvgChVio_dFFPreVPost_By%s_ExRamp%dExBl%dS_%s%s%dSpost_%dSPrePost.svg', condition2plot{cp}, pulse_num, EXCLUDE{cp}, GroupBy{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, time_window))
            else
                saveas(gcf, sprintf('%s_Ex%s_AvgChVio_dFFPreVPost_By%s_ExRamp%dExBl%dS_%s%s%dSpost_%dSPrePost.tif', condition2plot{cp}, EXCLUDE{cp}, GroupBy{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, time_window))
                saveas(gcf, sprintf('%s_Ex%s_AvgChVio_dFFPreVPost_By%s_ExRamp%dExBl%dS_%s%s%dSpost_%dSPrePost.svg', condition2plot{cp}, EXCLUDE{cp}, GroupBy{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS, PostStimTime, time_window))
            end
            close
            cd(start_dir)
        end
        clear struct2plot PlotTS cond TITLE YLABEL YLIMS DataLabels DataLabels_PVAL...
            PVALS BetweenConditions_pval adj_pval v idx_1 idx_2 tickLabels adj_pval_HolmBonferroni
    end
    clear e cmap cmap_num condition
end

clear condition2plot EXCLUDE COLORS COLOR_Scheme GroupBy cp AVG_dFF_PrePost...
    EBL ERRORBARS EVENTS ExcludeBaseline ExcludeRampingCells...
    PostStimTime pulse_cond pulse_num RampingCellThresh ResBasedOn ResNonRes...
    save_dir start_dir STAT time_window TSERIES UseAreaThresh WTS pulse_cond...
    pulse_num PlotViolin PlotScatter PAIRED
%% Section D. Extended Data Fig. 2c: scatter plots comparing mean baseline dF/F
%  in uncaging cells prior to GABA and glutamate uncaging 
%  Run section A prior to running this section.
%  Calls the functions Scatter_Median
%
%  For Extended Data Fig. 2c: condition2plot = {'WT'}; PAIRED = [1]; EXCLUDE = {''};
%              GroupBy = {'tseries'}; EVENTS = {'UncagingEvents'}; STAT = 'nonparametric';

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig2c';

% Choose which conditions to compare
condition2plot = {'WT'}; % 'R1', 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' Which conditions should be plotted
PAIRED = [1]; %1 if the data within that set is paired or 0 if it's unpaired (if 0, no plots will be returned)
EXCLUDE = {''};%{'WT'}; %Within that dataset, should any condition be left out?
GroupBy = {'tseries'}; %'tseries'
EVENTS = {'UncagingEvents'}; % {'UncagingEvents', 'NeighboringEvents'};
STAT = 'nonparametric';

clear condition
for cp = 1:length(condition2plot)
    if strcmp(condition2plot{cp}, 'AllConditions')
        temp_condition = fieldnames(mydata);
        if ~isempty(EXCLUDE{cp})
            idx = 1:length(temp_condition);
            exclude_idx = find(cell2mat(cellfun(@(x) contains(x, EXCLUDE{cp}), temp_condition, 'UniformOutput', 0)));
            idx = setdiff(idx, exclude_idx);
            condition = temp_condition(idx, 1);
        else
            condition = temp_condition;
        end  
        clear temp_condition idx exclude_idx
%         condition = fieldnames(mydata);
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

    for e = 1:length(EVENTS)
        % 1. Average baseline dF/F 
        struct2plot = struct();
        PlotTS = cell(1,length(condition));
        PVALS = NaN(1,length(condition));
        for cond = 1:length(condition)
            if strcmp(GroupBy{cp}, 'tseries')
                struct2plot.(condition{cond}) =  AVG_dFF_PrePost.(EVENTS{e}).(condition{cond})(:,1);
                PlotTS{1,cond} = TSERIES.(condition{cond});
            end
        end
        
        YLABEL = 'mean baseline dF/F';
        if strcmp(EVENTS{e}, 'UncagingEvents')
            YLIMS = [-1 2];
        elseif strcmp(EVENTS{e},'NeighboringEvents')
            if strcmp(GroupBy{cp}, 'tseries')
                YLIMS = [-0.5 0.5];
                %                 YLIMS = [-1 1];
            end
        end
        
        if PAIRED(cp) == 1
            if strcmp(STAT, 'parametric')
                [~, p_val] = ttest(struct2plot.(condition{1}), struct2plot.(condition{2})); %paired t-test
            elseif strcmp(STAT, 'nonparametric')
                p_val = signrank(struct2plot.(condition{1}), struct2plot.(condition{2})); %Wilcoxon signed rank test
            end
            
            TITLE = {sprintf('%s %s %s %d-100ms pulses', indicator, experiment_type, condition2plot{cp}, pulse_num),...
                sprintf('Mean baseline dF/F by %s: p = %.15f', GroupBy{cp}, p_val),...
                sprintf('%s %s %dS pre', EVENTS{e}, WTS, time_window), sprintf('Exclude Ramping Cells: %d, %s', ExcludeRampingCells, EBL)};
            
            DataLabels = cellfun(@(x) strrep(x, '_', ' '), condition, 'UniformOutput', 0);
            temp_PairedData = [struct2plot.(condition{1}), struct2plot.(condition{2})];
            
            temp_PairedData_selectTS = temp_PairedData(PlotTS{1,1}, :);
            temp_PairedData_NaNremoved = temp_PairedData_selectTS(all(~isnan(temp_PairedData_selectTS), 2), :);
          
            Scat_XLabel = {sprintf('%s', strrep(condition{1}, '_', ' ')), YLABEL};
            Scat_YLabel = {sprintf('%s', strrep(condition{2}, '_', ' ')), YLABEL};

            %scatter plot of paired data showing median and 25th and 75th percentile)
            Scatter_Median([struct2plot.(condition{1}), struct2plot.(condition{2})], PlotTS{1,1}, TITLE, Scat_XLabel, Scat_YLabel)
            if strcmp(EVENTS{e}, 'UncagingEvents')
                axis([0 1 0 1])
            elseif strcmp(EVENTS{e}, 'NeighboringEvents')
                axis([0 1 0 1])
            end
            
            cd(save_dir)
            saveas(gcf, sprintf('%s%dPul_AvgBLScatMedIQR_dFF%dSPre_By%s_ExRamp%dExBl%dS_%s%s.tif', condition2plot{cp}, pulse_num, time_window, GroupBy{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS))
            saveas(gcf, sprintf('%s%dPul_AvgBLScatMedIQR_dFF%dSPre_By%s_ExRamp%dExBl%dS_%s%s.svg', condition2plot{cp}, pulse_num, time_window, GroupBy{cp}, ExcludeRampingCells, ExcludeBaseline, EVENTS{e}, WTS))
            close
            cd(start_dir)
            clear Scat_XLabel Scat_YLabel rho corr_pval temp_PairedData_NaNremoved temp_PairedData_selectTS temp_PairedData p_val
        end
        clear TITLE DataLabels DataLabels_PVAL adj_pval struct2plot PlotTS cond TITLE YLABEL YLIMS DataLabels DataLabels_PVAL...
            PVALS BetweenConditions_pval adj_pval v idx_1 idx_2 tickLabels adj_pval_HolmBonferroni
    end
    clear e cmap cmap_num condition
end

clear condition2plot PAIRED EXCLUDE GroupBy cp STAT cp AVG_dFF_PrePost...
    EBL ERRORBARS EVENTS ExcludeBaseline ExcludeRampingCells...
    PostStimTime pulse_cond pulse_num RampingCellThresh ResBasedOn ResNonRes...
    save_dir start_dir STAT time_window TSERIES UseAreaThresh WTS pulse_cond...
    pulse_num p_val