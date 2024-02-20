%% Event features pre- and post-stim (area, duration and propagation)
% Fig. 1f-h: dose-response curves for various event features (avg. change in feature post-pre)
% Modified from 'DistOfEvtsPrePost_SpecTimeWindow_NoRegions_MultCond_MC_20210907.m'
% Ext. Data Fig. 1g: CDFs of event features pre- and post-stim. Modified
% from 'DistributionOfBLorPostStimEvents_MultCond_MC_20230815.m

%  1. Load 'AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat' workspace
%  2. Run'Fig1_PreppingDataStruct.m':
%     Section A: experiment_type = 'BathApp'; method_AgonistEntry =
%     'MaxCurve'
%     Section B: Calculating parameters/features for population-wide
%     activity
%  3. Run the following code

% Comparing events pre- and post- agonist, population-level activity for a
% flexible number of conditions for defined time period pre- and post-stim
% and separating treatment values if desired
% A. Record event indices in pre- and post-stim time window: 
%    time_window = 120
% B. Dose response curves of designated feature of events post - pre,
%    comparing the change in response to Baclofen & t-ACPD. Calls function the
%    DoseResponseCurve.m
%    For Fig. 1f: param_name = 'Area';
%    For Fig. 1g: param_name = 'Duration';
%    For Fig. 1h: param_name = 'Growing Propagation';
% C. CDFs of chosen event features pre- and post-stim (Extended Data Fig.
%    1g). Calls function cdf_plot and breakyaxis:
%    time_window = 120;
%    param_name = {'area_PreVPost', 'duration_PreVPost', 'propGrowOverall_PreVPost'};
% Michelle Cahill 20240109
%% A. Record event indicies occurring in specified pre- and post-stim time
% window record value of individual event features occurring in specified
% time window + avgerage of each event feature pre- and post-stim
% area, duration, amplitude, propGrowOverall, propShrinkOverall
time_window = 120; %time window (in seconds) taken before the stim and after the stim
TW = sprintf('TimeWindow%ds', time_window);

condition = fieldnames(mydata);

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            clear event_idx_pre event_idx_post
            uncaging_frame = mydata(tseries).(condition{cond}).uncagingframe;
            num_frames = round(time_window / mydata(tseries).(condition{cond}).SecPerFrame); % convert the time window from seconds to frames
            
            %find the indices of events occuring pre vs post simulation
            pre_idx = [];
            post_idx = [];
            
            for event = 1:length(mydata(tseries).(condition{cond}).res.ftsFilter.basic.area) %loop through each event, assigning it to pre or post stim
                event_time = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event);
                if event_time < uncaging_frame && event_time > (uncaging_frame-1-num_frames) %break the indices into pre vs post stimulus
                    pre_idx = [pre_idx, event];
                elseif event_time >= uncaging_frame && event_time < (uncaging_frame+num_frames)
                    post_idx = [post_idx, event];
                end
                clear event_time
            end
            clear event
            event_idx_pre = pre_idx;
            event_idx_post = post_idx;
            
            mydata(tseries).(condition{cond}).(TW).event_idx_pre = event_idx_pre;
            mydata(tseries).(condition{cond}).(TW).event_idx_post = event_idx_post;
            
            %calculate number of events in time period
            eventnumber_PreVPost = nan(1,2); %rows are cells/regions, column 1 pre-stim column 2 post-stim
            eventnumber_PreVPost(1, 1) = numel(event_idx_pre);  %1st column is # of events in the time_window pre-stim
            eventnumber_PreVPost(1, 2) = numel(event_idx_post);% # of events in the time_window post-stim
            
            mydata(tseries).(condition{cond}).(TW).numevents_PreVPost = eventnumber_PreVPost;
            clear eventnumber_PreVPost
            
            %area of each event divided between pre and post stim for each cell
            area_PreVPost = cell(1, 2);
            avgarea_PreVPost = nan(1, 2);
            
            area_PreVPost(1, 1) = {mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_pre)};
            area_PreVPost(1, 2) = {mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_post)};
            avgarea_PreVPost(1, 1) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_pre));
            avgarea_PreVPost (1, 2) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.basic.area(event_idx_post));
            
            mydata(tseries).(condition{cond}).(TW).area_PreVPost = area_PreVPost;
            mydata(tseries).(condition{cond}).(TW).avgarea_PreVPost = avgarea_PreVPost;
            clear area_PreVPost avgarea_PreVPost
            
            %duration (in seconds) of each event divided between pre and post stim for each cell
            duration_PreVPost = cell(1, 2);
            avgduration_PreVPost = nan(1, 2);
            frames_pre = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t1(event_idx_pre) - mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event_idx_pre);
            duration_PreVPost(1, 1) = {frames_pre * mydata(tseries).(condition{cond}).SecPerFrame};
            avgduration_PreVPost(1, 1) = mean(frames_pre * mydata(tseries).(condition{cond}).SecPerFrame);
            
            frames_post = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t1(event_idx_post) - mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(event_idx_post);
            duration_PreVPost(1, 2) = {frames_post * mydata(tseries).(condition{cond}).SecPerFrame};
            avgduration_PreVPost(1, 2) = mean(frames_post * mydata(tseries).(condition{cond}).SecPerFrame);
            mydata(tseries).(condition{cond}).(TW).duration_PreVPost = duration_PreVPost;
            mydata(tseries).(condition{cond}).(TW).avgduration_PreVPost = avgduration_PreVPost;
            clear duration_PreVPost avgduration_PreVPost frames_post frames_pre
            
            %dffMax (amplitude) of each event divided between pre and post stim for each cell
            dffMax_PreVPost = cell(1, 2);
            avgdffMax_PreVPost = nan(1, 2);
            dffMax_PreVPost(1, 1) = {mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_pre)};
            dffMax_PreVPost(1, 2) = {mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_post)};
            avgdffMax_PreVPost(1, 1) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_pre));
            avgdffMax_PreVPost (1, 2) = mean(mydata(tseries).(condition{cond}).res.ftsFilter.curve.dffMax(event_idx_post));
            mydata(tseries).(condition{cond}).(TW).dffMax_PreVPost = dffMax_PreVPost;
            mydata(tseries).(condition{cond}).(TW).avgdffMax_PreVPost = avgdffMax_PreVPost;
            clear dffMax_PreVPost avgdffMax_PreVPost
            
            %propGrowOverall of each event divided between pre and post stim for each cell
            propGrowOverall_PreVPost = cell(1, 2);
            avgpropGrowOverall_PreVPost = nan(1, 2);
            
            propGrowOverall_PreVPost(1, 1) = {sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall(event_idx_pre,:),2)'};
            propGrowOverall_PreVPost(1, 2) = {sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall(event_idx_post,:),2)'};
            avgpropGrowOverall_PreVPost(1, 1) = mean(sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall(event_idx_pre,:),2));
            avgpropGrowOverall_PreVPost (1, 2) = mean(sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propGrowOverall(event_idx_post,:),2));
            
            mydata(tseries).(condition{cond}).(TW).propGrowOverall_PreVPost = propGrowOverall_PreVPost;
            mydata(tseries).(condition{cond}).(TW).avgpropGrowOverall_PreVPost = avgpropGrowOverall_PreVPost;
            clear propGrowOverall_PreVPost avgpropGrowOverall_PreVPost
            
            %propShrinkOverall of each event divided between pre and post stim for each cell
            propShrinkOverall_PreVPost = cell(1, 2);
            avgpropShrinkOverall_PreVPost = nan(1, 2);
            
            propShrinkOverall_PreVPost(1, 1) = {sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propShrinkOverall(event_idx_pre,:),2)'};
            propShrinkOverall_PreVPost(1, 2) = {sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propShrinkOverall(event_idx_post,:),2)'};
            avgpropShrinkOverall_PreVPost(1, 1) = mean(sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propShrinkOverall(event_idx_pre,:),2));
            avgpropShrinkOverall_PreVPost (1, 2) = mean(sum(mydata(tseries).(condition{cond}).res.ftsFilter.propagation.propShrinkOverall(event_idx_post,:),2));
            
            mydata(tseries).(condition{cond}).(TW).propShrinkOverall_PreVPost = propShrinkOverall_PreVPost;
            mydata(tseries).(condition{cond}).(TW).avgpropShrinkOverall_PreVPost = avgpropShrinkOverall_PreVPost;
            clear propShrinkOverall_PreVPost avgpropShrinkOverall_PreVPost
            
            clear uncaging_frame num_frames pre_idx post_idx event_idx_pre event_idx_post eventnumber_PreVPost
        end
    end
end
clear time_window TW cond tseries
%% Define the time window to use
% condition = {'Baclofen', 'tACPD'};
% TreatmentVals = [5 25 50 100];
time_window = 120; %time in seconds for pre- and post-stim time window, [] for all time pre- and post-stim
if ~isempty(time_window)
    TW = sprintf('TimeWindow%ds', time_window);
end
%% B. Plotting dose-response curves for designated feature/parameter, combining across regions
% Fig. 1f-h
% 1. For each condition & each treatment level, calculate avg change for
%    chosen event feature per t-series/recording
% 2. Dose reponse curves for specified conditions avg post - avg pre of an
%    event feature in the specified time window (t-test & Holm-Bonferroni
%    correction)
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig1F_H';

param_name = 'Growing Propagation'; % What parameter to plot? 'Area', 'Duration', 'Amplitude', 'Growing Propagation' or 'Shrinking Propagation'
ERRORBARS = 'sem'; %'sem' or 'std'

if strcmp(param_name, 'Area')
    param = 'area_PreVPost'; %Which parameter would you like to aggregate?
    param_avg = 'avgarea_PreVPost';
    param_unit = '(um^2)'; %What are the units of this parameter? put in (), e.g. '(s)'
    YLIMS = [0 600]; %[-10 10]; ylims for the avg pre & post per condition/treatment level
    DR_YLIMS = [-200 300]; %ylims for the dose response curves
    BIN_WIDTH = 10;
elseif strcmp(param_name, 'Duration')
    param = 'duration_PreVPost'; %Which parameter would you like to aggregate?
    param_avg = 'avgduration_PreVPost';
    param_unit = '(s)'; %What are the units of this parameter? put in (), e.g. '(s)'
    YLIMS = [0 20]; %[-20 40] %[-10 10];
    DR_YLIMS = [-5 7]; %[-5 10] ylims for the dose response curves
    BIN_WIDTH = 1;
elseif strcmp(param_name, 'Amplitude')
    param = 'dffMax_PreVPost'; %Which parameter would you like to aggregate?
    param_avg = 'avgdffMax_PreVPost';
    param_unit = '(dF/F max)'; %What are the units of this parameter? put in (), e.g. '(s)'
    YLIMS = [0 10]; %[-20 20];
    DR_YLIMS = [-2 4]; %[-5 10] ylims for the dose response curves
    BIN_WIDTH = 0.25;
elseif strcmp(param_name, 'Growing Propagation')
    param = 'propGrowOverall_PreVPost'; %Which parameter would you like to aggregate?
    param_avg = 'avgpropGrowOverall_PreVPost';
    param_unit = '(um)'; %What are the units of this parameter? put in (), e.g. '(s)'
    YLIMS = [0 30]; %[-20 20];
    DR_YLIMS = [-10 20]; %[-5 10] ylims for the dose response curves
    BIN_WIDTH = 1;
elseif strcmp(param_name, 'Shrinking Propagation')
    param = 'propShrinkOverall_PreVPost'; %Which parameter would you like to aggregate?
    param_avg = 'avgpropShrinkOverall_PreVPost';
    param_unit = '(um)'; %What are the units of this parameter? put in (), e.g. '(s)'
    YLIMS = [-30 0]; %[-20 20];
    DR_YLIMS = [-20 10]; %[-4000 2000] ylims for the dose response curves
    BIN_WIDTH = 1;
end

CMAP.Baclofen = [0.5 0.5 .5; 1 0 1]; %grey & magenta
CMAP.tACPD = [0.5 0.5 .5; 0 1 1]; %grey & cyan
CMAP.Baclofen_CBX = [0.5 0.5 .5; 1 0 1]; %grey & magenta
CMAP.tACPD_CBX = [0.5 0.5 .5; 0 1 1]; %grey & cyan

% 1. Calculate avg change for chosen event feature for each t-series 
ChangeEvtFt = struct(); %record (avg post - avg pre) for each tseries
for cond = 1:length(condition)
    tsavg_comb = nan(length(mydata), 2);
    ts_comb = struct();
    ts_comb.pre = [];
    ts_comb.post = [];
    
    %average of evt feature pre- and post-stim & average change post-pre
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond})) &&...
                ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
            if isempty(time_window)
                tsavg_comb(tseries,:) = mydata(tseries).(condition{cond}).(param_avg);
            else
                tsavg_comb(tseries,:) = mydata(tseries).(condition{cond}).(TW).(param_avg);
            end
        end
    end
    ChangeEvtFt.(condition{cond}) = tsavg_comb(:,2) - tsavg_comb(:,1);
end
clear CMAP cond

% 2. Dose response curves
if ~isempty(TreatmentVals)
    condition2plot = {{'Baclofen', 'tACPD'}};
    TTEST_Type = {'paired'};
    for cp = 1:length(condition2plot)
        CONDITION = condition2plot{cp}; %{'Baclofen', 'tACPD'}; %{'Baclofen_CBX', 'tACPD_CBX'}
        TITLE = {'Cyto-GCaMP 16X Bath App:',sprintf('Event %s', param_name),...
            sprintf('%dS pre & post', time_window)};
        YLABEL = strcat('\Delta in mean event', sprintf('%s', lower(param_name)));
        %         YLIMS = [-200 300];
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
        DoseResponseCurves(ChangeEvtFt, pulse_x_idx, TreatmentVals, CONDITION,...
            TITLE, YLABEL, DR_YLIMS, XLABEL, XLIMS,...
            DATACOLORS, ERRORBARS, 2, 1, TTEST_Type{cp})
        
        WhichDataSets = string(strrep(join(CONDITION), ' ', ''));
        cd(save_dir)
        saveas(gcf, sprintf('DoseResponseEvt%s_%dSPrePost_%s.tif', param_name, time_window, WhichDataSets))
        saveas(gcf, sprintf('DoseResponseEvt%s_%dSPrePost_%s.svg', param_name, time_window, WhichDataSets))
        close
        cd(start_dir)
        clear CONDITION TITLE YLABEL XLABEL XLIMS DATACOLORS WhichDataSets cond
    end
end
clear cp condition2plot ChangeEvtFt

clear param param_avg param_name param_unit YLIMS DR_YLIMS TTEST_Type TV BIN_WIDTH...
    ans ERRORBARS save_dir start_dir time_window ts_comb tsavg_comb tseries TW
%% C. CDFs of chosen event features pre- and post-stim: Extended Data Fig. 1g
% 1. Aggregating together chosen event features across all
% recordings for each conidition and treatment value
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig1g';

condition = fieldnames(mydata);
BaselineComp = struct();
pulse_num = [5 25 50 100];
time_window = 120; % time in seconds to pull events from pre- and post-stim, [] to plot all events
if ~isempty(time_window)
    TW = sprintf('TimeWindow%ds', time_window);
end

param_name = {'area_PreVPost', 'duration_PreVPost', 'propGrowOverall_PreVPost'}; %,'numevents_PreVPost',  'propGrowOverall_PreVPost', 'propShrinkOverall_PreVPost', 'dffMax_PreVPost';
for param = 1:length(param_name)
    for cond = 1:length(condition)
        for TV = pulse_num
            concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
            temp_pre = [];
            temp_post = [];
            for tseries = pulse_x_idx.(concentration_TS_idx)
                if ~isempty(mydata(tseries).(condition{cond}))
                    if isempty(time_window)
                        param_PrePost = mydata(tseries).(condition{cond}).(param_name{param});
                    else
                        param_PrePost = mydata(tseries).(condition{cond}).(TW).(param_name{param});
                    end
                    if iscell(param_PrePost)
                        param_PrePost_vert = cellfun(@(x) reshape(x, [ ], 1), param_PrePost, 'UniformOutput', 0); %for each region list the event params in a column
                        pre_acrossRegions = cell2mat(param_PrePost_vert(:,1))'; % combine all regions into a single array- pre-stim events
                        post_acrossRegions = cell2mat(param_PrePost_vert(:,2))'; % combine all regions into a single array- post-stim events
                    else
                        param_PrePost_vert = param_PrePost;
                        pre_acrossRegions = param_PrePost_vert(:,1)'; % combine all regions into a single array- pre-stim events
                        post_acrossRegions = param_PrePost_vert(:,2)'; % combine all regions into a single array- post-stim events
                    end
                    temp_pre = [temp_pre, pre_acrossRegions(~isinf(pre_acrossRegions)&~isnan(pre_acrossRegions))]; %take only baseline events
                    temp_post = [temp_post,  post_acrossRegions(~isinf(post_acrossRegions)&~isnan(post_acrossRegions))]; %take only post-stim events
                    clear param_PrePost param_PrePost_vert pre_acrossRegions post_acrossRegions
                end
            end
            BaselineComp.Pre.(param_name{param}).(concentration_TS_idx) = temp_pre';
            BaselineComp.Post.(param_name{param}).(concentration_TS_idx) = temp_post';
            clear temp_pre temp_post concentration_TS_idx
        end
    end
end
clear cond tseries

% 2. Plot pre- and post-stim event parameter distributions as violin plots 
% and cdfs by treatment val/pulse number with all conditions on the same axes
% compare distributions with 1-way ANOVA
PLOT_BASELINE_TOGETHER = 0; %To decide the order of violin plots (1 to plot all pre-stim distributions first, 0 to plot pre- and post-stim distributions for each condition together)
STATS = 'ANOVA'; %'ANOVA' or 'RankSum'

for param = 1:length(param_name)
    if strcmp(param_name{param}, 'ActivePxPerFrame_NormAllPxFOV')
        YLABEL = '% px active/all px in FOV';
        BREAK_Y_AXES = [];
    elseif contains(param_name{param}, 'num')
        YLABEL = '# events/min';
        BREAK_Y_AXES = [];
    elseif contains(param_name{param}, 'area')
        YLABEL = 'event area (um^2)';
        YLIMS = [0 1000];
        BREAK_Y_AXES = [0.1 0.6];
    elseif contains(param_name{param}, 'duration')
        YLABEL = 'event duration (s)';
        YLIMS = [0 30];
        BREAK_Y_AXES = [0.1 0.3];
    elseif contains(param_name{param}, 'dff')
        YLABEL = 'event amplitude (dF/F max)';
        YLIMS = [0 15];
        BREAK_Y_AXES = [];
    elseif contains(param_name{param}, 'propGrowOverall')
        YLABEL = 'event propagation distance (um)';
        YLIMS = [0 25];
        BREAK_Y_AXES = [0.1 0.6];
    end
    for TV = pulse_num 
        temp_condition = fieldnames(BaselineComp.Pre.(param_name{param}));
        COND_TV = [sprintf("Baclofen%d", TV), sprintf("tACPD%d", TV)];
        idx = find(cell2mat(cellfun(@(x) endsWith(x, COND_TV), temp_condition, 'UniformOutput', 0)));
        condition = temp_condition(idx, 1);
        
        struct2plot = struct();
        if PLOT_BASELINE_TOGETHER
            for cond = 1:length(condition) %plot pre-stim distributions first 
                temp_FNPre = sprintf('Pre_%s', condition{cond});
                struct2plot.(temp_FNPre) = BaselineComp.Pre.(param_name{param}).(condition{cond});
            end
            for cond = 1:length(condition) %plot pre- and post-stim distributions next to one another
                temp_FNPost = sprintf('Post_%s', condition{cond});
                struct2plot.(temp_FNPost) = BaselineComp.Post.(param_name{param}).(condition{cond});
            end
        else
            for cond = 1:length(condition) %plot pre- and post-stim distributions together
                temp_FNPre = sprintf('Pre_%s', condition{cond});
                struct2plot.(temp_FNPre) = BaselineComp.Pre.(param_name{param}).(condition{cond});
                temp_FNPost = sprintf('Post_%s', condition{cond});
                struct2plot.(temp_FNPost) = BaselineComp.Post.(param_name{param}).(condition{cond});
            end
        end
        clear temp_condition idx COND_TV cond temp_FNPost temp_FNPre
        
        TITLE = {sprintf('%s: %s %duM', strrep(indicator, '_', ' '), erase(strrep(param_name{param}, '_', ' '), 'PreVPost'), TV),...
            sprintf('%dS Pre & Post', time_window)};
        
        %a. create violin plots of the distributions
        violinplot_STATS(struct2plot, [], TITLE, YLABEL,...
            [], [], STATS)
        
        cd(save_dir)
        saveas(gcf, sprintf('%s%duM_%dSPrePostEvent%s_violinplot_%s.tif', indicator, TV, time_window, param_name{param}, STATS))
        saveas(gcf, sprintf('%s%duM_%dSPrePostEvent%s_violinplot_%s.svg', indicator, TV, time_window, param_name{param}, STATS))
        close
        cd(start_dir)
        
        XLABEL = YLABEL;
        temp_YLABEL = 'probability';
        
        temp_FN = fieldnames(struct2plot);
        DATA_LABELS = cellfun(@(x) strrep(x, '_', ' '), temp_FN, 'UniformOutput', 0);
        COLORS = [201/255 148/255 199/255;...
            221/255 28/255 119/255;...
            161/255 217/255 155/255;...
            49/255 163/255 84/255]; %pink, raspberry, light green, dark green
        
        % plot cdfs
        struct2plot_cellarray = cell(1,length(temp_FN));
        for f = 1:length(temp_FN)
            struct2plot_cellarray{1, f} = struct2plot.(temp_FN{f});
        end
        cdf_plot(struct2plot_cellarray, COLORS, XLABEL, TITLE, DATA_LABELS)
        % uncomment this if statement to break y-axis (this hides the
        % figure legend)
%         if ~isempty(BREAK_Y_AXES)
%             ylim([0 1])
%             breakyaxis(BREAK_Y_AXES)
%         end
        
        cd(save_dir)
        saveas(gcf, sprintf('%s%duM_%dSPrePostEvent%s_CDF_%s.tif', indicator, TV, time_window, param_name{param}, STATS))
        saveas(gcf, sprintf('%s%duM_%dSPrePostEvent%s_CDF_%s.svg', indicator, TV, time_window, param_name{param}, STATS))
        close
        cd(start_dir)
     
        clear strut2plot_Full time_period condition struct2plot cond TITLE...
            GROUPINGS STRUC2PLOT_array XLIMS BINWIDTH BW DATA_LABELS XLABEL temp_YLABEL COLORS...
            struct2plot_cellarray f temp_FN
    end
    clear TV YLABEL
end
clear start_dir save_dir PlotBaseline param_name STATS cmap cmap_num...
    PLOT_BASELINE_TOGETHER STATS param h ans BaselineComp BREAK_Y_AXES...
    time_window TW YLIMS