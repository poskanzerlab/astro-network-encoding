%% Fig 1d and Extended Data Fig. 1f: Scatter plots of event features color-coded by onset latency
% For a flexible number of conditions, features to plot on x & y and
% flexible number of treatment levels (if provided)
% modified from 'AllTS_Area_Duration_Amp_aggregated_20191031_NoReg_raincloudCDF.m'
%
%  1. Load 'AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat' workspace
%  2. Run'Fig1_PreppingDataStruct.m':
%     Section A: experiment_type = 'BathApp'; method_AgonistEntry = 'MaxCurve'
%     Section B: Run as is
%  3. Run the following code
%
% Section A: time_window = []
% Section B:
%     For Fig. 1d: time_window = []; poststim_period = []; PlotPre = 1; PlotPost = 1
%                  ft1 = 'area'; ft2 = 'duration';
%     For Extended Data Fig. 1f left: time_window = []; poststim_period = [0 60];
%                                 PlotPre = 0; PlotPost = 1; ft1 = 'area'; ft2 = 'duration';
%     For Extended Data Fig. 1f right: time_window = []; poststim_period = [150 300]; 
%                                   PlotPre = 0; PlotPost = 1; ft1 = 'area'; ft2 = 'duration';
%
% Michelle Cahill 20231213
%% Section A: calculate onset latency for each event in the given time window
condition = fieldnames(mydata);
time_window = []; %time in seconds to plot events pre and post-stim, [] for all events pre- and post-stim
TW = sprintf('TimeWindow%ds', time_window);

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            if isempty(time_window)
                event_idx_pre = mydata(tseries).(condition{cond}).event_idx_pre;
                event_idx_post = mydata(tseries).(condition{cond}).event_idx_post;
            else
                event_idx_pre = mydata(tseries).(condition{cond}).(TW).event_idx_pre;
                event_idx_post = mydata(tseries).(condition{cond}).(TW).event_idx_post;
            end
            
            pre_latency = nan(1, length(event_idx_pre));
            post_latency = nan(1, length(event_idx_post));
            
            for evt = event_idx_pre
                evt_idx = find(event_idx_pre == evt);
                if ~isempty(mydata(tseries).(condition{cond}).tLocs{1, evt})
                    pre_latency(1, evt_idx) = min(mydata(tseries).(condition{cond}).tLocs{1,evt}) * mydata(tseries).(condition{cond}).SecPerFrame; %latency from start of tseries
                end
                clear evt_idx
            end
            clear evt
            
            for evt = event_idx_post
                evt_idx = find(event_idx_post == evt);
                if ~isempty(mydata(tseries).(condition{cond}).tLocs{1, evt})
                    post_latency(1, evt_idx) = ((min(mydata(tseries).(condition{cond}).tLocs{1,evt}))- mydata(tseries).(condition{cond}).uncagingframe)...
                        * mydata(tseries).(condition{cond}).SecPerFrame; %onset latency from stim
                end
                clear evt_idx
            end
            clear evt
            
            if isempty(time_window)
                mydata(tseries).(condition{cond}).latency_PreVPost = {pre_latency, post_latency};
            else
                mydata(tseries).(condition{cond}).(TW).latency_PreVPost = {pre_latency, post_latency};
            end
            clear pre_latency post_latency event_idx_pre event_idx_post
        end
    end
end
clear tseries cond
%% Section B: For 2 chosen features in a specific time period, aggregate together events and create a scatter plot
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig1D';

condition = fieldnames(mydata);
TreatmentVals = [5 25 50 100];
% condition = fieldnames(mydata);

% determine which features to plot
ft1 = 'area'; %'area', 'duration', 'dffMax', 'propGrowOverall'
ft2 = 'duration';
PlotPre = 0;
PlotPost = 1;
time_window = []; %time in seconds to plot events pre and post-stim, [] for all events pre- and post-stim
TW = sprintf('TimeWindow%ds', time_window);
poststim_period = []; %time in seconds to plot events post-stim (uses full color spectrum based on time_window variable), [] for all events in time_window

% 1. For the desired features aggregate together events from the same
%    condition/treatment level
ft1_name = sprintf('%s_PreVPost', ft1);
ft2_name = sprintf('%s_PreVPost', ft2);
fts_pre_post = struct();
fts_pre_post_sorted = struct();
MaxLatencyPost = [];
for cond = 1:length(condition)
    if isempty(TreatmentVals)
        fts_pre_post.(condition{cond}).pre = [];
        fts_pre_post.(condition{cond}).post = [];
        for tseries = 1:length(mydata)
            if ~isempty(mydata(tseries).(condition{cond}))
                if isempty(time_window)
                    temp_fts_pre = nan(length(mydata(tseries).(condition{cond}).event_idx_pre), 3); %rows are events, column 1: ft1, column2:ft2, column3: latency
                    temp_fts_post = nan(length(mydata(tseries).(condition{cond}).event_idx_post), 3);
                    
                    temp_fts_pre(:, 1) = mydata(tseries).(condition{cond}).(ft1_name){1,1}';
                    temp_fts_pre(:, 2) = mydata(tseries).(condition{cond}).(ft2_name){1,1}';
                    temp_fts_pre(:, 3) = mydata(tseries).(condition{cond}).latency_PreVPost{1,1}';
                    
                    temp_fts_post(:, 1) = mydata(tseries).(condition{cond}).(ft1_name){1,2}';
                    temp_fts_post(:, 2) = mydata(tseries).(condition{cond}).(ft2_name){1,2}';
                    temp_fts_post(:, 3) = mydata(tseries).(condition{cond}).latency_PreVPost{1,2}';
                else
                    temp_fts_pre = nan(length(mydata(tseries).(condition{cond}).(TW).event_idx_pre), 3); %rows are events, column 1: ft1, column2:ft2, column3: latency
                    temp_fts_post = nan(length(mydata(tseries).(condition{cond}).(TW).event_idx_post), 3);
                    
                    temp_fts_pre(:, 1) = mydata(tseries).(condition{cond}).(TW).(ft1_name){1,1}';
                    temp_fts_pre(:, 2) = mydata(tseries).(condition{cond}).(TW).(ft2_name){1,1}';
                    temp_fts_pre(:, 3) = mydata(tseries).(condition{cond}).(TW).latency_PreVPost{1,1}';
                    
                    temp_fts_post(:, 1) = mydata(tseries).(condition{cond}).(TW).(ft1_name){1,2}';
                    temp_fts_post(:, 2) = mydata(tseries).(condition{cond}).(TW).(ft2_name){1,2}';
                    temp_fts_post(:, 3) = mydata(tseries).(condition{cond}).(TW).latency_PreVPost{1,2}';
                end
                fts_pre_post.(condition{cond}).pre = [fts_pre_post.(condition{cond}).pre; temp_fts_pre];
                fts_pre_post.(condition{cond}).post = [fts_pre_post.(condition{cond}).post; temp_fts_post];
                clear temp_fts_pre temp_fts_post
            end
        end
        [fts_pre_post_sorted.(condition{cond}).post, ~] = sortrows(fts_pre_post.(condition{cond}).post, 3, {'ascend'}); %sort events by onset time
        MaxLatencyPost = [MaxLatencyPost, max(fts_pre_post.(condition{cond}).post(:,3))];
    else
        for TV = TreatmentVals
            concentration_TS_idx = sprintf('%s%d',condition{cond}, TV);
            fts_pre_post.(concentration_TS_idx).pre = [];
            fts_pre_post.(concentration_TS_idx).post = [];
            for tseries = pulse_x_idx.(concentration_TS_idx)
                if ~isempty(mydata(tseries).(condition{cond}))
                    if isempty(time_window)
                        temp_fts_pre = nan(length(mydata(tseries).(condition{cond}).event_idx_pre), 3); %rows are events, column 1: ft1, column2:ft2, column3: latency
                        temp_fts_post = nan(length(mydata(tseries).(condition{cond}).event_idx_post), 3);
                        
                        temp_fts_pre(:, 1) = mydata(tseries).(condition{cond}).(ft1_name){1,1}';
                        temp_fts_pre(:, 2) = mydata(tseries).(condition{cond}).(ft2_name){1,1}';
                        temp_fts_pre(:, 3) = mydata(tseries).(condition{cond}).latency_PreVPost{1,1}';
                        
                        temp_fts_post(:, 1) = mydata(tseries).(condition{cond}).(ft1_name){1,2}';
                        temp_fts_post(:, 2) = mydata(tseries).(condition{cond}).(ft2_name){1,2}';
                        temp_fts_post(:, 3) = mydata(tseries).(condition{cond}).latency_PreVPost{1,2}';
                    else
                        temp_fts_pre = nan(length(mydata(tseries).(condition{cond}).(TW).event_idx_pre), 3); %rows are events, column 1: ft1, column2:ft2, column3: latency
                        temp_fts_post = nan(length(mydata(tseries).(condition{cond}).(TW).event_idx_post), 3);
                        
                        temp_fts_pre(:, 1) = mydata(tseries).(condition{cond}).(TW).(ft1_name){1,1}';
                        temp_fts_pre(:, 2) = mydata(tseries).(condition{cond}).(TW).(ft2_name){1,1}';
                        temp_fts_pre(:, 3) = mydata(tseries).(condition{cond}).(TW).latency_PreVPost{1,1}';
                        
                        temp_fts_post(:, 1) = mydata(tseries).(condition{cond}).(TW).(ft1_name){1,2}';
                        temp_fts_post(:, 2) = mydata(tseries).(condition{cond}).(TW).(ft2_name){1,2}';
                        temp_fts_post(:, 3) = mydata(tseries).(condition{cond}).(TW).latency_PreVPost{1,2}';
                    end
                    fts_pre_post.(concentration_TS_idx).pre = [fts_pre_post.(concentration_TS_idx).pre; temp_fts_pre];
                    fts_pre_post.(concentration_TS_idx).post = [fts_pre_post.(concentration_TS_idx).post; temp_fts_post];
                    clear temp_fts_pre temp_fts_post
                end
            end
            [fts_pre_post_sorted.(concentration_TS_idx).post, ~] = sortrows(fts_pre_post.(concentration_TS_idx).post, 3, {'ascend'}); %sort events by onset time
            MaxLatencyPost = [MaxLatencyPost, max(fts_pre_post.(concentration_TS_idx).post(:,3))];
            clear tseries concentration_TS_idx
        end
        clear TV
    end
end
clear cond

% 2. Determine the color scale based in the timing of the last event
clear cmap
cmap = flipud(jet(ceil(max(MaxLatencyPost))));
COLOR = struct();
CONDITION = fieldnames(fts_pre_post);
for cond = 1:length(CONDITION)
    COLOR.(CONDITION{cond}) = nan(length(fts_pre_post_sorted.(CONDITION{cond}).post), 3);
    for evt = 1:length(fts_pre_post_sorted.(CONDITION{cond}).post)
        COLOR.(CONDITION{cond})(evt,:) = cmap(round(fts_pre_post_sorted.(CONDITION{cond}).post(evt,3)), :);
    end
    clear evt
end
clear cond

% 3. Create 2D scatter plot & save
fts_comb = {ft1, ft2};
fts_PlotDets = struct();
for f= 1:length(fts_comb)
    if strcmp(fts_comb{f}, 'area')
        fts_PlotDets(f).param_unit = '(um^2)'; %What are the units of this parameter? put in (), e.g. '(s)'
        fts_PlotDets(f).SCALE = 'log';
        fts_PlotDets(f).LIMS = [0 10^4]; %[-10 10]; ylims for the avg pre & post per condition/treatment level
    elseif strcmp(fts_comb{f}, 'duration')
        fts_PlotDets(f).param_unit = '(s)'; %What are the units of this parameter? put in (), e.g. '(s)'
        fts_PlotDets(f).SCALE = 'linear';
        fts_PlotDets(f).LIMS = [0 100]; %[-20 40] %[-10 10];
    elseif strcmp(fts_comb{f}, 'dffMax')
        fts_PlotDets(f).param_unit = '(dF/F max)'; %What are the units of this parameter? put in (), e.g. '(s)'
        fts_PlotDets(f).SCALE = 'linear';
        fts_PlotDets(f).LIMS = [0 35]; %[-20 40] %[-10 10];
        elseif strcmp(fts_comb{f}, 'propGrowOverall')
        fts_PlotDets(f).param_unit = '(um)'; %What are the units of this parameter? put in (), e.g. '(s)'
        fts_PlotDets(f).SCALE = 'log';
        fts_PlotDets(f).LIMS = [0 10^6]; %[-20 40] %[-10 10];
    end
end

for cond = 1:length(CONDITION)
    if ~isempty(poststim_period)
        PostStimEvts = find(fts_pre_post_sorted.(CONDITION{cond}).post(:,3) >= poststim_period(1) & fts_pre_post_sorted.(CONDITION{cond}).post(:,3) <= poststim_period(2));
    else
        PostStimEvts = 1:size(fts_pre_post_sorted.(CONDITION{cond}).post,1);
    end
    figure()
    cmap = flipud(jet(ceil(max(MaxLatencyPost))));
    if PlotPost == 1
        temp_post = scatter(fts_pre_post_sorted.(CONDITION{cond}).post(PostStimEvts,1),...
            fts_pre_post_sorted.(CONDITION{cond}).post(PostStimEvts,2), 15, COLOR.(CONDITION{cond})(PostStimEvts,:), 'filled', 'o');
        set (temp_post, 'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'none');
%         set (temp_post, 'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.65);
        hold on
    end
    if PlotPre == 1
        temp_pre = scatter(fts_pre_post.(CONDITION{cond}).pre(:,1),...
            fts_pre_post.(CONDITION{cond}).pre(:,2), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        set (temp_pre, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
%         set (temp_pre, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
    end
    
    set(gca, 'XScale', fts_PlotDets(1).SCALE);
    set(gca, 'YScale', fts_PlotDets(2).SCALE);
    xlabel(sprintf('%s %s', ft1, fts_PlotDets(1).param_unit))
    ylabel(sprintf('%s %s', ft2, fts_PlotDets(2).param_unit))
    axis([fts_PlotDets(1).LIMS fts_PlotDets(2).LIMS]) %([0 10^6 0 80]), ([0 10^6 0 50])
    axis square
    grid off
    NT = strrep(CONDITION{cond}, '_', ' ');
    if PlotPost == 1 && PlotPre == 1
        prepost = 'pre- and post-stim';
        prepost_File = 'PrePost';
    elseif PlotPost == 1 && PlotPre == 0
        prepost = 'post-stim';
        prepost_File = 'Post';
    elseif PlotPost == 0 && PlotPre == 1
        prepost = 'pre-stim';
        prepost_File = 'Pre';
    end
    title({sprintf('Cyto-GCaMP: %s %s vs %s', NT, ft1, ft2),sprintf('%s events', prepost), sprintf('time window %dS', time_window)})
    colormap(cmap)
    c = colorbar('southoutside', 'Ticks', [0,0.5,1],...
        'TickLabels', {'0',sprintf('%.0f', max(MaxLatencyPost)/2),sprintf('%.0f', max(MaxLatencyPost))});
    c.Label.String = 'event onset latency (s)';
    
    cd(save_dir)
    if isempty(poststim_period)
        saveas(gcf, sprintf('Scatter%sV%s_%s_%dS%s_0_%dSPost.tif', ft1, ft2, CONDITION{cond}, time_window, prepost_File, time_window))
        saveas(gcf, sprintf('Scatter%sV%s_%s_%dS%s_0_%dSPost.svg', ft1, ft2, CONDITION{cond}, time_window, prepost_File, time_window))
    else
        saveas(gcf, sprintf('Scatter%sV%s_%s_%dS%s_%d_%dSPost.tif', ft1, ft2, CONDITION{cond}, time_window, prepost_File, poststim_period(1), poststim_period(2)))
        saveas(gcf, sprintf('Scatter%sV%s_%s_%dS%s_%d_%dSPost.svg', ft1, ft2, CONDITION{cond}, time_window, prepost_File, poststim_period(1), poststim_period(2)))
    end
    close
    cd(start_dir)
    clear c temp_post temp_pre NT prepost prepost_File PostStimEvts
end

clear cond fts_comb fts_PlotDets f CONDITION COLOR cmap ft1 ft2 PlotPre PlotPost...
    time_window TW ft1_name ft2_name fts_pre_post fts_pre_post_sorted MaxLatencyPost...
    save_dir start_dir poststim_period condition