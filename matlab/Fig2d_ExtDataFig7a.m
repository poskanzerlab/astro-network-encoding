%% Fig. 2d & Ext. Data Fig. 7a: Detected GluSnFR events following RuBi-glutamate uncaging 
% Plotting AQuA events selected as 'favorite' during event detection
% for a flexible number of conditions and pulse numbers
% From AggregatingAndPlottingFavorites_MultCond_MC20230811.m
%
% For Fig. 2d & Ext. Data Fig. 7a:
%  1. Load 'AQuA_GluSnFR_2PGluUncaging.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 1;
%     Section B: data_sets = {};
%  3. Run the following sections of code (A-C) with the variables as listed at the top of each section  
%
% Michelle Cahill 20240131
%% Section A. Separate favorite event indices by uncaging site for each tseries 
%  & record the average distance of these events from their respective uncaging 
%  sites (min and max distances)

condition = fieldnames(mydata);
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond})) && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
            fav_idx = mydata(tseries).(condition{cond}).res.btSt.evtMngrMsk; %pull the indices of favorite events
            idx_by_uncagingsite = cell(1, 6); % create a cell array to hold the indices of events that occurred at/during that uncaging episode
            distAvg_by_uncagingsite = cell(1, 6);
            maxdistAvg_by_uncagingsite = cell(1, 6);
            [~, LM_idx] = sort(mydata(tseries).(condition{cond}).res.fts.region.landMark.name); %using the landmark names, order the landmark indices so that landmark names will be in ascending order
            LM_name = 1;
            for LM = LM_idx
                evt_idx = fav_idx(mydata(tseries).(condition{cond}).res.fts.region.landmarkDist.distMin(fav_idx, LM) < 15 &... %events must be within 5um of the uncaging landmark
                    mydata(tseries).(condition{cond}).res.fts.loc.t0(1, fav_idx)'>= (mydata(tseries).(condition{cond}).uncagingframe(LM_name)-5) &... %events must start no sooner than 5 frames before the uncaging frame
                    mydata(tseries).(condition{cond}).res.fts.loc.t0(1, fav_idx)'<= (mydata(tseries).(condition{cond}).uncagingframe_end(LM_name)+5)); %events must start no later than 5 frames after the end of uncaging
                idx_by_uncagingsite{1, LM_name} = evt_idx;
                distAvg_by_uncagingsite{1, LM_name} = mydata(tseries).(condition{cond}).res.fts.region.landmarkDist.distAvg(evt_idx, LM);
                maxdistAvg = cellfun(@(x) nanmean(x,1), mydata(tseries).(condition{cond}).res.fts.region.landmarkDist.distMaxPerFrame, 'UniformOutput' ,0); %average the max dist across frames for all LM and all evts
                maxdistAvg_LM = cell2mat(cellfun(@(x) x(1,LM), maxdistAvg, 'UniformOutput', 0)); %For all evts, pull out the value for the LM in question
                maxdistAvg_by_uncagingsite{1, LM_name} = maxdistAvg_LM(1, evt_idx); %Take the values for only the evts that are identified for that LM/uncaging site/period
                LM_name = LM_name + 1;
                clear evt_idx maxdistAvg maxdistAvg_LM
            end
            mydata(tseries).(condition{cond}).EvtIdxByLM = idx_by_uncagingsite;
            mydata(tseries).(condition{cond}).AvgDistByLM = distAvg_by_uncagingsite;
            mydata(tseries).(condition{cond}).AvgMaxDistByLM = maxdistAvg_by_uncagingsite;
            clear LM_name LM LM_idx fav_idx idx_by_uncagingsite evt_idx distAvg_by_uncagingsite maxdistAvg_by_uncagingsite
        elseif ~isempty(mydata(tseries).(condition{cond})) %If no events were detected, but the t-series exists, record 0 events at each landmark/uncaging train
            temp_num_uncSites = mydata(tseries).(condition{cond}).num_bursts;
            mydata(tseries).(condition{cond}).EvtIdxByLM = cell(1,temp_num_uncSites);
            mydata(tseries).(condition{cond}).AvgDistByLM = cell(1,temp_num_uncSites);
            mydata(tseries).(condition{cond}).AvgMaxDistByLM = cell(1,temp_num_uncSites);
            clear temp_num_uncSites
        end
    end
    clear tseries
end

clear cond
%% Section B. Aggregate together parameter of choice for chosen pulse numbers 
% For Fig. 2d:
%    pulse_num = [10]; param_name = {'area', 'duration', 'min_distance', 'max_distance'};
% For Ext. Data Fig. 7a:
%    pulse_num = [10]; param_name = {'area', 'number_events'};

pulse_num = [10]; %List the different number of pulses per train included in these datasets  %[1 5 10];
param_name = {'area', 'duration', 'min_distance', 'max_distance', 'dff', 'number_events'}; %which event parameters to aggregate? Options include: 'area', 'duration', 'min_distance', 'max_distance', 'dff', 'number_events'

condition = fieldnames(mydata);

Fav_Events = struct();
for cond = 1:length(condition)
    for PULSE = 1:length(pulse_num)
        pulse_cond = strcat(condition{cond}, sprintf('_%d',pulse_num(PULSE)));
        if isfield(pulse_x_idx, pulse_cond)
            TSERIES = intersect(pulse_x_idx.(pulse_cond), 1:length(mydata));
            for param = 1:length(param_name)
                agg_param = [];
                for tseries = TSERIES
                    if ~isempty(mydata(tseries).(condition{cond})) && ~isempty(mydata(tseries).(condition{cond}).res.ftsFilter.basic)
                        FavEventIdx = mydata(tseries).(condition{cond}).res.btSt.evtMngrMsk;
                        if strcmp(param_name{param}, 'area')
                            FT = mydata(tseries).(condition{cond}).res.fts.basic.area(FavEventIdx);
                        elseif strcmp(param_name{param}, 'duration')
                            t0 = mydata(tseries).(condition{cond}).res.fts.loc.t0(FavEventIdx);
                            t1 = mydata(tseries).(condition{cond}).res.fts.loc.t1(FavEventIdx);
                            FT = (t1 - t0 +1) * mydata(tseries).(condition{cond}).SecPerFrame;
                            clear t0 t1
                        elseif strcmp(param_name{param}, 'min_distance')
                            FT = cell2mat(mydata(tseries).(condition{cond}).AvgDistByLM')';
                        elseif strcmp(param_name{param}, 'max_distance')
                            FT = cell2mat(mydata(tseries).(condition{cond}).AvgMaxDistByLM);
                        elseif strcmp(param_name{param}, 'dff')
                            FT = mydata(tseries).(condition{cond}).res.fts.curve.dffMax(FavEventIdx);
                        elseif strcmp(param_name{param}, 'number_events')
                            FT = cellfun(@numel, mydata(tseries).(condition{cond}).EvtIdxByLM, 'UniformOutput', 1);
                        end
                        agg_param = [agg_param, FT];
                    elseif ~isempty(mydata(tseries).(condition{cond})) && strcmp(param_name{param}, 'number_events') %Even if no events are detected, record 0 events per uncaging site
                        FT = cellfun(@numel, mydata(tseries).(condition{cond}).EvtIdxByLM, 'UniformOutput', 1);
                        agg_param = [agg_param, FT];
                    end
                    clear FT FavEventIdx
                end
                Fav_Events.(param_name{param}).(pulse_cond) = agg_param;
                clear tseries agg_param
            end
            clear TSERIES param
        end
        clear pulse_cond
    end
end
clear cond PULSE pulse_num param_name
%% Section C. Plot parameters of choice for chosen conditions as a violin plot
% comparing conditions with a rank sum test if there are two conditions and
% an ANOVA when comparing > 2 conditions
%
% For Fig. 2d:
%    param_name = {'area', 'duration', 'min_distance', 'max_distance'};
%    param_units = {'um^2', 's', 'um', 'um'};
%    YLIMITS = [0 100; 0 2; 0 10; 0 25];
%    condition2plot = {'orig'}; (this corresponds to the single round glutamate uncaging dataset)
%    EXCLUDE = {''}; COLORS = {'Greens'}; COLOR_Scheme = {'seq'};
%
% For Ext. Data Fig. 7a:
%    param_name = {'area', 'number_events'}; param_units = {'um^2', ''};
%    YLIMITS = [0 80; 0 6];
%    condition2plot = {'AllConditions'}; EXCLUDE = {''}; 
%    COLORS = {'Greens'}; COLOR_Scheme = {'seq'};
%
% Calls the functions: cbrewer.m & violinplot_STATS.m

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig7a';

param_name = {'area', 'duration', 'min_distance', 'max_distance', 'number_events'}; %'dff' %Which event parameters to plot? Each will be it's own violin plot. Options include: 'area', 'duration', 'min_distance', 'max_distance', 'dff', 'number_events'
param_units = {'um^2', 's', 'um', 'um', ''}; % What are the units for the parameters listed above? Options include: 'um^2', 's', 'um', 'um', 'dFF', ''
YLIMITS = [0 100; 0 2; 0 10; 0 25; 0 6]; % For each parameter listed above, what should the min and max y-values be for the violin plots? Each parameter is a row column 1: min y-val, column 2: max y-val

% Choose which conditions to compare: test difference between conditions
condition2plot = {'AllConditions'}; %{'WT_R1'}; % 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', '10' 'AllConditions'; Which conditions should be plotted
EXCLUDE = {''}; %{'WT'}; %Within that dataset, should any condition be left out?
COLORS = {'Greens'}; %'RdPu', 'Blues', 'YlGnBu' %'RdPu', 'Blues', 'MgCy', 'MgCy', 'MgCy'
COLOR_Scheme = {'seq'}; %'seq'

clear condition
for cp = 1:length(condition2plot)
    if strcmp(condition2plot{cp}, 'AllConditions')
        temp_condition = fieldnames(pulse_x_idx);
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
        temp_condition = fieldnames(pulse_x_idx);
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
        cmap_num = [1, 2]; %[2,1];
    end
    
    
    for param = 1:length(param_name)
        struct2plot = struct();
        for cond = 1:length(condition)
            if ~isempty(Fav_Events.(param_name{param}).(condition{cond}))
                struct2plot.(condition{cond}) = Fav_Events.(param_name{param}).(condition{cond})';
            end
        end
        
        if length(fieldnames(struct2plot))>2 %if there are more than 2 conditions with values, run an ANOVA
            STATS = 'ANOVA';
        elseif length(fieldnames(struct2plot))==2 %if there are only 2 conditions with values, run run a rank sum test
            STATS = 'RankSum';
        else %if there is only one condition, don't run any stats
            STATS = 'None';
        end
        
        TITLE = {sprintf('%s %s %s', indicator, experiment_type, strrep(condition2plot{cp}, '_', ' ')),...
            sprintf('Event %s', strrep(param_name{param}, '_', ' '))};
        YLABEL = sprintf('%s (%s)', strrep(param_name{param}, '_', ' '), param_units{param});
        YLIMS = YLIMITS(param, :);

        violinplot_STATS(struct2plot, '', TITLE, YLABEL,...
            cmap, cmap_num, STATS)
        
        if strcmp(STATS, 'ANOVA')
            subplot(4,1,1:2)
            ylim(YLIMS)
        elseif strcmp(STATS, 'RankSum')
            subplot(8,1,1:7)
            ylim(YLIMS)
        elseif strcmp(STATS, 'None')
            ylim(YLIMS)
        end

        cd(save_dir)
        saveas(gcf, sprintf('%s_%s_%s_FavoritedEvts_violin_%s.tif', indicator, condition2plot{cp}, param_name{param}, STATS))
        saveas(gcf, sprintf('%s_%s_%s_FavoritedEvts_violin_%s.svg', indicator, condition2plot{cp}, param_name{param}, STATS))
        close
        cd(start_dir)
        
        clear struct2plot PULSE TITLE YLABEL YLIMS GROUPS DataLabels v PulseNumString...
            cond idx_1 idx_2 tickLabels  
    end
    clear cmap cmap_num condition STATS
end

clear start_dir save_dir pulse_num param_name param_units YLIMITS COLOR...
    condition cmap cmap_num param cond condition2plot EXCLUDE COLORS COLOR_Scheme cp