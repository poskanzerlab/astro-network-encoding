%% Extended Data Fig. 7b: Distance of neighboring astrocytes from the uncaging site
% Distance of active neighboring astrocytes from the uncaging site (neighboring cells that have at least 1 AQuA event within PostStimActivityWindow)
% Exclude uncaging cell and ramping cells
% Measures distance (in um) from the center of each neighboring astrocyte
% to the center of the uncaging site (currently set up for AQuA res files
% with one lankmark identified)
% From DistanceNeighboringCellsFromLmk_MC20230809.m
%
% Extended Data Fig. 7b compares the distance of neighboring astrocytes in
% two separate datasets: 1. single round of glutamate uncaging 2. multiple
% rounds of glutamate uncaging
%
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  (single round glutamate uncaging dataset is 'RuBiGlu_WT')
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%  3. Run the following sections of code:
%     A. Record the distance of each active neighboring 
%        cell from the uncaging site and the average distance per t-series:
%        condition = {'RuBiGlu_WT'};
%     B. Save the data for this dataset: running section as is
%  Clear the workspace
%
%  4. Load 'AQuA_CytoGCaMP_2PMultiRoundGluUncaging_70AU_ExtDataFig7.mat' workspace 
%     (multi-round glutamate uncaging dataset)
%  5. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B, C & D: as is
%  6. Run the following sections of code:
%     A. Record the distance of each active neighboring 
%        cell from the uncaging site and the average distance per t-series:
%        condition = {''};
%     B. Save the data for this dataset: running section as is
%
%  7. Run Section C loading 'DIST_NeighborFromLmk' structures from both datasets 

% Michelle Cahill 20240116
%% A. For each chosen condition record the distance of each active neighboring 
%     cell from the uncaging site and the average distance per t-series
%
% For single round glutamate uncaging: ExcludeRampingCells = 1;
%     RampingCellThresh = 'pval_1'; PostStimActivityWindow = 150; 
%     pulse_num = 10; condition = {'RuBiGlu_WT'};
%
% For multiple rounds glutamate uncaging: ExcludeRampingCells = 1;
%     RampingCellThresh = 'pval_1'; PostStimActivityWindow = 150; 
%     pulse_num = 10; condition = {'RuBiGlu_R1'};

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig7b';

ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run
PostStimActivityWindow = 150; %[] for all time in the recording, number in S for activity to be evaluate (all pre-stim frames & post-stim frames <= activity window)
pulse_num = 10;

condition = {'RuBiGlu_R1'}; %{'RuBiGlu_WT'}; {'RuBiGlu_R1'}; %fieldnames(mydata);
TSERIES = struct(); %Indicate the t-series from each condition to include
for cond = 1:length(condition)
    if exist('pulse_num', 'var')
        if exist('data_sets', 'var') && contains(data_sets{1,2}, 'Cx43')
            pulse_cond = strcat(extractAfter(condition{cond}, '_'), sprintf('_%d',pulse_num)); %For WT, CX43 and CBX dataset
        else
            pulse_cond = strcat(condition{cond}, sprintf('_%d',pulse_num)); %For laser uncaging control, receptor antagonist and multiple round of glutamate uncaging dataset
        end
    else
        pulse_cond = extractAfter(condition{cond}, '_');
    end
        TSERIES.(condition{cond}) = intersect(pulse_x_idx.(pulse_cond), 1:length(mydata));
clear pulse_cond    
end    
clear cond 

% For each condition record the distance of each active neighboring cell from the uncaging site and the average distance per t-series
DIST_NeighborFromLmk = struct(); 
for cond = 1:length(condition)
    ActiveNeighbors_DistMicrons_distribution = [];
    ActiveNeighbors_DistMicrons_AvgPerTS = [];
    for tseries = TSERIES.(condition{cond}) %1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            
            %1. Find the indices of neighboring regions with some activity
            UncagingID = find(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk == 1); %Identifying the uncaging cell
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
            EX_IDs = unique([UncagingID, ramping_cell_idx']); %the uncaging cell and ramping cells will be excluded
            All_reg = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk, 1);
            INCLUD_IDs = All_reg(~ismember(All_reg, EX_IDs));
            
            if isempty(PostStimActivityWindow) %To look for activity at any point during the recording
                IDRegWithEvents = find(sum(mydata(tseries).(condition{cond}).numevents_PreVPost, 2)); %Indices of regions with >= 1 AQuA event at any point during the recording
            else %To look for activity happening at or before PostStimActivityWindow cutoff (including all pre-stim time)
                PostStimActivityWindow_frames = ceil(PostStimActivityWindow / mydata(tseries).(condition{cond}).SecPerFrame);
                EvtsWithinXFrames = mydata(tseries).(condition{cond}).uncagingframe + PostStimActivityWindow_frames;
                ActiveReg = cell2mat(cellfun(@(x) sum(x <= EvtsWithinXFrames), mydata(tseries).(condition{cond}).event_times_frames, 'UniformOutput', false));
                IDRegWithEvents = find(ActiveReg); %Indices of regions with >= 1 AQuA event from the start of the recording up to PostStimeActivityWindow post-stim
                clear PostStimActivityWindow_frames EvtsWithinXFrames ActiveReg
            end
            
            %Find the overlap of region indices (neighboring cells that are
            %not excluded due to ramping 'INCLUD_IDs' and cells with some
            %events 'IDRegWithEvents'
            ActiveNeighborIDs = intersect(INCLUD_IDs, IDRegWithEvents');
            
            %2. Calculate the distance of each region centroid from the
            %uncaging centroid
            Region_Centroids = mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.center; %Regions are rows, column 1: X-coordinate, column 2: Y-coordinate
            Uncaging_Centroid = mydata(tseries).(condition{cond}).res.ftsFilter.region.landMark.center;
            Dist_Px = pdist2([Region_Centroids], [Uncaging_Centroid]);
            Dist_Microns = (Dist_Px*mydata(tseries).(condition{cond}).res.opts.spatialRes)';
            ActiveNeighbors_DistMicrons = Dist_Microns(ActiveNeighborIDs);
            
            ActiveNeighbors_DistMicrons_distribution = [ActiveNeighbors_DistMicrons_distribution, ActiveNeighbors_DistMicrons];
            ActiveNeighbors_DistMicrons_AvgPerTS = [ActiveNeighbors_DistMicrons_AvgPerTS, mean(ActiveNeighbors_DistMicrons, 'omitnan')];
            
            clear UncagingID ramping_cell_idx EX_IDs All_reg INCLUD_IDs IDRegWithEvents ActiveNeighborIDs Region_Centroids Uncaging_Centroid...
                Dist_Px Dist_Microns ActiveNeighbors_DistMicrons
        end
    end
    DIST_NeighborFromLmk.Distribution_Microns.(condition{cond}) = ActiveNeighbors_DistMicrons_distribution;
    DIST_NeighborFromLmk.AvgPerTS_Microns.(condition{cond}) = ActiveNeighbors_DistMicrons_AvgPerTS;
    clear tseries ActiveNeighbors_DistMicrons_distribution ActiveNeighbors_DistMicrons_AvgPerTS
end
clear cond
%% B. Save the data for the current dataset
if strcmp(condition{1}, 'RuBiGlu_WT') % Run this to save data for the single round of glutamate uncaging
    DIST_NeighborFromLmk_SingleRoundGluUncaging = DIST_NeighborFromLmk;
    cd(save_dir)
    save('DIST_NeighborFromLmk_SingleRoundGluUncaging.mat', 'DIST_NeighborFromLmk_SingleRoundGluUncaging')
    cd(start_dir)
elseif strcmp(condition{1},'RuBiGlu_R1') % Run this to save data for the multiple rounds of glutamate uncaging
    DIST_NeighborFromLmk_MultiRoundGluUncaging = DIST_NeighborFromLmk;
    cd(save_dir)
    save('DIST_NeighborFromLmk_MultiRoundGluUncaging.mat', 'DIST_NeighborFromLmk_MultiRoundGluUncaging')
    cd(start_dir)
end

clear ans condition DIST_NeighborFromLmk DIST_NeighborFromLmk_SingleRoundGluUncaging DIST_NeighborFromLmk_MultiRoundGluUncaging
%% C. Combine and compare distance structs from different datasets
% For Ext. Data Fig. 7b:
% 1. Load in DIST_NeighborFromLmk structs from both datasets
%    DIST_NeighborFromLmk_SingleRoundGluUncaging.mat &
%    DIST_NeighborFromLmk_MultiRoundGluUncaging.mat
% 2. Run this section with the following variables:
%    condition2plot = {'WT_R1'}; EXCLUDE = {''};
%    COLORS = {'Blues'}; COLOR_Scheme = {'seq'};
%  Single round glutamate uncaging FOVs will be labeled 'RuBiGlu WT'
%  Multi-round glutamate uncaging FOVs will be labeled 'RuBiGlu R1'
%
%  Comparisons between 2 distributions will be done using a rank sum test
%  Comparisons with > 2 distributions will be done using an ANOVA
%
%  Calls the following functions: cbrewer.m & violinplot_STATS.m

start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig7b';

%Combine the structures into a single structure
DIST_NeighborFromLmk = struct();
temp_struct = [DIST_NeighborFromLmk_SingleRoundGluUncaging, DIST_NeighborFromLmk_MultiRoundGluUncaging]; %Include the structures to combine (with same field names)
VARS = fieldnames(temp_struct);

for v = 1:length(VARS)
    DIST_NeighborFromLmk.(VARS{v}) = struct();
    for num_struct = 1:length(temp_struct)
        temp_fields = fieldnames(temp_struct(num_struct).(VARS{v}));
        for f = 1:length(temp_fields)
            DIST_NeighborFromLmk.(VARS{v}).(temp_fields{f}) = temp_struct(num_struct).(VARS{v}).(temp_fields{f})';
        end
    end
end

% Violin plots testing difference between designated conditions
% Choose which conditions to compare: test difference between conditions
% and the difference of each condition from mean = 0 
condition2plot = {'WT_R1'}; %{'WT_R1'}; This will select fields containing either 'WT' or 'R1'; Options: 'RuBiGABA', 'RuBiGlu', 'WT', 'CBX', 'Cx43fl', 'AllConditions' Which conditions should be plotted
EXCLUDE = {''}; %{'WT'}; %Within that dataset, should any condition be left out?
COLORS = {'Blues'}; %'RdPu', 'Blues', 'YlGnBu' %'RdPu', 'Blues', 'MgCy', 'MgCy', 'MgCy'
COLOR_Scheme = {'seq'}; %'seq'

clear condition
for cp = 1:length(condition2plot)
    if strcmp(condition2plot{cp}, 'AllConditions')
        temp_condition = fieldnames(DIST_NeighborFromLmk.Distribution_Microns);
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
        temp_condition = fieldnames(DIST_NeighborFromLmk.Distribution_Microns);
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
    
    if length(condition) > 2
        STATS = 'ANOVA';
    else
        STATS = 'RankSum';
    end
    
    % Violin plots of distance of each neighboring cell per
    % condition
    struct2plot = struct();
    for cond = 1:length(condition)
        struct2plot.(condition{cond}) = DIST_NeighborFromLmk.Distribution_Microns.(condition{cond});
    end
    
    YLABEL = 'distance (um)';
    TITLE = {sprintf('%s %s', indicator, experiment_type),...
        'Distance between uncaging site &'...
        sprintf('neighboring cells with activity (<= %dS Post)', PostStimActivityWindow),...
        sprintf('Exclude Ramping Cells: %d', ExcludeRampingCells)};
    
%     violinplot_ANOVA(struct2plot, '', TITLE,...
%         cmap, cmap_num, YLABEL)
    violinplot_STATS(struct2plot, '', TITLE, YLABEL,...
    cmap, cmap_num, STATS)
    
    cd(save_dir)
    saveas(gcf, sprintf('%s_Ex%s_DistUncSiteActiveNeigh_Within%dSPost_ExRamp%d_%s.tif', condition2plot{cp}, EXCLUDE{cp}, PostStimActivityWindow, ExcludeRampingCells, STATS))
    saveas(gcf, sprintf('%s_Ex%s_DistUncSiteActiveNeigh_Within%dSPost_ExRamp%d_%s.svg', condition2plot{cp}, EXCLUDE{cp}, PostStimActivityWindow, ExcludeRampingCells, STATS))
    close
    cd(start_dir)
    
    clear struct2plot cond TITLE YLABEL...
        idx_1 idx_2 tickLabels cmap cmap_num condition STATS
end

clear temp_struct VARS v num_struct temp_fields f condition2plot EXCLUDE COLORS COLOR_Scheme cp...
    save_dir start_dir
%% 
clear condition DIST_NeighborFromLmk ExcludeRampingCells PostStimActivityWindow...
    pulse_num RampingCellThresh save_dir start_dir TSERIES...
    DIST_NeighborFromLmk_SingleRoundGluUncaging DIST_NeighborFromLmk_MultiRoundGluUncaging