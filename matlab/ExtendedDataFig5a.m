%% Extended Data Fig. 5a: heat maps of post-stim AQuA events color-coded by onset latency after uncaging
%  All events from the entire FOV 0-120s following uncaging will be
%  plotted, along with the uncaging site (drawn in AQuA as landmark 1)
% 
% Written for multiple conditions, excluding ramping cells and
% choosing to plot events from only the uncaging cell or from all cells
% Choose which post-stim period to include
% From SpatialMaps_AQuAEvtsColoredByOnsetLatency_MC20220819.m
%
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%  3. Run the following code:
%      time_window_start = 0; time_window_end = 120; ScalingFactor = 1;
%      WC = 'all'; ExcludeRampingCells = 1; RampingCellThresh = 'pval_1';
%      uncaging_point = 1; soma = 0; cond = 1:2
%   Note: running two conditions with 28 FOV each takes about 16min
%   WT FOV/t-series 14 is the representative FOV in Ext. Data Fig. 5a
%
% Michelle Cahill 20240115
%%
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig5a';
fprintf(sprintf('START %s', datetime))
time_window_start = 0; %time from uncaging to begin (0s would be starting from the uncaging frame)
time_window_end = 120; %[30 60 149] %time from uncaging
ScalingFactor = 1; %set the transparency (alpha) for events 0 = transparent 1 = opaque

WC = 'all'; %'uncaging' to only plot events from the uncaging cell, 'all' to plot events from all regions
ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run

uncaging_point = 1; %0 no uncaging point was defined, 1: uncaging point was defined; This will draw landmark 1 as drawn in AQuA
soma = 0; %0 no soma or second landmark was defined, 1: soma or second landmark was defined; This will draw landmark 2 as drawn in AQuA

condition = fieldnames(mydata);

for cond = 1:2%length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            %which frames correspond to the time window?
            temp_uncagingFrame = mydata(tseries).(condition{cond}).uncagingframe;
            temp_SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
            start_frame = temp_uncagingFrame + ceil(time_window_start / temp_SecPerFrame);
            end_frame = temp_uncagingFrame + ceil(time_window_end / temp_SecPerFrame);
            total_frames = numel(start_frame:end_frame);
            clear temp_uncagingFrame temp_SecPerFrame
            
            xSize = mydata(tseries).(condition{cond}).res.opts.sz(2); %pixels in x (minus what AQuA takes off of borders) +1 to account for x-vals of 0
            ySize = mydata(tseries).(condition{cond}).res.opts.sz(1); %pixels in y (minus what AQuA takes off of borders) +1 to account for y-vals of 0
            ActivePxPerFrame3D = zeros(ySize, xSize, total_frames); %matrix of 0s to represent pxs in each frame 1.rows:y-coordinate, 2.columns:x-coordinate, 3.z:frames
            
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
            
            if strcmp(WC, 'uncaging')
                region_idx = find(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.incluLmk == 1);
            elseif strcmp(WC, 'all')
                region_idx = 1:size(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx, 2);
            end
            region_idx = region_idx(~ismember(region_idx, ramping_cell_idx));
            
            for region = region_idx
                for evt = find(mydata(tseries).(condition{cond}).res.ftsFilter.region.cell.memberIdx(:,region) == 1)'
                    if (mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(evt) >= start_frame) && (mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(evt) <= end_frame)
                        %find all x,y pairs for pixels in that event
                        xyPairs(:,1) = mydata(tseries).(condition{cond}).xLocs{1,evt}; %x
                        xyPairs(:,2) = mydata(tseries).(condition{cond}).yLocs{1,evt}; %y
                        FRAME = mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0(evt) - start_frame + 1; %onset frame for the event (compared to start frame)
                        for pp = 1:numel(xyPairs)/2
                            % put a 1 at each x,y pair in the frame that
                            % the event begins
                            clear x y
                            x = xyPairs(pp,1);
                            y = xyPairs(pp,2);
                            ActivePxPerFrame3D(y, x, FRAME) = 1;
                        end
                        clear pp xyPairs FRAME x y
                    end
                end
            end
            clear xSize ySize region evt ramping_cell_idx region_idx
            
            ActivePxPerFrame3D = flip(ActivePxPerFrame3D,1);
            if uncaging_point == 1
                uc_bound_X = mydata(tseries).(condition{cond}).res.ftsFilter.region.landMark.border{1,1}(:,2); %gather x and y coordinates for the boundaries of landmark 1
                uc_bound_Y = mydata(tseries).(condition{cond}).res.ftsFilter.region.landMark.border{1,1}(:,1);
                uc_bound_reorder = boundary(uc_bound_X, uc_bound_Y, 0.1);
            end
            if soma == 1
                soma_bound_X = mydata(tseries).(condition{cond}).res.ftsFilter.region.landMark.border{2,1}(:,2); %gather x and y coordinates for the boundaries of landmark 2
                soma_bound_Y = mydata(tseries).(condition{cond}).res.ftsFilter.region.landMark.border{2,1}(:,1);
                Soma_bound_reorder = boundary(soma_bound_X, soma_bound_Y, 0.1);
            end
            
            clims = [1 size(ActivePxPerFrame3D,3)];
            col = spring(size(ActivePxPerFrame3D,3)); %cbrewer('seq', 'Oranges', size(ActivePxPerFrame3D,3)); %flipud(cbrewer('seq', 'YlOrRd', size(ActivePxPerFrame3D,3)+50)); %spring(size(ActivePxPerFrame3D,3)); %flipud(jet(size(ActivePxPerFrame3D,3))); %red for earlier frames, blue for later frames
%             col = col(51:end,:);
            MicronsPerPx = mydata(tseries).(condition{cond}).res.opts.spatialRes;
            TICKMARKS = (50/MicronsPerPx):(50/MicronsPerPx):size(ActivePxPerFrame3D,2);
            TICKMARKLABELS = {};
            for t = 1:length(TICKMARKS)
                TICKMARKLABELS = [TICKMARKLABELS, sprintf('%d', t*50)];
            end
            figure()
            colormap(col)
            hold on
            for F = size(ActivePxPerFrame3D,3):-1:1
                Slice = ActivePxPerFrame3D(:,:,F);
                imagesc(Slice*F, 'AlphaData', Slice*ScalingFactor, clims)
                clear Slice
            end
            
            if uncaging_point == 1
                plot(uc_bound_X(uc_bound_reorder), uc_bound_Y(uc_bound_reorder), 'k-', 'LineWidth', 2);
            end
            if soma == 1
                plot(soma_bound_X(Soma_bound_reorder), soma_bound_Y(Soma_bound_reorder), 'k-', 'LineWidth', 2);
            end
            c = colorbar;
            box off
            c.Label.String = 'onset latency (s)';
            c.Ticks = [1, round(size(ActivePxPerFrame3D,3)/2), size(ActivePxPerFrame3D,3)];
            c.TickLabels = {sprintf('%d', time_window_start), sprintf('%d', time_window_end/2), sprintf('%d', time_window_end)};
            xlim([0 size(ActivePxPerFrame3D,2)])
            ylim([0 size(ActivePxPerFrame3D,1)])
            set(gca, 'YDir', 'reverse')
            xticks(TICKMARKS)
            xticklabels(TICKMARKLABELS)
            xlabel('microns')
            yticks(TICKMARKS)
            yticklabels(TICKMARKLABELS)
            ylabel('microns')
            axis('square')
            title({sprintf('%s TS%d %s cell(s), %d-100ms pulses', strrep(condition{cond}, '_', ' '), tseries, WC, mydata(tseries).(condition{cond}).num_pulses),...
                sprintf('%d-%ds post-stim', time_window_start, time_window_end)})
            
            cd(save_dir)
            saveas(gcf, sprintf('%s_AQuAEvtsByOnsetLatency_TS%d%sCell_%dpulse_%d_%dsPost.tif', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses, time_window_start, time_window_end))
            saveas(gcf, sprintf('%s_AQuAEvtsByOnsetLatency_TS%d%sCell_%dpulse_%d_%dsPost.svg', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses, time_window_start, time_window_end))
            close
            cd(start_dir)
            
            clear start_frame end_frame total_frames xSize ySize ActivePxPerFrame3D ramping_cell_idx region_idx...
                uc_bound_X uc_bound_Y uc_bound_reorder soma_bound_X soma_bound_Y Soma_bound_reorder clims...
                col MicronsPerPx TICKMARKS TICKMARKLABELS t F c
        end
    end
end
fprintf(sprintf('DONE! %s', datetime))

clear start_dir save_dir time_window_start time_window_end ScalingFactor WC ExcludeRampingCells...
    RampingCellThresh uncaging_point soma condition cond tseries