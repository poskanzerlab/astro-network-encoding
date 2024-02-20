%% Fig. 3b, right (overlay): representative spatial heatmap showing the 
%  change in the number of active frames/ px (post - pre). Using all events 
%  150s pre- and post-stim within the field-of-view. If landmarks were
%  drawn in AQuA, these can be drawn on top of the heatmaps.
%
%  Modified from SpatialHeatMaps_ActFramesPerPx_MC20210922.m
%  Every active pixel in each frame for each event is counted as 1. 
%  The sum of active frames for each px is calculated across the chosen time
%  window
% 
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%  3. Run the following code: time_windows = [149];
%     Note: to run all t-series for all WT, Cx43 and CBX datasets takes
%     ~15min. WT t-series 27 is the FOV shown in Fig. 3b
%
% Michelle Cahill 20240114
%%
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig3b';
fprintf(sprintf('START %s', datetime))
time_windows = [149];% the time window (s) in which to incldue events pre- and post-stim % [30 60 149];

WC = 'all'; %'uncaging' to only plot events from the uncaging cell, 'all' to plot events from all regions
ExcludeRampingCells = 1; %0 to plot all cells of that type, 1 to exclude ramping cells (determined by Max)
RampingCellThresh = 'pval_1'; %If different thresholds were originally run for ramping cells, which threshold should be used? '' if only 1 threshold was run

uncaging_point = 1; %0 no uncaging point was defined, 1: uncaging point was defined (this will plot landmark 1)
soma = 0; %0 no soma or second landmark was defined, 1: soma or second landmark was defined (this will plot landmark 2)

condition = fieldnames(mydata);

for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            xSize = mydata(tseries).(condition{cond}).res.opts.sz(2); %pixels in x (minus what AQuA takes off of borders) +1 to account for x-vals of 0
            ySize = mydata(tseries).(condition{cond}).res.opts.sz(1); %pixels in y (minus what AQuA takes off of borders) +1 to account for y-vals of 0
            ActivePxPerFrame3D = zeros(ySize, xSize, mydata(tseries).(condition{cond}).totalframes); %matrix of 0s to represent pxs in each frame 1.rows:y-coordinate, 2.columns:x-coordinate, 3.z:frames
            
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
                    if ~isempty(mydata(tseries).(condition{cond}).tLocs{1,evt})
                        clear frames zIdx ff
                        [frames, ~, zIdx] = unique(mydata(tseries).(condition{cond}).tLocs{1,evt}); % get frames at which it is active
                        for ff = 1:numel(frames) % for each frame at which that event is active
                            clear xyPairs IDX pp
                            % find all x,y pairs for that frame>
                            IDX = find(zIdx == ff);
                            xyPairs(:,1) = mydata(tseries).(condition{cond}).xLocs{1,evt}(IDX); %x
                            xyPairs(:,2) = mydata(tseries).(condition{cond}).yLocs{1,evt}(IDX); %y
                            for pp = 1:numel(xyPairs)/2
                                % sum a 1 at each x,y pair>
                                clear x y
                                x = xyPairs(pp,1);
                                y = xyPairs(pp,2);
                                ActivePxPerFrame3D(y, x, frames(ff)) = 1;
                            end
                        end
                    end
                end
            end
            clear xSize ySize  region evt frames zIdx ff xyPairs IDX pp x y ramping_cell_idx
            
            uncaging_frame = mydata(tseries).(condition{cond}).uncagingframe;
            ActivePxPerFrame3D = flip(ActivePxPerFrame3D,1);
            for time_window = time_windows
                if ~isempty(time_window) %if user designated a time window, calculate how many frames pre and post
                    num_frames = round(time_window / mydata(tseries).(condition{cond}).SecPerFrame); % convert the time window from seconds to frames
                    ActivePxPerFrame3D_pre = ActivePxPerFrame3D(:,:,(uncaging_frame-1-num_frames):(uncaging_frame-1));
                    ActPxMap_pre = sum(ActivePxPerFrame3D_pre,3);
                    clear ActivePxPerFrame3D_pre
                    ActivePxPerFrame3D_post = ActivePxPerFrame3D(:,:,uncaging_frame:(uncaging_frame+num_frames));
                    ActPxMap_post = sum(ActivePxPerFrame3D_post,3);
                    clear ActivePxPerFrame3D_post
                else
                    ActivePxPerFrame3D_pre = ActivePxPerFrame3D(:,:,1:(uncaging_frame-1));
                    ActPxMap_pre = sum(ActivePxPerFrame3D_pre,3);
                    clear ActivePxPerFrame3D_pre
                    ActivePxPerFrame3D_post = ActivePxPerFrame3D(:,:,uncaging_frame:end);
                    ActPxMap_post = sum(ActivePxPerFrame3D_post,3);
                    clear ActivePxPerFrame3D_post
                end

                ActPxMapdiffPostPre = ActPxMap_post - ActPxMap_pre;
                %         NetDiff_GABA(tseries,1) = sum(sum(ActPxMap_GABA_diffPostPre));
                 
                max_val = max(max(max(ActPxMap_pre)), max(max(ActPxMap_post)));
                clims = [0 max_val];
                figure()
                x0 = 100;
                y0 = 200;
                width = ((mydata(tseries).(condition{cond}).res.opts.sz(2) * 3) + 200)*3; %Take the width of the tseries * # of subplots + some amount of space for gaps between graphs
                height = (mydata(tseries).(condition{cond}).res.opts.sz(1)+10)*3;
                set(gcf,'position', [x0, y0, width/2, height/2])
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
                
                subplot(1,3,1)
                if sum(ActPxMap_pre, 'all') > 0
                    imagesc(ActPxMap_pre, clims)
                end
                hold on
                if uncaging_point == 1
                    plot(uc_bound_X(uc_bound_reorder), uc_bound_Y(uc_bound_reorder), 'k-', 'LineWidth', 2);
                end
                if soma == 1
                    plot(soma_bound_X(Soma_bound_reorder), soma_bound_Y(Soma_bound_reorder), 'k-', 'LineWidth', 2);
                end
                %         plot(mydata(tseries).(condition{cond}).res.ftsFilter.region.landMark.center(:,1), mydata(tseries).(condition{cond}).res.ftsFilter.region.landMark.center(:,2), 'k--')
                title(sprintf('%s TS%d %s cell(s), %d-100ms pulses, %ds pre/post: pre-stim', strrep(condition{cond}, '_', ' '), tseries, WC, mydata(tseries).(condition{cond}).num_pulses, time_window))
                % colormap limited to 20 colors>
                % col = hot(20); col(1,:) = [0 0 0]; % 20 colors
                col = flipud(autumn(20)); col(1,:) = [1 1 1];
                colormap(gca, col); c = colorbar;
                box off
                c.Label.String = '# of active frames per pixel';
                axis('square')
                
                subplot(1,3,2)
                if sum(ActPxMap_post, 'all') > 0
                    imagesc(ActPxMap_post, clims)
                end
                hold on
                if uncaging_point == 1
                    plot(uc_bound_X(uc_bound_reorder), uc_bound_Y(uc_bound_reorder), 'k-', 'LineWidth', 2);
                end
                if soma == 1
                    plot(soma_bound_X(Soma_bound_reorder), soma_bound_Y(Soma_bound_reorder), 'k-', 'LineWidth', 2);
                end
                title('post-stim')
                % colormap limited to 20 colors>
                %         col = hot(20); col(1,:) = [0 0 0]; % 20 colors
                %         col = flipud(hot(20)); col(1,:) = [0 0 0];
                col = flipud(autumn(20)); col(1,:) = [1 1 1];
                colormap(gca, col); c = colorbar;
                box off
                c.Label.String = '# of active frames per pixel';
                axis('square')
                
                subplot(1,3,3)
%                 if sum(ActPxMapdiffPostPre, 'all') > 0
                    imagesc(ActPxMapdiffPostPre, [-25, 25])
%                 end
                hold on
                if uncaging_point == 1
                    plot(uc_bound_X(uc_bound_reorder), uc_bound_Y(uc_bound_reorder), 'k-', 'LineWidth', 2);
                end
                if soma == 1
                    plot(soma_bound_X(Soma_bound_reorder), soma_bound_Y(Soma_bound_reorder), 'k-', 'LineWidth', 2);
                end
                title('change (post-stim - pre-stim)')
                % colormap limited to 20 colors>
                col = flipud(cbrewer('div', 'RdBu', 51, 'PCHIP'));
                col(26,:) = [1 1 1];
                %         col = jet(20); col(1,:) = [0 0 0]; % 20 colors
                colormap(gca, col); c = colorbar;
                box off
                c.Label.String = 'change in # of active frames per pixel';
                axis('square')
                
                
%                 cd(save_dir)
%                 if ~isempty(time_window)
%                     saveas(gcf, sprintf('%s_SpHeatMap_ActFramesPerPx_TS%d%sCell_%dpulse_TW%ds.tif', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses, time_window))
%                     saveas(gcf, sprintf('%s_SpHeatMap_ActFramesPerPx_TS%d%sCell_%dpulse_TW%ds.svg', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses, time_window))
%                 else
%                     saveas(gcf, sprintf('%s_SpHeatMap_ActFramesPerPx_TS%d%sCell_%dpulse.tif', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses))
%                     saveas(gcf, sprintf('%s_SpHeatMap_ActFramesPerPx_TS%d%sCell_%dpulse.svg', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses))
%                 end
                close
                cd(start_dir)
                
                %plot only the change (post-pre)
                x0 = 100;
                y0 = 200;
                width = ((mydata(tseries).(condition{cond}).res.opts.sz(2)*3) + 90); %Take the width of the tseries * # of subplots + some amount of space for gaps between graphs
                height = ((mydata(tseries).(condition{cond}).res.opts.sz(1)*3) +30);
                
                figure()
                set(gcf,'position', [x0, y0, width/3, height/3])
%                 if sum(ActPxMapdiffPostPre, 'all') > 0
                    imagesc(ActPxMapdiffPostPre, [-25, 25])
%                 end
                hold on
                if uncaging_point == 1
                    plot(uc_bound_X(uc_bound_reorder), uc_bound_Y(uc_bound_reorder), 'k-', 'LineWidth', 2);
                end
                if soma == 1
                    plot(soma_bound_X(Soma_bound_reorder), soma_bound_Y(Soma_bound_reorder), 'k-', 'LineWidth', 2);
                end
                title('change (post-stim - pre-stim)')
                title({sprintf('%s TS %d %d-100ms pulses %ds, %s cell', strrep(condition{cond}, '_', ' '), tseries, mydata(tseries).(condition{cond}).num_pulses, time_window, WC), 'change(post-stim - pre-stim)'})
                % colormap limited to 20 colors>
                col = flipud(cbrewer('div', 'RdBu', 51, 'PCHIP'));
                col(26,:) = [1 1 1];
                %         col = jet(20); col(1,:) = [0 0 0]; % 20 colors
                colormap(gca, col); c = colorbar;
                box off
                c.Label.String = 'change in # of active frames per pixel';
                axis('square')
                
                cd(save_dir)
                if ~isempty(time_window)
                    saveas(gcf, sprintf('%s_SpHeatMap_ChangeActFramesPerPx_TS%d%sCell_%dpulse_TW%ds.tif', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses, time_window))
                    saveas(gcf, sprintf('%s_SpHeatMap_ChangeActFramesPerPx_TS%d%sCell_%dpulse_TW%ds.svg', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses, time_window))
                    
                else
                    saveas(gcf, sprintf('%s_SpHeatMap_ChangeActFramesPerPx_TS%d%sCell_%dpulse.tif', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses))
                    saveas(gcf, sprintf('%s_SpHeatMap_ChangeActFramesPerPx_TS%d%sCell_%dpulse.svg', condition{cond}, tseries, WC, mydata(tseries).(condition{cond}).num_pulses))
                end
                close
                cd(start_dir)
                clear clims max_val x0 y0 width height col c num_frames...
                    ActivePxPerFrame3D_pre ActivePxPerFrame3D_post ActPxMap_pre ActPxMap_post ActPxMapdiffPostPre...
                    uc_bound_X uc_bound_Y uc_bound_reorder soma_bound_X soma_bound_Y Soma_bound_reorder
            end
            
            clear ActivePxPerFrame3D uncaging_frame time_window region_idx
        end
    end
end
clear start_dir save_dir uncaging_point soma tseries time_windows date cond...
    WC ExcludeRampingCells RampingCellThresh condition

fprintf(sprintf('DONE! %s', datetime))