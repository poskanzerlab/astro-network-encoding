function [XSecBinsPrePost_PermMean,XSecBinsChangeAvgBaseline_PermMean, XSecBinsChangeAvgBLPrePost_PermMean] = MeanFt_XSecStatBins_PermEventOnsetByRegion_NoLmk(mydata, CONDITION, TSERIES,...
    time_window, ExcludeBaseline, perm_num, param_name, ExcludeRampingCells, RampingCellThresh, NBinsPre, NBinsPost)
% MeanFt_XSecStatBins_PermEventOnsetByRegion_NoLmk performs a permutation test 
% for n-permutations (perm_num) on the chosen feature/parameter. 
% Currently written to shuffle events placed in x-sec time bins per region 
% in 'Aggregate_MC_ftsFilter_MultCond_XSecBins_20220524.m' and give a
% permuted mean per time bin for a chosen feature.
% For each tseries, the average of the chosen feature is calculated in 
% Xsec bins (determined by the variable time_window) for each chosen region.
% The average across all regions with some activity is reported.

% OUTPUTS:
% The function will return 2 matrices:
% XSecBinsPrePost_PermMean: rows will be the mean across t-series for one round of
% permutation, columns will be x-second bins both pre- and post-stim. N
% rows = perm_num
% XSecBinsChangeAvgBaseline_PermMean: rows will be the mean across t-series 
% for one round of permutation, columns will be average change of post-stim 
% x-second bins from average baseline (change calculated by region before 
% averaging across regions & t-series). N rows = perm_num
% XSecBinsChangeAvgBLPrePost_PermMean: rows will be the mean across t-series 
% for one round of permutation, columns will be average change of pre- and post-stim 
% x-second bins from average baseline (change calculated by region before 
% averaging across regions & t-series). N rows = perm_num

%INPUTS
% mydata: structed array with AQuA results files aggregated
% CONDITION: the specific condition/fieldname of mydata to run the permutation test
% for (ex: 'RuBiGABA_WT' for mydata.(CONDITION))
% TSERIES: indices of t-series to average across
% time_window: Duration of each x-sec static bin
% ExcludeBaseline: number of seconds to exclude from the start of each
% recording (should be a multiple of 'time_window')
% perm_num: number of permutations to carry out
% param_name: string to denote which feature/parameter to average:
% 'Number', 'Area', 'Duration' or 'Amplitude' (may want to rethink
% shuffling method for area, duration and amplitude to leave events in the
% original xsec bin and just shuffle the event parameter from events in
% that region)
% ExcludeRampingCells: 1 to exclude any cells identified, 0 to use all
% cells
% RampingCellThresh: if ExcludeRampingCells == 1, which threshold should be
% used to define ramping cells: 'pval_05', 'pval_1', 'pval_2'
% NBinsPre: number of pre-stim x-sec bins to include 
% NBinsPost: number of post-stim x-sec bins to include

% Michelle Cahill 20220317

XSecBinsPrePost_PermMean = NaN(perm_num, NBinsPre+NBinsPost);
XSecBinsChangeAvgBaseline_PermMean = NaN(perm_num, NBinsPost);
XSecBinsChangeAvgBLPrePost_PermMean = NaN(perm_num, NBinsPre+NBinsPost);
TB = sprintf('StaticTimeBins%ds', time_window);

for perm = 1:perm_num
    XSecBinsPrePostPerm = NaN(length(TSERIES), NBinsPre+NBinsPost); %average value per X-sec bins per t-series for this round of permutation
    XSecBinsChangeAvgBaselinePerm = NaN(length(TSERIES), NBinsPost); %average change from baseline per X-sec post-stim bins per t-series for this round of permutation
    XSecBinsChangeAvgBLPrePostPerm = NaN(length(TSERIES), NBinsPre+NBinsPost); %average change from baseline per X-sec pre- and post-stim bins per t-series for this round of permutation
    TS = 1;
    for tseries = TSERIES %For each t-series randomly assign each event to an Xsec bin for each region
        if ~isempty(mydata(tseries).(CONDITION))
            ExcludeBaselineBins = ExcludeBaseline / time_window;
            
            if ExcludeRampingCells
                if isfield(mydata(tseries).(CONDITION), 'RampingCellIdx')
                    if isempty(RampingCellThresh)
                        ramping_cell_idx = mydata(tseries).(CONDITION).RampingCellIdx;
                    else
                        ramping_cell_idx = mydata(tseries).(CONDITION).RampingCellIdx.(RampingCellThresh);
                    end
                else
                    ramping_cell_idx = [];
                end
            else
                ramping_cell_idx = [];
            end

            Cell_idx = 1:size(mydata(tseries).(CONDITION).res.ftsFilter.region.cell.memberIdx,2);
            Cell_idx = Cell_idx(~ismember(Cell_idx, ramping_cell_idx)); %exclude any ramping cells from the group (if ExcludeRampingCells = 1)
            EvtsPre = cell2mat(cellfun(@(x) ~isempty(x), mydata(tseries).(CONDITION).(TB).event_idx_pre, 'UniformOutput', 0));
            EvtsPost = cell2mat(cellfun(@(x) ~isempty(x), mydata(tseries).(CONDITION).(TB).event_idx_post, 'UniformOutput', 0));
            sum_evts = sum([EvtsPre,EvtsPost],2); %Sum the events pre & post in each cell
            Nonzero_cells = Cell_idx(sum_evts(Cell_idx) > 0); %Find the indices of neighboring cells with at least one event
            clear Cell_idx sum_events
            
            if ~isempty(Nonzero_cells)
                EVENT_IDX_PRE = mydata(tseries).(CONDITION).(TB).event_idx_pre(Nonzero_cells,ExcludeBaselineBins+1:end);
                EVENT_IDX_POST = mydata(tseries).(CONDITION).(TB).event_idx_post(Nonzero_cells,:);
                NumPreBins = size(EVENT_IDX_PRE, 2);
                NumPostBins = size(EVENT_IDX_POST, 2);
                TotalBins = NumPreBins + NumPostBins;
                EVENT_IDX_ALL = [EVENT_IDX_PRE, EVENT_IDX_POST];
                
                EVENT_IDX_PERM = cell(size(EVENT_IDX_ALL));
                for region = 1:size(EVENT_IDX_PERM, 1)
                    EVENT_IDX_COMB = cat(2, EVENT_IDX_ALL{region,:});
                    AssignBins = randi(TotalBins, 1, length(EVENT_IDX_COMB));
                    for BIN = 1:TotalBins
                        temp_IDX = EVENT_IDX_COMB(AssignBins == BIN);
                        EVENT_IDX_PERM{region, BIN} = temp_IDX;
                        clear temp_IDX
                    end
                    clear AssignBins EVENT_IDX_COMB BIN
                end
                clear region
                
                %Using the shuffled event indices, pull out the parameter of
                %interest and calculate the average per xsec bin for each
                %region
                if strcmp(param_name, 'Number')
                    EventFt_Perm = cell2mat(cellfun(@(x) numel(x), EVENT_IDX_PERM, 'UniformOutput', 0));
                elseif strcmp(param_name, 'Area')
                    Area_PerBin = cellfun(@(x) mydata(tseries).(CONDITION).res.ftsFilter.basic.area(x), EVENT_IDX_PERM, 'UniformOutput', 0);
                    EventFt_Perm = cell2mat(cellfun(@(x) mean(x, 'omitnan'), Area_PerBin, 'UniformOutput', 0));
                elseif strcmp(param_name, 'Duration')
                    duration = (mydata(tseries).(CONDITION).res.ftsFilter.loc.t1 - mydata(tseries).(CONDITION).res.ftsFilter.loc.t0) * mydata(tseries).(CONDITION).SecPerFrame;
                    Duration_PerBin = cellfun(@(x) duration(x), EVENT_IDX_PERM, 'UniformOutput', 0);
                    EventFt_Perm = cell2mat(cellfun(@(x) mean(x, 'omitnan'), Duration_PerBin, 'UniformOutput', 0));
                elseif strcmp(param_name, 'Amplitude')
                    dffMax_PerBin = cellfun(@(x) mydata(tseries).(CONDITION).res.ftsFilter.curve.dffMax(x), EVENT_IDX_PERM, 'UniformOutput', 0);
                    EventFt_Perm = cell2mat(cellfun(@(x) mean(x, 'omitnan'), dffMax_PerBin, 'UniformOutput', 0));
                elseif strcmp(param_name, 'PropFrequency')
                    Sum_GrowingProp = sum(mydata(tseries).(CONDITION).res.ftsFilter.propagation.propGrowOverall,2)';
                    GrowingProp_PerBin = cellfun(@(x) Sum_GrowingProp(x), EVENT_IDX_PERM, 'UniformOutput', 0);
                    EventFt_Perm  = cell2mat(cellfun(@(x) sum(x>1), GrowingProp_PerBin, 'UniformOutput', 0));
                end
                
                % Adjust for differences in pre and post bin number and
                % shift bins for pre- and post- to line up across t-series
                Pre_Size_Diff = NBinsPre - NumPreBins; %number of pre bins to include - number of pre bins in this t-series
                Post_Size_Diff = NBinsPost - NumPostBins; %number of post bins to include - number of post bins in this t-series
                
                AVG_PRE_ByRegion = mean(EventFt_Perm(:, 1:NumPreBins), 2, 'omitnan'); %average value at baseline for each region
                ChgAvgBL_PrePost = mean((EventFt_Perm - AVG_PRE_ByRegion), 1, 'omitnan'); %change from average BL pre- and post-stim bins averages across regions
                ChgAvgBL_Post = mean((EventFt_Perm(:, NumPreBins+1:end) - AVG_PRE_ByRegion), 1, 'omitnan'); %change from average BL post-stim bins averages across regions
                
                EventFt_Perm_Shifted = EventFt_Perm;
                if Pre_Size_Diff > 0
                    EventFt_Perm_Shifted = [NaN(size(EventFt_Perm_Shifted,1), Pre_Size_Diff), EventFt_Perm_Shifted];
                    ChgAvgBL_PrePost = [NaN(1, Pre_Size_Diff), ChgAvgBL_PrePost];
                else
                    EventFt_Perm_Shifted = EventFt_Perm_Shifted(:,1-Pre_Size_Diff:end);
                    ChgAvgBL_PrePost = ChgAvgBL_PrePost(:,1-Pre_Size_Diff:end);
                end
                
                if Post_Size_Diff > 0
                    EventFt_Perm_Shifted = [EventFt_Perm_Shifted, NaN(size(EventFt_Perm_Shifted,1), Post_Size_Diff)];
                    ChgAvgBL_PrePost = [ChgAvgBL_PrePost, NaN(1, Post_Size_Diff)];
                    ChgAvgBL_Post = [ChgAvgBL_Post, NaN(1, Post_Size_Diff)];
                else
                    EventFt_Perm_Shifted = EventFt_Perm_Shifted(:,1:end+Post_Size_Diff);
                    ChgAvgBL_PrePost = ChgAvgBL_PrePost(:,1:end+Post_Size_Diff);
                    ChgAvgBL_Post = ChgAvgBL_Post(:,1:end+Post_Size_Diff);
                end
                
                XSecBinsPrePostPerm(TS,:) = mean(EventFt_Perm_Shifted, 1, 'omitnan');
                XSecBinsChangeAvgBaselinePerm(TS,:) = ChgAvgBL_Post;
                XSecBinsChangeAvgBLPrePostPerm(TS,:) = ChgAvgBL_PrePost;
                
%                 tempChange = mean((EventFt_Perm(:, NumPreBins+1:end) - mean(EventFt_Perm(:, 1:NumPreBins), 2, 'omitnan')), 1, 'omitnan');
%                 if Post_Size_Diff > 0
%                     tempChange = [tempChange, NaN(1, Post_Size_Diff)];
%                 else
%                     tempChange = tempChange(:,1:end+Post_Size_Diff);
%                 end
%                 XSecBinsChangeAvgBaselinePerm(TS,:) = tempChange;
            end
            clear ExcludeBaselineBins EVENT_IDX_PRE EVENT_IDX_POST NumPreBins...
                NumPostBins TotalBins EVENT_IDX_ALL EVENT_IDX_PERM EventFt_Perm...
                tempChange Post_Size_Diff Pre_Size_Diff EventFt_Perm_Shifted...
                AVG_PRE_ByRegion ChgAvgBL_PrePost ChgAvgBL_Post
        end
        TS = TS+1;
    end
    
    XSecBinsPrePost_PermMean(perm, :) = mean(XSecBinsPrePostPerm, 1, 'omitnan');
    XSecBinsChangeAvgBaseline_PermMean(perm,:) = mean(XSecBinsChangeAvgBaselinePerm, 1, 'omitnan');
    XSecBinsChangeAvgBLPrePost_PermMean(perm,:) = mean(XSecBinsChangeAvgBLPrePostPerm, 1, 'omitnan');
    clear IndPerm IndPerm_ChangeFromBaseline
end