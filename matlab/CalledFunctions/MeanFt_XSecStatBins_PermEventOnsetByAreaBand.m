function [XSecBinsPrePost_PermMean,XSecBinsChangeAvgBaseline_PermMean] = MeanFt_XSecStatBins_PermEventOnsetByAreaBand(mydata, CONDITION, TSERIES,...
    time_window, ExcludeBaseline, perm_num, param_name, WhichEvents, NBinsPre, NBinsPost,...
    AB, lmk_idx, bin)
% MeanFt_XSecStatBins_PermEventOnsetByAreaBand modified from
% MeanFt_XSecStatBins_PermEventOnsetByRegion performs a permutation test
% for n-permutations (perm_num) on the chosen feature/parameter.
% Currently written to shuffle events placed in x-sec time bins for the given area band
% in 'ShollLikeAnalysis_LocAtEvtStart_MultCond_XSecBins_MC20220322.m' and give a
% permuted mean per time bin for a chosen feature across chosen t-series.
% For each tseries, the average of the chosen feature is calculated in

% OUTPUTS:
% The function will return 2 matrices:
% XSecBinsPrePost_PermMean: rows will be the mean across t-series for one round of
% permutation, columns will be x-second bins both pre- and post-stim. N
% rows = perm_num
% XSecBinsChangeAvgBaseline_PermMean: rows will be the mean across t-series
% for one round of permutation, columns will be average change of post-stim
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
% WhichEvents: 'UncagingEvents', 'NeighboringEvents', 'AllEvents'
% If ramping cells are to be excluded, that's to be done when originally
% placing event indices into area bands
% NBinsPre: number of pre-stim x-sec bins to include
% NBinsPost: number of post-stim x-sec bins to include
% AB: sprintf('AreaBands%dum', area_band);
% lmk_idx: number of the landmark defined in AQuA for which events are being measured from for area band
% bin: number of the area band to shuffle

% Michelle Cahill 20220323

XSecBinsPrePost_PermMean = NaN(perm_num, NBinsPre+NBinsPost);
XSecBinsChangeAvgBaseline_PermMean = NaN(perm_num, NBinsPost);
TB = sprintf('StaticTimeBins%ds', time_window);

for perm = 1:perm_num
    XSecBinsPrePostPerm = NaN(length(TSERIES), NBinsPre+NBinsPost); %average value per X-sec bins per t-series for this round of permutation
    XSecBinsChangeAvgBaselinePerm = NaN(length(TSERIES), NBinsPost); %average change from baseline per X-sec post-stim bins per t-series for this round of permutation
    TS = 1;
    for tseries = TSERIES %For each t-series randomly assign each event to an Xsec bin for each region
        if ~isempty(mydata(tseries).(CONDITION))
            ExcludeBaselineBins = ExcludeBaseline / time_window;
            EVENT_IDX_PRE = mydata(tseries).(CONDITION).(TB).(AB)(lmk_idx).(WhichEvents).evt_idx_by_band_pre(bin,ExcludeBaselineBins+1:end);
            EVENT_IDX_POST = mydata(tseries).(CONDITION).(TB).(AB)(lmk_idx).(WhichEvents).evt_idx_by_band_post(bin,:);
            NumPreBins = size(EVENT_IDX_PRE, 2);
            NumPostBins = size(EVENT_IDX_POST, 2);
            TotalBins = NumPreBins + NumPostBins;
            EVENT_IDX_ALL_TEMP = [EVENT_IDX_PRE, EVENT_IDX_POST];
            EVENT_IDX_ALL = cellfun(@(x) x', EVENT_IDX_ALL_TEMP, 'UniformOutput', 0);
            
            EVENT_IDX_PERM = cell(size(EVENT_IDX_ALL));
            for region = 1:size(EVENT_IDX_PERM, 1) %There will only be 1 row, as all events from one area band are reported together
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
            end
            
            XSecBinsPrePostPerm(TS,:) = mean(EventFt_Perm, 1, 'omitnan');
            XSecBinsChangeAvgBaselinePerm(TS,:) = mean((EventFt_Perm(:, NumPreBins+1:end) - mean(EventFt_Perm(:, 1:NumPreBins), 2, 'omitnan')), 1, 'omitnan');
        end
        clear ExcludeBaselineBins EVENT_IDX_PRE EVENT_IDX_POST NumPreBins...
            NumPostBins TotalBins EVENT_IDX_ALL EVENT_IDX_PERM EventFt_Perm
        TS = TS+1;
    end
    XSecBinsPrePost_PermMean(perm, :) = mean(XSecBinsPrePostPerm, 1, 'omitnan');
    XSecBinsChangeAvgBaseline_PermMean(perm,:) = mean(XSecBinsChangeAvgBaselinePerm, 1, 'omitnan');
    clear IndPerm IndPerm_ChangeFromBaseline
end