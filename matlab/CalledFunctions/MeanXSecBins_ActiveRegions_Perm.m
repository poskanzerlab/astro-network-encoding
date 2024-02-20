function [XsecBin_MeanPerm] = MeanXSecBins_ActiveRegions_Perm(mydata, condition, CONDITION, NUM_REGIONS, TSERIES, Xsec, NormBy, perm_num,...
    NumNaNBins, Trace2Plot, AVERAGE_ACROSS, ExcludeBaseline, ExcludeEnd)
% MeanXSecBinsRegions_Perm performs a permutation test for n-permutations (perm_num) on active pixels per frame
% data. The function returns a nx1000 matrix, with each row being the mean
% across a t-series per frame/Xsec bin for one round of permutation.
% mean taken across active regions of all t-series, regions within the same
% t-series are shuffled together
% Check this matrix against the experimental mean to see how many times the
% permuted mean is greater than the experimental mean to obtain a p-value
% for each frame or x-sec bin
% Each t-series is shuffled independently of one another in
% each round of permutations
% mydata: structed array with AQuA results files aggregated
% condition: {'condition1', 'condition2'} the conditions from which to
% determine active cells
% CONDITION: the specific condition/fieldname of mydata to run the permutation test
% for (ex: 'RuBiGABA_WT' for mydata.(CONDITION))
% NUMREGIONS: number of active cells in that condition
% TSERIES: indices of t-series to average across
% XSec: The length of each bin (1, 5, 10, etc.) 1S is usually just 1 frame
% after rounding
% NormBy: a string identifying which normalizing method the active frames
% per pixel data should use (ex: '_NormBaselineActPx')
% perm_num: number of permutations to carry out
% NumNaNBins: number of dummy bins to put values in (match 2nd dimension of
% matrices in XsecBin_AllActPx used for experimental mean)

% Currently not used
% ExcludeBaseline: time in seconds to remove from the beginning of all
% recordings
% ExcludeEnd: time in seconds to remove from the end of all recordings
% Michelle Cahill 20220124 modified from MeanXSecBinsNoRegions_perm.m 20210518
XsecBin_MeanPerm = NaN(perm_num, NumNaNBins);
for perm = 1:perm_num
    if strcmp(AVERAGE_ACROSS, 'TSERIES')
        XsecBin_AllActPx = NaN(numel(TSERIES),2000);
    elseif strcmp(AVERAGE_ACROSS, 'REGIONS')
        XsecBin_AllActPx = NaN(NUM_REGIONS,2000);
    end
    
    Region_count = 0;
    %populate the matrix with % active px values for each tseries
    for tseries = TSERIES
        if ~isempty(mydata(tseries).(CONDITION)) && ~isempty(mydata(tseries).(CONDITION).res.ftsFilter.basic)
            
            if contains(Trace2Plot, 'dff')
                PerFieldActive = mydata(tseries).(CONDITION).(EBL).(Trace2Plot);
            elseif contains(Trace2Plot, 'dFF')
                % Finding cells with at lest 1 AQuA event within across all
                % conditions in that FOV
                CellsWithActivityCond1 = sum(mydata(tseries).(condition{1}).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                CellsWithActivityCond2 = sum(mydata(tseries).(condition{2}).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                CellsWithActivity = sum([CellsWithActivityCond1; CellsWithActivityCond2], 1);
                % Finding cells with at lest 1 AQuA event within this
                % condition
                %                 CellsWithActivity = sum(mydata(tseries).(CONDITION).res.ftsFilter.region.cell.memberIdx, 1, 'omitnan'); %Count the number of AQuA events in each cell
                ACTIVE_REGIONS = find(CellsWithActivity);
                if strcmp(AVERAGE_ACROSS, 'TSERIES')
                    PerFieldActive = mean(mydata(tseries).(CONDITION).ImageJByRegion.(Trace2Plot)(:,ACTIVE_REGIONS),2, 'omitnan')';
                elseif strcmp(AVERAGE_ACROSS, 'REGIONS')
                    PerFieldActive = mydata(tseries).(CONDITION).ImageJByRegion.(Trace2Plot)(:,ACTIVE_REGIONS)';
                end
                clear CellsWithActivity ACTIVE_REGIONS CellsWithActivityCond1 CellsWithActivityCond2
            else
                PerFieldActive = mydata(tseries).(CONDITION).(Trace2Plot);
            end
            
            SecPerFrame = mydata(tseries).(CONDITION).SecPerFrame;
            NumFrames_Xsec = round(Xsec/SecPerFrame); %The number of frames in x seconds, listed above as 5sec (Xsec = 5)
            FrameOfIncrease = mydata(tseries).(CONDITION).uncagingframe;
            %             EXCLUDEFRAMES_START = ceil(ExcludeBaseline/SecPerFrame);
            %             EXCLUDEFRAMES_END = ceil(ExcludeEnd/SecPerFrame);
            %             FrameOfIncrease = mydata(tseries).(CONDITION).uncagingframe-EXCLUDEFRAMES_START;
            
            %             PerFieldActive = PerFieldActive(:, 1+EXCLUDEFRAMES_START:end-EXCLUDEFRAMES_END);
            Shuffled_PerFieldActive = PerFieldActive(:, randperm(size(PerFieldActive, 2)));
            
            for reg = 1:size(Shuffled_PerFieldActive, 1)
                Region_count = Region_count + 1;
                FramesCroppedBeginning = mod(FrameOfIncrease-1, NumFrames_Xsec);
                if FramesCroppedBeginning == 0
                    PerFieldActive_Pre = Shuffled_PerFieldActive(reg, 1 : FrameOfIncrease-1);
                else
                    PerFieldActive_Pre = Shuffled_PerFieldActive(reg, FramesCroppedBeginning+1 : FrameOfIncrease-1);
                end
                
                FramesCroppedEnd = mod((length(PerFieldActive)- FrameOfIncrease), NumFrames_Xsec);
                PerFieldActive_Post = Shuffled_PerFieldActive(reg, FrameOfIncrease : end-(FramesCroppedEnd+1));
                
                if contains(Trace2Plot, 'ActivePx')
                    PerFieldActivePerXSec_Pre = sum(reshape(PerFieldActive_Pre, NumFrames_Xsec, []),1);
                    PerFieldActivePerXSec_Post = sum(reshape(PerFieldActive_Post, NumFrames_Xsec, []),1);
                else
                    PerFieldActivePerXSec_Pre = mean(reshape(PerFieldActive_Pre, NumFrames_Xsec, []),1);
                    PerFieldActivePerXSec_Post = mean(reshape(PerFieldActive_Post, NumFrames_Xsec, []),1);
                end
                
                XsecBin_AllActPx(Region_count, ((NumNaNBins/2)-length(PerFieldActivePerXSec_Pre)+1):(NumNaNBins/2)) = PerFieldActivePerXSec_Pre;
                XsecBin_AllActPx(Region_count, ((NumNaNBins/2)+1):((NumNaNBins/2)+length(PerFieldActivePerXSec_Post))) = PerFieldActivePerXSec_Post;
            end
            
            
            %             PerFieldActive = mydata(tseries).(CONDITION).(sprintf('ActivePxPerFrame%s', NormBy));
            %             Shuffled_PerFieldActive = PerFieldActive(randperm(size(PerFieldActive, 2)));
            % %             Shuffled_PerFieldActive = nan(size(PerFieldActive));
            % %             for r = 1:size(PerFieldActive,1)
            % %                 Shuffled_PerFieldActive(r,:) = PerFieldActive(r, randperm(size(PerFieldActive, 2)));
            % %             end
            %             FrameOfIncrease = mydata(tseries).(CONDITION).uncagingframe;
            %             SecPerFrame = mydata(tseries).(CONDITION).SecPerFrame;
            %             NumFrames_Xsec = round(Xsec/SecPerFrame); %The number of frames in x seconds, listed above as 5sec (Xsec = 5)
            %             FramesCroppedBeginning = mod(FrameOfIncrease-1, NumFrames_Xsec);
            %             FramesCroppedEnd = mod((length(PerFieldActive)- FrameOfIncrease + 1), NumFrames_Xsec);
            %
            % %             if FramesCroppedBeginning == 0
            % %                 PerFieldActive_Pre = Shuffled_PerFieldActive(1 : FrameOfIncrease-1);
            % %             else
            %                 PerFieldActive_Pre = Shuffled_PerFieldActive(FramesCroppedBeginning+1 : FrameOfIncrease-1);
            % %             end
            %             PerFieldActivePerXSec_Pre = sum(reshape(PerFieldActive_Pre, NumFrames_Xsec, []),1);
            %
            %             PerFieldActive_Post = Shuffled_PerFieldActive(FrameOfIncrease : end-(FramesCroppedEnd+1));
            %             PerFieldActivePerXSec_Post = sum(reshape(PerFieldActive_Post, NumFrames_Xsec, []),1);
            %
            %             tseries_idx = find(TSERIES == tseries);
            %             XsecBin_AllActPx(tseries_idx, ((NumNaNBins/2)-length(PerFieldActivePerXSec_Pre)+1):(NumNaNBins/2)) = PerFieldActivePerXSec_Pre;
            %             XsecBin_AllActPx(tseries_idx, ((NumNaNBins/2)+1):((NumNaNBins/2)+length(PerFieldActivePerXSec_Post))) = PerFieldActivePerXSec_Post;
            %
            clear PerFieldActive FrameOfIncrease SecPerFrame NumFrames_Xsec FramesCroppedBeginning FramesCroppedEnd tseries_idx
        end
    end
    XsecBin_MeanPerm(perm,:) = nanmean(XsecBin_AllActPx, 1);
    clear XsecBin_AllActPx tseries
end

