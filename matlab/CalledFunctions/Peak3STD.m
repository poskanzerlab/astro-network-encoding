function [First_PostStim_PeakFrame] = Peak3STD(FluorescenceTrace,StimFrame, BaselinePeriod)
%Peak3STD finds the frame of the 1st peak greater than or equal to 3STD
%above the baseline mean following a stimulation
% Michelle Cahill 20211005
%   FluorescenceTrace: 
%   StimFrame: a number n. Find the first peak at or after n
%   BaselinePeriod: [a b] a: start frame of baseline period, b: end frame
%   of baselineperiod

% 1. Calculate 3 standard deviations above baseline mean
BL = FluorescenceTrace(BaselinePeriod(1):BaselinePeriod(2));
BL_Mean = nanmean(BL);
BL_STD = nanstd(BL);
PeakThreshold = BL_Mean + (3*BL_STD);

% 2. Find peaks that occur post-stim
[PEAKS, LOCS] = findpeaks(FluorescenceTrace);
POSTSTIM_LOCS = LOCS(LOCS >= StimFrame);
POSTSTIM_PEAKS = PEAKS(LOCS >= StimFrame);

% 3. Find the frame of the first peak greater than or equal to 3STD above
% the baseline mean
% X = POSTSTIM_PEAKS(POSTSTIM_PEAKS >= PeakThreshold);
temp_loc = POSTSTIM_LOCS(POSTSTIM_PEAKS >= PeakThreshold);
if ~isempty(temp_loc)
    First_PostStim_PeakFrame = min(temp_loc);
else
    First_PostStim_PeakFrame = nan;
end
end

