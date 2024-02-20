function [FirstFrame_AboveBaseline] = FirstFrame_AboveBaselineMean(FluorescenceTrace, BaselinePeriod, z, STDorSEM, SearchAfterFrameX)
%FirstFrame_AboveBaselineMean finds the first frame greater than or equal
%to Z STD or SEM above the baseline mean
% Michelle Cahill 20211111 adapted from Peak3STD written by MC 20211005
%   FluorescenceTrace: 
%   BaselinePeriod: [a b] a: start frame of baseline period, b: end frame
%   of baselineperiod
%   z: the number of std or sem above the mean to set the threshold
%   STDorSEM: 'std' or 'sem' to determine how to set the threshold
%   SearchAfterFrameX: the frame number to begin looking for increases. Any
%   increases before this frame will be ignored

% 1. Calculate Z standard deviations or standard error of mean above baseline mean
BL = FluorescenceTrace(BaselinePeriod(1):BaselinePeriod(2));
BL_Mean = nanmean(BL);
if strcmp(STDorSEM, 'std')
    VARIATION = nanstd(BL);
elseif strcmp(STDorSEM, 'sem')
    VARIATION = nanstd(BL)/sqrt(length(BL)-1);
end
FrameThreshold = BL_Mean + (z*VARIATION);

% 2. Find frames above the threshold that occur after the cutoff
AboveThresholdFrames = find(FluorescenceTrace >= FrameThreshold);

% 3. Find the first frame above the threshold occurring after the manually
% designated frame (SearchAfterFrameX)
FramesAfterCutoff = AboveThresholdFrames(AboveThresholdFrames > SearchAfterFrameX);
FirstFrame_AboveBaseline = min(FramesAfterCutoff);

end
