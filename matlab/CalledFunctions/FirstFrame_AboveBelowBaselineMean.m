function [FirstFrame_AboveBelowBaseline, IncDecNoChange] = FirstFrame_AboveBelowBaselineMean(FluorescenceTrace, BaselinePeriod, z, STDorSEM, SearchAfterFrameX)
%FirstFrame_AboveBaselineMean finds the first frame greater than or equal
%to Z STD or SEM above the baseline mean
% Michelle Cahill 202111202 adapted from FirstFrame_AboveBaselineMean
% written by MC 20211111
%   FluorescenceTrace:
%   BaselinePeriod: [a b] a: start frame of baseline period, b: end frame
%   of baselineperiod
%   z: the number of std or sem above the mean to set the threshold
%   STDorSEM: 'std' or 'sem' to determine how to set the threshold
%   SearchAfterFrameX: the frame number to begin looking for increases or decreases. Any
%   changes before this frame will be ignored

% 1. Calculate Z standard deviations or standard error of mean above baseline mean
BL = FluorescenceTrace(BaselinePeriod(1):BaselinePeriod(2));
BL_Mean = nanmean(BL);
if strcmp(STDorSEM, 'std')
    VARIATION = nanstd(BL);
elseif strcmp(STDorSEM, 'sem')
    VARIATION = nanstd(BL)/sqrt(length(BL)-1);
end
FrameThresholdAbove = BL_Mean + (z*VARIATION);
FrameThresholdBelow = BL_Mean - (z*VARIATION);

% 2. Find frames above and below the threshold that occur after the cutoff
AboveThresholdFrames = find(FluorescenceTrace >= FrameThresholdAbove);
BelowThresholdFrames = find(FluorescenceTrace <= FrameThresholdBelow);

% 3. Find the first frame above or below the threshold occurring after the manually
% designated frame (SearchAfterFrameX)
%    If more frames are above threshold IncDecNoChange == 1
%    If more frames are below threshold IncDecNoChange == -1
%    If no frames are above or below threshold IncDecNoChange == 0
AboveFramesAfterCutoff = AboveThresholdFrames(AboveThresholdFrames > SearchAfterFrameX);
BelowFramesAfterCutoff = BelowThresholdFrames(BelowThresholdFrames > SearchAfterFrameX);

NumFramesAbove = numel(AboveFramesAfterCutoff);
NumFramesBelow = numel(BelowFramesAfterCutoff);

if NumFramesAbove > NumFramesBelow
    IncDecNoChange = 1;
    FirstFrame_AboveBelowBaseline = min(AboveFramesAfterCutoff);
elseif NumFramesBelow > NumFramesAbove
    IncDecNoChange = -1;
    FirstFrame_AboveBelowBaseline = min(BelowFramesAfterCutoff);
elseif NumFramesBelow == 0 && NumFramesAbove == 0
    IncDecNoChange = 0;
    FirstFrame_AboveBelowBaseline = NaN;
else
    IncDecNoChange = NaN;
    FirstFrame_AboveBelowBaseline = NaN;
end
end
