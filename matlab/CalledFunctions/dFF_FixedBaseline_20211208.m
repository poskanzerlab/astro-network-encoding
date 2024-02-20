function [dFFTraces] = dFF_FixedBaseline_20211208(RawTraces, BaselinePeriod)
% Calculating dF/F using a fixed baseline window on the raw traces
% Michelle Cahill 20211208 
%   RawTraces: m x n matrix of raw fluorescence. m = number of frames; n =
%   number of individual traces. Df/f will be calculated for each
%   individual trace (column)
%   BaselinePeriod: [x, y] x:starting baseline frame, y:ending
%   baseline frame

numFrames = size(RawTraces,1);
regions = size(RawTraces,2);

%F0 = mean baseline period frames
% We'd like to get
% (Fi - F0)/F0

dFFTraces = NaN(numFrames, regions);

%df/f fixed baseline 
for reg = 1:regions %loop through each column
    trace = RawTraces(:,reg); %pull out current trace for ease of indexing later on.
    BaselineMean = mean(trace(BaselinePeriod(1):BaselinePeriod(2)));
    dFFTraces(:,reg) = (trace - BaselineMean) ./ BaselineMean;
end
end
