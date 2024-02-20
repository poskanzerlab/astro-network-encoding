function [dFFTraces] = dFF_lowerXPercent_20220127(RawTraces, SecPerFrame, MovingWindow, LowerXPercent)
% dFF_lower calculates the dF/F from raw traces, taking the lower X% of
% frames in the moving window to calculated the mean fluorescence
% Michelle Cahill 20220127 (adapted from Nicole DelRosso and from dFF_lower_20211202.m)
%   RawTraces: m x n matrix of raw fluorescence. m = number of frames; n =
%   number of individual traces. Df/f will be calculated for each
%   individual trace (column)
%   SecPerFrame: temporal resolution
%   MovingWindow: # in seconds for the width of the moving window


%here, write in your movie parameters
time_res = SecPerFrame;
options.nSecs.value = MovingWindow;
%pathname = '~/Documents/MATLAB';%wherever your tsmovavg.m file is stored
%cd(pathname)
nrecordings = size(RawTraces,1);
len = size(RawTraces,2);


%F0 = mean of the lowest X% of n previous frames
% We'd like to get
% (Fi - F0)/F0

nFrames=ceil(options.nSecs.value/time_res);
if mod(nFrames,2)==0
    nFrames=nFrames+1;  %make sure we take an odd window size for the median filter.
end

mid=ceil(nFrames/2);
dFFTraces = [];

%df/f lower 
for ii = 1:len %loop through each column
    y=RawTraces(:,ii); %pull out current trace for ease of indexing later on.
    foo=mean(y(1:nFrames));
    for j=1:nFrames %set fo = mean of first nFrames for the first nFrames of X.
        dFFTraces(j,ii)=(y(j) - foo)/ foo;
    end
    for j=nFrames+1:nrecordings
        fo=y(j-nFrames:j-1); %look at previous nFrames
        fo=sort(fo);
        XPERCENT_FRAMES = floor(length(fo) * (LowerXPercent/100));
        fo=fo(1:XPERCENT_FRAMES); % keep the lower X% of previous nFrames
        fo=mean(fo); %take their mean
        dFFTraces(j,ii) = (y(j) - fo)./ fo; %subtract and divide raw data by mean
    end
end
end

