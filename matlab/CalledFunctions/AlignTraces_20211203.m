function [ALIGNED_TRACES, AlignOnset2Frame] = AlignTraces_20211203(TRACES,OnsetFrame)
%AlignTraces_20211203 takes a matrix of fluorescence traces and the frame of
%some stimulus for each sample. Individual traces are aligned to the frame
%of the stimulus
% Michelle Cahill 20211203
%  TRACES: matrix of fluorescence traces, rows are individual samples and
%  columns are frames
%  OnsetFrame: for every sample, what is the frame to align to
%  Output:
%  ALIGNED_TRACES: rows are individual samples and columns are frames. The
%  onset frame or time of the stimulus is the middle column
%  (AlignOnset2Frame)
%  AlignOnset2Frame: The column in ALIGNED_TRACES taht indicates stimulus
%  onset
[num_traces, num_frames] = size(TRACES);
SizeMatrix = round(num_frames, -1) * 2; %the matrix will be 2X the number of frames(rounded to the nearest 10)
AlignOnset2Frame = SizeMatrix/2 + 1;

ALIGNED_TRACES = NaN(num_traces, SizeMatrix);
for t = 1:num_traces
    trace = TRACES(t, :);
    OF = OnsetFrame(t);
    START = AlignOnset2Frame - OF + 1;
    END = START + num_frames - 1;
    ALIGNED_TRACES(t,START:END) = trace;
    clear trace OF START END
end

