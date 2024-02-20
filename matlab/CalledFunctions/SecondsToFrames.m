function [Num_frames] = SecondsToFrames(Num_seconds,SecPerFrame)
%SecondsToFrames converts seconds to frames
%   Num_seconds: the number of seconds
%   SecPerFrame: the conversion factor - the number of seconds per frame
%   Num_frames: the number of frames (rounded up to the nearest whole
%   number)
Num_frames = ceil(Num_seconds/SecPerFrame);
end

