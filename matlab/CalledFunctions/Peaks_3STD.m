function [ peak_frames, peak_dFFMax, peak_latency ] = Peaks_3STD( dFF_trace, uncaging_frame, SecPerFrame, num_std, peakthresh, MinPeakWidth)
%Peaks_3STD finds peaks that are >= 3std above the baseline mean, and
%records the frames and the amplitude of these peaks, split into those
%peaks that occur before and after a certain frame
%   dFF_trace: Nx1 matrix of a normalized trace
%   uncaging_frame: the first frame of some manipulation (to split pre- and post-stim frames)
%   num_std: How many standard deviations above the mean should the peak be
%   to be considered?
%   peakthresh: how many standard deviations above neighboring samples
%   should a peak be?
%   MinPeakWidth: How wide should a peak be to be detected?
%
%   peak_frames: 1x2 cell array, with column 1: frames of peaks pre-stim
%   column 2: frames of peaks post-stim
%   peak_dFFMax: 1x2 cell array, with column 1: dFF_Max of peaks pre-stim
%   column 2: dFF_Max of peaks post-stim
%   peak_latency: 1x2 cell array, with column 1: latency (in sec) of peak
%   since start of trace, column 2: latency (in sec) of peak
%   since uncaging/stim
% modified on 20220126 to add 'MinPeakDistance': 10, 'Threshold':
% X*baseline_std & 20220128 'MinPeakWidth': user_defined

baseline_mean = nanmean(dFF_trace(1:(uncaging_frame-1),1));
baseline_std = nanstd(dFF_trace(1:(uncaging_frame-1),1));
[pks,frame] = findpeaks(dFF_trace, 'MinPeakHeight', (baseline_mean+(num_std*baseline_std)),...
    'MinPeakDistance', 10, 'Threshold', peakthresh*baseline_std, 'MinPeakWidth', MinPeakWidth);
findpeaks(dFF_trace, 'MinPeakHeight', (baseline_mean+(num_std*baseline_std)),...
    'MinPeakDistance', 10, 'Threshold', peakthresh*baseline_std, 'MinPeakWidth', MinPeakWidth); %to display which peaks are considered
peak_frames_pre = frame(frame < uncaging_frame);
peak_dFFMax_pre = pks(frame < uncaging_frame);
peak_latency_pre = peak_frames_pre * SecPerFrame;
peak_frames_post= frame(frame >= uncaging_frame);
peak_dFFMax_post = pks(frame >= uncaging_frame);
peak_latency_post = (peak_frames_post - uncaging_frame) * SecPerFrame; 
peak_frames = {peak_frames_pre, peak_frames_post};
peak_dFFMax = {peak_dFFMax_pre, peak_dFFMax_post};
peak_latency = {peak_latency_pre, peak_latency_post};

end

