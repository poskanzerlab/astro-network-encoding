function [ interpulse_duration_ms, num_pulses, num_bursts, pulseduration_ms, uncaging_start_ms, uncaging_end_ms, uncaging_duration_ms, uncaging_frame_start, uncaging_frame_end, max_voltage ] = square_pulse_detection( r, frameRate )
%square_pulse_detection detects square pulses from a voltage recording
%(column 1: time points, column 2: voltage) and the frame rate of the video
%associated with the voltage recording
%   The function detects a train of square pulses. This version is for
%   recordings in which there is only one train of pulses
% Output arguments:
% interpulse_duration_ms: time between pulses (in ms if
% voltage recording was done on the order of ms) [within a burst; between
% bursts]
% num_pulses: number of pulses within a burst
% nu_bursts: number of bursts
% pulseduration_ms: duration of a single pulse (ms, assuming all
% pulses are the same duration)
% uncaging_start_ms: time of the start of the first square pulse (ms) for
% each burst
% uncaging_end_ms: time of the end of the last square pulse (ms) for each
% burst
% square pulse (ms) uncaging_frame_start: time of the start of the first
% square pulse (frame) for each burst 
% uncaging_frame_end: time of the end of the last square pulse (frame) for
% each burst

if isa(r, 'struct')
uncaging_recording = table2array(readtable(r.name));
elseif isa(r, 'double') %added 20240111 to work with recordings already in array format
    uncaging_recording = r;
end

max_voltage = max(uncaging_recording(:,2));
uncagingtime_ms = uncaging_recording((uncaging_recording(:,2) > (max_voltage-0.02)),1); %Find times (in ms) where the voltage is < 0.02 from the max voltage
continuous_pulse = diff(uncagingtime_ms); %calculate if there were multiple pulses (values in this case will jump above 0.1 or 1)
[~, pulseend_idx] = findpeaks(continuous_pulse, 'MinPeakHeight', 1); %find the indices of breaks in the pulse of 'continuous_pulse', this corresponds to the end of each pulse in uncagingtime_ms
if ~isempty(pulseend_idx) %if multiple pulses were detected
    interpulse_duration_ms = unique(continuous_pulse(pulseend_idx)); %find the interpulse_duration (assumes this value is consistent between pulses)
    if interpulse_duration_ms(1) > 150 %if there is greater than 150ms between pulses, consider the pulses to be part of separate bursts
        interpulse_duration_ms = [nan; interpulse_duration_ms];
        num_pulses = 1;
        pulsestart_idx = [1 ; pulseend_idx+1]; %find the indices of the start of each pulse for uncagingtime_ms
        pulseend_idx = [pulseend_idx ; length(continuous_pulse)]; %include end of continuous pulse for end of the last pulse
        num_bursts = numel(pulseend_idx); % Count the number of bursts detected in the voltage recording
        pulseduration_ms = ceil((uncagingtime_ms(pulseend_idx) - uncagingtime_ms(pulsestart_idx))/10)*10; %calculate the duration of each pulse (rounded to the nearest ms)
        pulseduration_ms = unique(pulseduration_ms);
        uncaging_start_ms = uncagingtime_ms(pulsestart_idx);
        uncaging_end_ms = uncagingtime_ms(pulseend_idx);
        uncaging_duration_ms = uncaging_end_ms - uncaging_start_ms;
        uncaging_frame_start = ceil((uncaging_start_ms/1000) / frameRate)+1;
        uncaging_frame_end = ceil((uncaging_end_ms/1000) / frameRate)+1;
    else %otherwise consider multiple pulses within a single burst
%         IDX_burstbreak = find(continuous_pulse == interpulse_duration_ms(1)); %Find the indices separating bursts
        IDX_burstbreak = find(continuous_pulse == interpulse_duration_ms(2)); %Find the indices separating bursts index of 2 was used prior to 20200908
        num_pulses = numel(find(continuous_pulse(1:IDX_burstbreak(1)) == interpulse_duration_ms(1))) + 1; %Count the # of pulses taking place in one burst (asssumes all bursts have same num of pulses)
        num_bursts = numel(IDX_burstbreak) + 1;
        pulsestart_idx = [1 ; pulseend_idx+1]; %find the indices of the start of each pulse for uncagingtime_ms
        pulseend_idx = [pulseend_idx ; length(continuous_pulse)]; %include end of continuous pulse for end of the last pulse
        pulseduration_ms = ceil((uncagingtime_ms(pulseend_idx) - uncagingtime_ms(pulsestart_idx))/10)*10; %calculate the duration of each pulse (rounded to the nearest ms)
        pulseduration_ms = unique(pulseduration_ms);
        
        [~,burstend_idx] = findpeaks(continuous_pulse, 'MinPeakHeight', 1000); %find indices of breaks between bursts, this corresponds to the end of the burst in uncagingtime_ms
        burststart_idx = [1; burstend_idx+1];
        burstend_idx = [burstend_idx; length(continuous_pulse)];
        uncaging_start_ms = uncagingtime_ms(burststart_idx);
        uncaging_end_ms = uncagingtime_ms(burstend_idx);
        uncaging_duration_ms = uncaging_end_ms - uncaging_start_ms;
        uncaging_frame_start = ceil((uncaging_start_ms/1000) / frameRate)+1;
        uncaging_frame_end = ceil((uncaging_end_ms/1000) / frameRate)+1;
    end
else
    interpulse_duration_ms = nan;
    num_pulses = 1;
    num_bursts = 1;
    pulseduration_ms = uncagingtime_ms(end) -uncagingtime_ms(1);
    uncaging_start_ms = uncagingtime_ms(1);
    uncaging_end_ms = uncagingtime_ms(end);
    uncaging_duration_ms = uncaging_end_ms - uncaging_start_ms;
    uncaging_frame_start = ceil((uncaging_start_ms/1000) / frameRate)+1;
    uncaging_frame_end = ceil((uncaging_end_ms/1000) / frameRate)+1;
end
end

