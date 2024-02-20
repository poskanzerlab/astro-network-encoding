function [ interpulse_duration_ms, num_pulses, pulseduration_ms, uncaging_start_ms, uncaging_end_ms, uncaging_duration_ms, uncaging_frame_start, uncaging_frame_end, max_voltage ] = square_pulse_detection( r, frameRate )
%square_pulse_detection detects square pulses from a voltage recording
%(column 1: time points, column 2: voltage) and the frame rate of the video
%associated with the voltage recording
%   The function detects a train of square pulses. This version is for
%   recordings in which there is only one train of pulses
% Output arguments: 
% interpulse_duration_ms: time between pulses (in ms if
% voltage recording was done on the order of ms) 
% num_pulses: number of pulses 
% pulseduration_ms: duration of a single pulse (ms, assuming all
% pulses are the same duration) 
% uncaging_start_ms: time of the start of the first square pulse (ms) 
% uncaging_end_ms: time of the end of the last
% square pulse (ms) uncaging_frame_start: time of the start of the first
% square pulse (frame) uncaging_frame_end: time of the end of the last
% square pulse (frame)

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
                pulsestart_idx = [1 ; pulseend_idx+1]; %find the indices of the start of each pulse for uncagingtime_ms
                pulseend_idx = [pulseend_idx ; length(continuous_pulse)]; %include end of continuous pulse for end of the last pulse
                num_pulses = numel(pulseend_idx); % Count the number of pulses detected in the voltage recording
                pulseduration_ms = round(uncagingtime_ms(pulseend_idx) - uncagingtime_ms(pulsestart_idx)); %calculate the duration of each pulse (rounded to the nearest ms)
                pulseduration_ms = unique(pulseduration_ms);
            else
                interpulse_duration_ms = nan;
                num_pulses = 1;
                pulseduration_ms = uncagingtime_ms(end) -uncagingtime_ms(1);
            end
            
            uncaging_start_ms = uncagingtime_ms(1);
            uncaging_end_ms = uncagingtime_ms(end);
            uncaging_duration_ms = uncaging_end_ms - uncaging_start_ms;
            
            if size(r, 1) > 350000
                %6/12/2023 the following is used to calculate uncaging frame
                %start and end (used for all uncaging recordings longer
                %than 5 min)
                uncaging_frame_start = floor((uncaging_start_ms/1000) / frameRate);
                uncaging_frame_end = floor((uncaging_end_ms/1000) / frameRate);
            else
                % Prior to 6/12/203 the following was used to calculated
                % uncaging frame start and end (used for all 5 min
                % recordings)
                uncaging_frame_start = ceil((uncaging_start_ms/1000) / frameRate)+1;
                uncaging_frame_end = ceil((uncaging_end_ms/1000) / frameRate)+1;
            end

end

