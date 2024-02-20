%% Load voltage recording of uncaging laser Pockels cell for each recording included in mydata
% Assumes the filename listed in each recording of my data links to a 
% folder thatincldues a voltage recording .csv file

% Michelle Cahill 20240111
%%
start_dir = pwd;
condition = fieldnames(mydata);
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
       if ~isempty(mydata(tseries).(condition{cond}))
           FN = mydata(tseries).(condition{cond}).file; %AQuA res folder for the recording
           cd(FN) %go to the folder for that recording
           VRName = dir('*VoltageRecording_001.csv'); %load in the electrical recording for the uncaging laser if this file exists
           VR = table2array(readtable(VRName.name));
           mydata(tseries).(condition{cond}).UncLaserVoltageRecording = VR;
           % Check that data loaded in correct format
%            [interpulse_duration_ms, num_pulses, pulseduration_ms, uncaging_start_ms, uncaging_end_ms, uncaging_duration_ms, uncaging_frame_start, uncaging_frame_end, max_voltage ] = square_pulse_detection_20200408( mydata(tseries).(condition{cond}).UncLaserVoltageRecording, mydata(tseries).(condition{cond}).res.opts.frameRate );
           cd(start_dir)
           clear FN VR VRName
       end
    end
end
clear cond condition tseries

%%
clear interpulse_duration_ms num_pulses pulseduration_ms uncaging_start_ms uncaging_end_ms uncaging_duration_ms uncaging_frame_start uncaging_frame_end max_voltage start_dir