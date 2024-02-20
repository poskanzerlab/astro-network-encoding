%% For Fig. 1C top: plot event count per frame across an entire FOV
%
%  1. Load 'AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat' workspace
%  2. Run'Fig1_PreppingDataStruct.m':
%     Section A: experiment_type = 'BathApp'; method_AgonistEntry = 'MaxCurve'
%     Section B: Run as is
%  3. Run the following code (t-series 15 was used for the traces shown in
%     Fig. 1C top)
%
% Modified from first section of
% EventCountPerFrame_AndBinned_ForLongerTseries_20190904.m
%
% Michelle Cahill 20231213
%%
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\Fig1C\Top';

% Calculate the event count per frame for each t-series (raw and normalized
% by total number of events)
% Plot each t-series, aligning stim frame to 0 seconds
for cond = 1:length(condition)
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            clear EventCountPerFrame
            EventCountPerFrame = zeros(1, mydata(tseries).(condition{cond}).totalframes); %1 x # of frames
            num_evts = length(mydata(tseries).(condition{cond}).res.ftsFilter.loc.t0);
            for evt = 1:num_evts %for each event, find the frames in which the event took place and add them to the correct row & column in the EventCount matrix
                active_frames = unique(mydata(tseries).(condition{cond}).tLocs{1,evt}); %Pull out the active frames for that event
                EventCountPerFrame(1, min(active_frames):max(active_frames)) = EventCountPerFrame(1, min(active_frames):max(active_frames)) + 1; %add 1 to every frame of that cell that that event occurred
                clear active_frames
            end
            clear evt
            mydata(tseries).(condition{cond}).EventCountPerFrame = EventCountPerFrame;
            mydata(tseries).(condition{cond}).EventCountPerFrame_NormAllEvts = EventCountPerFrame ./ num_evts;
            clear num_evts EventCountPerFrame
            
            OnsetFrame = mydata(tseries).(condition{cond}).uncagingframe;
            SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
            
            % Shifting the Onset Frame to 0 and converting from frames to seconds
            XVALS = 1:mydata(tseries).(condition{cond}).totalframes;
            XVALS = (XVALS - OnsetFrame) * SecPerFrame;
            
            if strcmp(experiment_type, 'BathApp')
                DyeTrace = mydata(tseries).(condition{cond}).Alexa594.RawTraceXY(:,2);
            end
            
            traces2plot = {'EventCountPerFrame', 'EventCountPerFrame_NormAllEvts'};
            for t = 1:length(traces2plot)
                figure() %For each tseries plot the # of events per frame with the fluorescence trace for Alexa 594 beneath, for two wash in conditions
                if strcmp(experiment_type, 'BathApp')
                    subplot(2,1,1)
                end
                plot(XVALS, mydata(tseries).(condition{cond}).(traces2plot{t}))
                ylabel('count')
                xlabel('time (s)')
                title(sprintf('%s %s(%duM): %s t-series %d', experiment_type, condition{cond}, mydata(tseries).(condition{cond}).concentration, strrep(traces2plot{t}, '_', ' '), tseries))
                
                if strcmp(experiment_type, 'BathApp')
                    subplot(2,1,2)
                    hold on
                    plot(XVALS, DyeTrace, 'k') %Fluorescence traces from ImageJ
                    title('Alexa594 entering bath')
                end
                
                cd(save_dir)
                saveas(gcf, sprintf('%s%duM_%s_tseries%d.tif', condition{cond}, mydata(tseries).(condition{cond}).concentration, traces2plot{t}, tseries))
                saveas(gcf, sprintf('%s%duM_%s_tseries%d.svg', condition{cond}, mydata(tseries).(condition{cond}).concentration, traces2plot{t}, tseries))
                close
                cd(start_dir)
            end
            clear t traces2plot XVALS DyeTrace OnsetFrame SecPerFrame
        end
    end
end
clear tseries cond save_dir start_dir