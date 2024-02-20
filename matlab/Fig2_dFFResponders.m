%% Identifying responding uncaging cells for Fig 2
%  List the indices of tseries in each condition in which the uncaging cell was a "responder"
%  "responder" cells are any cells with post-stim dFF values 3std above baseline
%  From 'dFF_Responders_MC20220908.m'
%
%  *note that this is a different designation for responders than identified in
%  Fig. 4 and related extended data figures (which is based on event
%  frequency increases)*
%
%  1. Load 'AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat' workspace
%  2. Run the following sections of 'Fig2_3_PreppingDataStruct.m':
%     Section A: experiment_type = 'uncaging'; MultUncagingReps = 0;
%     Section B: data_sets = {'WT', 'Cx43', 'CBX'};
%     Section C & D as is
%  3. Run section A of 'Fig2h_3c_f.m':
%     ExcludeRampingCells = 1; RampingCellThresh = 'pval_1'; ExcludeBaseline = 60; 
%     WE = 'UncagingEvents';
%  4. Run the following code using the following parameters:
%     WE = 'UncagingEvents'; ExcludeBaseline = 60; time_window = [];
%
%  Calls function 'SecondsToFrames.m' if time_window isn't empty
%  Michelle Cahill 20240112
%%
condition = fieldnames(mydata);
WE = 'UncagingEvents'; %'UncagingEvents', 'NeighboringEvents' 

ExcludeBaseline = 60; %number of seconds at the beginning of baseline to exclude, 0 to include all baseline period
EBL = sprintf('Exclude%dSBL', ExcludeBaseline);

time_window = []; %time (in seconds) to look for a response post-stim, [] for all post-stim time 

dFFResponders = struct();

for cond = 1:length(condition)
    temp_responders = [];
    for tseries = 1:length(mydata)
        if ~isempty(mydata(tseries).(condition{cond}))
            UncagingFrame = mydata(tseries).(condition{cond}).uncagingframe;
            SecPerFrame = mydata(tseries).(condition{cond}).SecPerFrame;
            if ~isempty(time_window)
                Num_Frames = SecondsToFrames(time_window, SecPerFrame);
                FramesAboveThresh = mydata(tseries).(condition{cond}).(EBL).(WE).dffMatFilter2_Thresh3std(UncagingFrame:UncagingFrame+Num_Frames);
            else
                FramesAboveThresh = mydata(tseries).(condition{cond}).(EBL).(WE).dffMatFilter2_Thresh3std(UncagingFrame:end);
            end
            if sum(FramesAboveThresh) > 0
                temp_responders = [temp_responders, tseries];
            end
            clear FramesAboveThresh UncagingFrame Num_Frames SecPerFrame
        end
    end
    temptime = sprintf('poststim%ds', time_window);
    dFFResponders.(temptime).(WE).(condition{cond}) = temp_responders;
    clear temp_responders temptime
end

clear tseries cond condition pulse_num WE ExcludeBaseline EBL time_window