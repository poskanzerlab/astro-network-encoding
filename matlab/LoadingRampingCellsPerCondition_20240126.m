%% Section 1. Single round 2P uncaging data WT, Receptor Antagonist and Laser
% Uncaging Control
start_dir = pwd;
ramping_dir = 'G:\Michelle_Event_Detection\EventDet_2PGABAGluUncSameCell_localnet_Cx43fl\PostEvtDet\WTCX43CBX_FullDatasets\ramping_2022-02-25_14-44-10';
% file_names = {'ramping_2021-05-24_002_Control.csv', 'ramping_2021-05-24_002_GABA-CGP.csv', 'ramping_2021-05-24_002_Glu-LY.csv', 'ramping_2021-05-24_002_RuBi.csv'};
% file_names = {'ramping_2021-05-24_002_CBX.csv', 'ramping_2021-05-24_002_Cx43.csv', 'ramping_2021-05-24_002_RuBi.csv'};
file_names = {'ramping_2022-02-25_14-44-10_GABA.csv', 'ramping_2022-02-25_14-44-10_Glu.csv',...
    'ramping_2022-02-25_14-44-10_GABA-CGP.csv','ramping_2022-02-25_14-44-10_Glu-LY.csv',...
    'ramping_2022-02-25_14-44-10_Control-Laser.csv',...
    };
cd(ramping_dir)

%Load .csv ramping cell information into the workspace
ramping_cells = struct();
for f = 1:length(file_names)
    temp_cond = extractBefore(extractAfter(file_names{f}, 'ramping_2022-02-25_14-44-10_'), '.csv');
    temp_cond = strrep(temp_cond, '-', '_');
    ramping_cells.(temp_cond) = readtable(file_names{f});
    clear temp_cond
end
clear f file_names ramping_dir
cd(start_dir)
%% Section 2. Multiple rounds 2P glutamate uncaging (70 A.U.)
start_dir = pwd;
ramping_dir = 'G:\Michelle_Event_Detection\EventDet_2PUncagingRepeatTrials_70AU\PostEventDet\RampingCell_IDs';
file_names = {'ramping_2023-08-07_11-25-43_Glu-repeat-70.csv',...
    };
cd(ramping_dir)

%Load .csv ramping cell information into the workspace
ramping_cells = struct();
for f = 1:length(file_names)
    temp_cond = extractBefore(extractAfter(file_names{f}, 'ramping_2023-08-07_11-25-43_'), 'repeat-70.csv');
    temp_cond = strrep(temp_cond, '-', '_');
    ramping_cells.(temp_cond) = readtable(file_names{f});
    clear temp_cond
end
clear f file_names ramping_dir
cd(start_dir)

clear start_dir
%% Section 3. Multiple rounds 2P glutamate uncaging (25 A.U.)
start_dir = pwd;
ramping_dir = 'G:\Michelle_Event_Detection\EventDet_2PUncagingRepeatTrials\PostEventDet\RampingCell_IDs';
file_names = {'ramping_2023-08-09_14-36-43_Glu-repeat.csv',...
    };
cd(ramping_dir)

%Load .csv ramping cell information into the workspace
ramping_cells = struct();
for f = 1:length(file_names)
    temp_cond = extractBefore(extractAfter(file_names{f}, 'ramping_2023-08-09_14-36-43_'), 'repeat.csv');
    temp_cond = strrep(temp_cond, '-', '_');
    ramping_cells.(temp_cond) = readtable(file_names{f});
    clear temp_cond
end
clear f file_names ramping_dir
cd(start_dir)

clear start_dir