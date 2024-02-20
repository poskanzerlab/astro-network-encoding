function [ YVal_LAremoved ] = RemovingLightArtifacts( YVals, start_frames, end_frames)
%% RemovingLightArtifacts removes the light artifact in a fluorescence trace, given that the light artifact is the brightest thing in the t-series
% YVals = the raw fluorescence trace
% start_frames = # of frames to remove before the brightest frame
% end_frames = # of frames to remove after the brightest frame
% Michelle Cahill 20210816

%%
EndBaseline = round(length(YVals)/2) - 20;
[~, LA_max_frame] = max(YVals(1:EndBaseline));
LA_AppStart = LA_max_frame - start_frames;
LA_AppEnd = LA_max_frame + end_frames;

YVal_LAremoved = YVals;
YVal_LAremoved(LA_AppStart:LA_AppEnd) = YVals(LA_AppStart-1);

end
