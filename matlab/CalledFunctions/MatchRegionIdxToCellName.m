function [OrigRegionIdxStruct,UpdatedRegionIdxStruct] = MatchRegionIdxToCellName(RegionIdxStruct)
% MatchRegionIdxToCellName takes the AQuA res file output
% res.ftsFilter.region.cell and reorders the cell indices of each variable
% to match any manual re-naming of regions in cell.name
% The original .cell struct will be maintained in OrigRegionIdxStruct
% The updated .cell struct will have regions re-ordered to match cell.name
% in OrigRegionIdxStruct.
% To track that regions have been reordered:
% ...cell.IdxMatchManualRenaming = 'Yes'
% ...cell.name will be re-ordered
% Michelle Cahill 20230807
OrigRegionIdxStruct = RegionIdxStruct;
UpdatedRegionIdxStruct = struct();

OrigNames = RegionIdxStruct.name;
OrigIdx = cell2mat(cellfun(@(x) str2double(x), OrigNames, 'UniformOutput', false)); %Convert cell name characters into an array of doubles
StructFieldsCellIdxRows = {'mask','center','border','centerBorderAvgDist','incluLmk'};
for f = 1:length(StructFieldsCellIdxRows)
    temp_OrigVar = RegionIdxStruct.(StructFieldsCellIdxRows{f});
    if ~isempty(temp_OrigVar)
        UpdatedRegionIdxStruct.(StructFieldsCellIdxRows{f}) = temp_OrigVar(OrigIdx, :);
    else
        UpdatedRegionIdxStruct.(StructFieldsCellIdxRows{f}) = temp_OrigVar;
    end
    clear temp_OrigVar
end
clear f StructFieldsCellIdxRows

StructFieldsCellIdxCol = {'memberIdx', 'dist2border', 'dist2borderNorm', 'name'};
for f = 1:length(StructFieldsCellIdxCol)
    temp_OrigVar = RegionIdxStruct.(StructFieldsCellIdxCol{f});
    UpdatedRegionIdxStruct.(StructFieldsCellIdxCol{f}) = temp_OrigVar(:, OrigIdx);
    clear temp_OrigVar
end
clear f StructFieldsCellIdxCol

UpdatedRegionIdxStruct.IdxMatchManualRenaming = 'Yes';
end
