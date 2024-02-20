function [LABELS] = MultiLineLabels(LabelMatrix)
%MultiLineLabels takes a cell array of labels (columns are different
%labels, rows will be different lines of the labels) and converts it into
%individual labels

BaseString = '%s';
AdditionalString = '\\newline%s';
AddStringRep = repmat(AdditionalString, 1, size(LabelMatrix,1)-1);
EndString = '\n';
NEW_STRING = strcat(BaseString, AddStringRep, EndString);

LABELS = strtrim(sprintf(NEW_STRING, LabelMatrix{:})); %allows for multiline xlabels

% strtrim(sprintf('%s\\newline%s\\newline%s\n', DataLabels_PVAL{:})); %allows for multiline xlabels
end

