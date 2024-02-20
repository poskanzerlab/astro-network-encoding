%% Extended Data Fig. 1a and 3a: Ribosomal-mRNA expression in visual cortex astrocytes from the Farhy-Tselnicker et. al publicly available dataset
% (NCBI Gene Expression Omnibus, GSE161398)

% 1. Load 'GSE161398_FPKM_MasterTable_Development.xlsx' as a table
% 2. Run the following code

% Extended Data Fig. 1a: Ribosomal-mRNA expression of GABA and glutamate
% GPCRs in visual cortex astrocytes at ages p14 and p28
%  Run Section A
%  Section B: Run separately for p14 and p28 data
%             GENE_NAME = [GABAB_names, mGluR_names];
%             AGE = 'P14' and 'P28';
%             protein = 'GPCRs';

% Extended Data Fig. 3a: Ribosomal-mRNA expression of connexin isoforms in
% visual cortex astrocytes at age p28.
%  Run Section A
%  Section B: 
%             GENE_NAME = [Cx_names];
%             AGE = 'P28';
%             protein = 'Connexins';

% modified from 'RNAseq_GeneExpression_MC20221101.m
% Michelle Cahill 20240110
%% A. Create lists of gene names to include in plots
% Extended Data Fig. 1a
GABAB_names = {'Gabbr1', 'Gabbr2'}; %GABA B receptor subunits gene names
mGluR_names = {'Grm1', 'Grm2', 'Grm3', 'Grm4', 'Grm5', 'Grm6', 'Grm7', 'Grm8'}; %mGluR subunit gene names

% Extended Data Fig. 3a
Cx_names = {'Gja1', 'Gja3', 'Gja4', 'Gja5', 'Gja6', 'Gja8', 'Gja10',...
    'Gjb1','Gjb2', 'Gjb3', 'Gjb4', 'Gjb5', 'Gjb6', 'Gjb7',...
    'Gjc1', 'Gjc2', 'Gjc3', 'Gjd2', 'Gjd3', 'Gjd4', 'Gje1'}; % connexin isoform gene names
%% B. Create box plots for the chosen age group and gene names. 
% Calculating the ratio of FPKM for gene of interest / FPKM for GFAP to 
% normalize for potential difference in the sequencing depth of replicates 
start_dir = pwd;
save_dir = 'G:\DifferentialEncodingOfNTs\NatureRevision1\GeneratedGraphs\ExtDataFig3a';

RNAseq_FullTable = GSE161398FPKMMasterTableDevelopment;

GENE_NAME = [Cx_names]; %[GABAB_names, mGluR_names]; %from the list of gene names created in Section A
AGE = 'P28'; %'P14' or 'P28'
protein = 'Connexins';% 'GPCRs' or 'Connexins'

temp_VariableNames = RNAseq_FullTable.Properties.VariableNames;
temp_SelectVariableNames = temp_VariableNames(contains(temp_VariableNames, AGE));
SAMPLES = temp_SelectVariableNames(~contains(temp_SelectVariableNames, 'IN')); %values given as FPKM
INPUT = temp_SelectVariableNames(contains(temp_SelectVariableNames, 'IN'));
clear temp_VariableNames temp_SelectVariableNames

% Pull the data from the genes and age of interest
GENES = find(sum(RNAseq_FullTable.Gene == GENE_NAME, 2)); %find the rows for the genes of interest
GFAP = find(RNAseq_FullTable.Gene == 'Gfap');
ColNames = RNAseq_FullTable.Properties.VariableNames;
ColIdx = find(ismember(ColNames, ['Gene', SAMPLES, INPUT]));
clear ColNames 

FPKM_val_table = RNAseq_FullTable([GENES;GFAP], ColIdx); %GFAP expression is the final column
temp = table2cell(FPKM_val_table);
FPKM_val_table = cell2table(temp', 'RowNames', FPKM_val_table.Properties.VariableNames,'VariableNames', FPKM_val_table.Gene);

FPKM_val = str2double(table2array(FPKM_val_table(2:end,1:end-1))); %take away the column of GFAP expression
GFAP_expression_BySample = str2double(table2array(FPKM_val_table(2:end,end))); %isolate the column of GFAP expression
FPKM_val_log = log2(FPKM_val);
relative_FPKM_val = FPKM_val(1:end-1,:) ./ FPKM_val(end,:);
relative_FPKM_val_log = log2(relative_FPKM_val);
FPKM_val_NormGFAP = FPKM_val ./ GFAP_expression_BySample;
clear temp

% Relative Astrocyte Expression (to GFAP expression for each sample)
figure()
boxplot(FPKM_val_NormGFAP(1:end-1,:), 'Labels', FPKM_val_table.Properties.VariableNames(1:end-1))
ylabel('FPKM gene of interest / FPKM GFAP')
xtickangle(90)
title(sprintf('Visual Cortex (%s): %s astrocyte expression', AGE, strrep(protein, '_', ' ')))
cd(save_dir)
saveas(gcf, sprintf('VisualCortex%s_%sExpression_RelativeAstrocytesToGFAPPerSample.tif', AGE, protein))
saveas(gcf, sprintf('VisualCortex%s_%sExpression_RelativeAstrocytesToGFAPPerSample.svg', AGE, protein))
close
cd(start_dir)
%%
clear start_dir save_dir GENE_NAME SAMPLES INPUT AGE protein GENES ColIdx FPKM_val FPKM_val_log...
relative_FPKM_val relative_FPKM_val_log FPKM_val_table RNAseq_FullTable FPKM_val_NormGFAP GFAP_expression_BySample GFAP