# astro-network-encoding - MATLAB scripts

Accompanying the paper "Network-level encoding of local neurotransmitters in cortical astrocytes" (DOI: TODO)

## Prerequisites

Code written and run using MATLAB 2018b (Michelle Cahill).

The following toolboxes were installed:
* Simulink version 9.2
* 5G Toolbox version 1.0
* Bioinformatics Toolbox version 4.11
* Communications Toolbox version 7.0
* Curve Fitting Toolbox version 3.5.8
* DSP System Toolbox version 9.7
* Filter Design HDL Coder Version 3.1.4
* Fixed-Point Designer Version 6.2
* Parallel Computing Toolbox Version 6.13
* Signal Processing Toolbox Version 8.1
* Statistics and Machine Learning Toolbox Version 11.4
* Symbolic Math Toolbox Version 

## Instructions

For each figure panel created in MATLAB, open the indicated starting script and follow the directions in the opening comments for the specific figure panel. These directions will include:

1. Which dataset(s) to load (also indicated in the table below)
2. Which scripts need to be run prior to running the current script (such as `Fig1_PreppingDataStruct.m` or `Fig2_3_PreppingDataStruct.md`)
3. Which sections should be run and which variables should be chosen in the current script 

### Scripts for figure panels

| Figure Panel | Data to load | Starting MATLAB script |
| ------------ | ------------ | ---------------------- |
| Fig. 1c (top) | `AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat` | `Fig1CTop` |
| Fig. 1c (bottom) | '`AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat`' | `Fig1CBottom` |
| Fig. 1d | `AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat` | `Fig1D_ExtDataFig1f` |
| Fig. 1e | `AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat` | `Fig1E` |
| Fig. 1f–h | `AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat` | `Fig1F_H_ExtendedDataFig1g` |
| Fig. 2d | `AQuA_GluSnFR_2PGluUncaging.mat` | `Fig2d_ExtDataFig7a` |
| Fig. 2h | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig2h_3c_f` |
| Fig. 2i–j | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig2i_j_ExtDataFig2a_c_3d` |
| Fig. 2l–o | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig2l_o_ExtDataFig2f` |
| Fig. 3b | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig3b` |
| Fig. 3c–f | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig2h_3c_f` |
| Fig. 3h | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig3h_ExtDataFig3i` |
| Ext. Data Fig. 1a | `GSE161398_FPKM_MasterTable_Development.xlsx` | `ExtendedDataFig1a_ExtDataFig3a` |
| Ext. Data Fig. 1b–c | `AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat` | `ExtendedDataFig1b_c` |
| Ext. Data Fig. 1d–e | `AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat` | `ExtendedDataFig1d_e` |
| Ext. Data Fig. 1f | `AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat` | `Fig1D_ExtDataFig1f` |
| Ext. Data Fig. 1g | `AQuA_CytoGCaMP_ReceptorAgonistBathApp_Fig1.mat` | `Fig1F_H_ExtendedDataFig1g` |
| Ext. Data Fig. 1j | `FIJI_pinkFlamindo_ReceptorAgonistBathApp_ExtFig1.mat` | `ExtendedDataFig1j` |
| Ext. Data Fig. 1k | `FIJI_pinkFlamindo_ReceptorAgonistBathApp_ExtFig1.mat` and `FIJI_AQuA_CytoGCaMP_BathAppBacLY_ExtFig1k.mat` | `ExtendedDataFig1k` |
| Ext. Data Fig. 2a | `AQuA_CytoGCaMP_2PUncaging_WTRecAntaLaserCtrl_ExtDataFig2_3.mat` | `Fig2i_j_ExtDataFig2a_c_3d` |
| Ext. Data Fig. 2b | `AQuA_CytoGCaMP_2PUncaging_WTRecAntaLaserCtrl_ExtDataFig2_3.mat` | `ExtendedDataFig2b_ExtDataFig3g_h_ExtDataFig7d` |
| Ext. Data Fig. 2c | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig2i_j_ExtDataFig2a_c_3d` |
| Ext. Data Fig. 2d | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `ExtendedDataFig2d_ExtDataFig5f_g` |
| Ext. Data Fig. 2e | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `ExtendedDataFig2e__ExtDataFig3e_f` |
| Ext. Data Fig. 2f | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig2l_o_ExtDataFig2f` |
| Ext. Data Fig. 3a | `GSE161398_FPKM_MasterTable_Development.xlsx` | `ExtendedDataFig1a_ExtDataFig3a` |
| Ext. Data Fig. 3d | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig2i_j_ExtDataFig2a_c_3d` |
| Ext. Data Fig. 3e–f | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `ExtendedDataFig2e__ExtDataFig3e_f` |
| Ext. Data Fig. 3g | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `ExtendedDataFig2b_ExtDataFig3g_h_ExtDataFig7d` |
| Ext. Data Fig. 3h | `AQuA_CytoGCaMP_2PUncaging_WTRecAntaLaserCtrl_ExtDataFig2_3.mat` | `ExtendedDataFig2b_ExtDataFig3g_h_ExtDataFig7d` |
| Ext. Data Fig. 3i | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `Fig3h_ExtDataFig3i` |
| Ext. Data Fig. 5a | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `ExtendedDataFig5a` |
| Ext. Data Fig. 5f–g | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `ExtendedDataFig2d_ExtDataFig5f_g` |
| Ext. Data Fig. 5h–i | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `ExtendedDataFig5h_i` |
| Ext. Data Fig. 5j | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` | `ExtendedDataFig5j` |
| Ext. Data Fig. 7a | `AQuA_GluSnFR_2PGluUncaging.mat` | `Fig2d_ExtDataFig7a` |
| Ext. Data Fig. 7b | `AQuA_CytoGCaMP_2Puncaging_WTCx43CBX_Fig2_3.mat` and `AQuA_CytoGCaMP_2PMultiRoundGluUncaging_70AU_ExtDataFig7.mat` | `ExtendedDataFig7b` |
| Ext. Data Fig. 7d (top) | `AQuA_CytoGCaMP_2PMultiRoundGluUncaging_70AU_ExtDataFig7.mat` | `ExtendedDataFig2b_ExtDataFig3g_h_ExtDataFig7d` |
| Ext. Data Fig. 7d (bottom) | `AQuA_CytoGCaMP_2PMultiRoundGluUncagingCtrl_25AU_ExtDataFig7.mat` | `ExtendedDataFig2b_ExtDataFig3g_h_ExtDataFig7d` |
| Extended Data Fig. 7e | `AQuA_CytoGCaMP_2PMultiRoundGluUncaging_70AU_ExtDataFig7.mat` | `ExtendedDataFig7e` |

## Notes

* p-values calculated using permutation testing will differ slightly each time permutation testing is re-run, given that the shuffling of items for each round of permutation is random.

* The MATLAB workspaces available (listed in ‘Data to Load’ column) contain structured arrays with AQuA res files from multiple recordings (for details see the DryadReadMe file). To load new AQuA res files into a structured array that can run through any of the code listed above, run `AggregateAQuAResFilesIntoDataStruct_mydata_20240131.m`.