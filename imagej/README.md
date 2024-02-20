# astro-network-encoding - ImageJ macros

Connexin 43 Puncta Count Pipeline (Extended Data Fig. 3c).

Accompanying the paper "Network-level encoding of local neurotransmitters in cortical astrocytes" (DOI: TODO)

## Prerequisites

Macros written and run using FIJI/ImageJ (version 1.53c) (Vincent Tse)

Users will need to install the following Fiji/ImageJ plugin for puncta detection: https://github.com/yu-lab-vt/SynQuant.

## Instructions

### Notes

* ImageJ macros in Extended Data Fig. 3C Code folder are numerically labeled in the order that they should be executed. Each script builds off the other (e.g. input->output->input), so it is important to follow the below steps. All new folders created will be stored in the initial input folder. Be sure to use the corresponding folders/files when two are required.
* This pipeline is intended for .TIFF z-stack images with 3-color channels acquired from confocal microscopy. Two channels should mark astrocyte morphology and the other channel should contain a signal with individual puncta spread across the field-of-view. These images should be stored in the same folder on local machine for the following scripts to work. All new output folders will be contained in this folder to enable flexible organization.

### Pipeline steps

1. `1_splitChannels.ijm` takes a folder containing multi-channel .TIFF z-stack images and splits them into single-channel stacks. Returns a new folder titled 'Split channels'. Run this once. Each single-channel stack will indicate channels as a pre-fix in the file name (e.g. "C1").
2. `2_stackToImages.ijm` takes a folder containing single-channel z-stack images and unstack them into individual z-plane images, which is the output of '1_splitChannels'. Returns a new folder titled 'Z-planes'. User needs to then manually organize z-plane images into new folders for each fluorescent signal (ie. separate folders for Cx43, GFP and RFP). 
3. `3_binarizeCx43Puncta.ijm` takes the folder containing connexin 43 z-plane images from '2_stackToImages' and binarizes them via SynQuant Fiji/ImageJ plugin. Returns a new folder titled 'Binary puncta' containing binarized connexin 43 z-plane images.
4. `4a_detectGFPcells.ijm` and `4b_detectRFP+cells.ijm` perform the same function but take the folders containing either GFP+ or RFP+ z-plane images from `2_stackToImages` to create ROIs for # GFP+ and RFP+ cells. For each, returns two new folders ending with "ROI" which contains the ROIs (.zip files opened via ROI Manager) for each z-plane and "area CSV" which contains area measurements for each ROI in each z-plane. The files in each folder are ordered and can be visually checked by assessing the numerical suffix in the file name, which indicate z-plane order.
5. Once GFP+ ROIs have been detected, `4c_detectRFP-cells.ijm` needs to be run to generate ROIs that are actually RFP- as some GFP+ cells are also RFP+. This takes two inputs: a folder containing the RFP+ images, and a folder containing the GFP+ ROIs. This requires manual changes to the ROIs in the ROI Mananger, in which user should delete the GFP+ ROIs from the ROI Manager that directly coincide with RFP+ cells. Returns a new folder titled 'RFP- ROI'. The 'GFP+ ROI' folder generated from Step 4 can be ignored moving forward. These new RFP- ROIs will need new area measurements, but this will be done after any necessary post-hoc adjustments are done (Step 6).
6. `5a_posthocAdjustROI.ijm` can be used to sequentially and manually adjust RFP- and RFP+ ROIs. Requires two input folders: a folder of either the GFP+ (i.e. RFP-) or RFP+ images, and the corresponding folder containing ROIs. Returns a new folder titled 'post-hoc ROIs' which will contain new ROIs in .zip file format for each z-plane image. Repeat for the other channel.
7. `5b_posthocNewAreaROI.ijm` can then be run to calculate the new area for all adjusted ROIs. Requires two input folders: a folder of either the GFP+ (i.e. RFP-) or RFP+ images, and the corresponding post-hoc ROI folder generated in Step 6. Returns a folder titled "ROI Area" with .csv files for each z-plane image.
8. `6a_countPuncta-cells.ijm` and `6b_countPuncta+cells.ijm` are used to count the number of connexin 43 puncta in either RFP- or RFP+ cells. Requires two folders: a folder containing the connexin 43 binary image, and either the RFP- or RFP+ ROIs in .zip file format. Returns a folder ending with "Puncta Count". Repeat for the other channel. The detection pipeline is finished and ready for data processing and analysis.
    * **Note**: The output from 6a is a .csv file for each Z-plane. Combine these into a single .csv file per sample/mouse and calculate the puncta count/astrocyte area. The normalized puncta count per mouse in RFP+ vs RFP- astrocytes can now be plotted and compared.