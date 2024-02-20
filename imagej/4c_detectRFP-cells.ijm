dir_images = getDirectory("Open folder of RFP+ images");
fileList = getFileList(dir_images);
dir_ROIs = getDirectory("Select RFP- ROIs");
roiList = getFileList(dir_ROIs);
output_dir = dir_images + File.separator + "RFP- ROI" + File.separator; //output folder will be where GFP z-plane images are
File.makeDirectory(output_dir);

for (i = 0; i < lengthOf(fileList); i ++) {
	current_imagePath = dir_images + fileList[i];
	current_roiSet = dir_ROIs + roiList[i];
	if (!File.isDirectory(current_imagePath)){
		open(current_imagePath);
		roiManager("Open", current_roiSet);
		roiManager("Show All");
		waitForUser("User Manual Selection", "Adjust ROIs. Press space to continue");
		selectWindow("ROI Manager");
		roiManager("Save", output_dir + getTitle() + "-post-hoc_ROIset.zip");
		selectWindow("ROI Manager");
		run("Close");
		}
	run("Close All");
}
