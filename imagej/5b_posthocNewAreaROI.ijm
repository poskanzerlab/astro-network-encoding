dir_images = getDirectory("Open folder of GFP/RFP+ images");
fileList = getFileList(dir_images);
dir_ROIs = getDirectory("Select ROIs corresponding to folder of images");
roiList = getFileList(dir_ROIs);
output_dir = dir_ROIs + File.separator + "ROI Area" + File.separator;
File.makeDirectory(output_dir);

for (i = 0; i < lengthOf(fileList); i ++) {
	current_imagePath = dir_images + fileList[i];
	current_roiSet = dir_ROIs + roiList[i];
	if (!File.isDirectory(current_imagePath)){
		open(current_imagePath);
		roiManager("Open", current_roiSet);
		roiManager("Show All");
		roiManager("Measure");
		fileName = fileList[i];
		saveAs("Results", output_dir+ fileName + "-Results.csv");
		selectWindow("ROI Manager");
		run("Close");
		selectWindow("Results");
		run("Close");
	}
	run("Close All");
}
