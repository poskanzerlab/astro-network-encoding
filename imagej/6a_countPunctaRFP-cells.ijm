dir_puncta = getDirectory("Select Binary Puncta");
fileList = getFileList(dir_puncta);
dir_ROIs = getDirectory("Select ROIs");
roiList = getFileList(dir_ROIs);
output_dir = dir_puncta + File.separator + "RFP- Puncta Count" + File.separator;
File.makeDirectory(output_dir);

for (i = 0; i < lengthOf(fileList); i ++) {
	current_imagePath = dir_puncta + fileList[i];
	current_roiSet = dir_ROIs + roiList[i];
	if (!File.isDirectory(current_imagePath)){
		open(current_imagePath);
		roiManager("Open", current_roiSet);
		roiManager("Show All");
		numROIs = roiManager("Count");
		for (j = 0; j < numROIs ; j++) {
			roiManager("Select", j);
			run("Set Measurements...", "area redirect=None decimal=3");
			run("Analyze Particles...", "display clear");
			fileName = fileList[i];
			roiName = Roi.getName();
			saveAs("Results", output_dir + fileName + "-ROI # " + roiName + "-GFP.csv");
			}
		roiManager("Deselect");
		roiManager("Delete");
		if (isOpen("Results")){
		selectWindow("Results");
		run("Close");
		}
	}
	run("Close All");
}
