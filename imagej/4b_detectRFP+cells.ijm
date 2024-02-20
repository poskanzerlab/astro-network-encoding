dir = getDirectory("Select a Folder");
fileList = getFileList(dir);
output_dir = dir + File.separator + "RFP+ area CSV" + File.separator;
output_dir2 = dir + File.separator + "RFP+ ROI" + File.separator;
File.makeDirectory(output_dir);
File.makeDirectory(output_dir2);

for (i = 0; i < lengthOf(fileList); i ++) {
	current_imagePath = dir + fileList[i];
	if (!File.isDirectory(current_imagePath)){
		open(current_imagePath);
		setOption("ScaleConversions", true);
		run("8-bit");
		run("Gaussian Blur...", "sigma=3");
		run("Auto Local Threshold", "method=Phansalkar radius=1000 parameter_1=0 parameter_2=0 white");
		run("Analyze Particles...", "size=100-Infinity circularity=0.00-0.60 clear add");
		if (roiManager("Count") > 0){
			roiManager("Measure");
			fileName = fileList[i];
			saveAs("Results", output_dir+ fileName + "-Results.csv"); //saves area of GFP +ve ROI to "GFP+ area"
			selectWindow("ROI Manager");
			roiManager("Save", output_dir2 + getTitle() + "-ROIset.zip"); //saves ROI of GFP +ve cells to "GFP+ ROI"
			selectWindow("ROI Manager");
			run("Close");
		}
	}
	run("Close All");
}