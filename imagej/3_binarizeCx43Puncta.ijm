dir = getDirectory("Select Cx43 Image Folder");
fileList = getFileList(dir);
output_dir = dir + File.separator + "Binary puncta" + File.separator;
File.makeDirectory(output_dir);

for (i = 0; i < lengthOf(fileList); i++) {
	current_imagePath = dir + fileList[i];
	if(!File.isDirectory(current_imagePath)){
		open(current_imagePath);
		run("Enhance Contrast", "saturated=0.35");
		//synquant parameters optimized from Yu lab 1/6/21 VT for Cx43 IHC image
		run("SynQuantVid ", "z-score=0 min=5 max=200 min_0=0.50 max_0=4 post-synapse=[current_imagePath] pre-synapse=Null way=Null dendrite=Null extended=0 z=1 zscore=8");
		newImage("Untitled", "8-bit black", 2048, 2048, 1);
		roiManager("Show All without labels");
		roiManager("Set Fill Color", "yellow");
		run("Flatten");
		run("8-bit");
		setAutoThreshold("Default dark");
		//run("Threshold...");
		//setThreshold(129, 255);
		setOption("BlackBackground", true);
		run("Convert to Mask");
		fileName = fileList[i];
		selectWindow("Untitled-1");
		rename(fileName);
		saveAs("tiff", output_dir+fileName);
		selectWindow("Synapse detection results");
		run("Close");
	}
	run("Close All");
}
