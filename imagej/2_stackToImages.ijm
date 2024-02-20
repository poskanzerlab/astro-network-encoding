dir = getDirectory("Select a Folder");
fileList = getFileList(dir);
output_dir = dir + File.separator + "Z-planes" + File.separator;
File.makeDirectory(output_dir);

for (i = 0; i < lengthOf(fileList); i++) {
	current_imagePath = dir+fileList[i];
	if(!File.isDirectory(current_imagePath)){
		open(current_imagePath);
		run("Stack to Images");

		ch_nbr = nImages;
		for (c = 1; c <= ch_nbr; c++){
			//fileName = fileList[i];
			selectImage(c);
			//currentImage_name = fileName + getTitle();
			currentImage_name = getTitle();
			saveAs("tiff", output_dir+currentImage_name);
		}
		run("Close All");
	}
}
