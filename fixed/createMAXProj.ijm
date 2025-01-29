run("Close All");
directory = "/Users/yezhu/Desktop/Lab file/Data/CadSwitch FIXED/240729 EcadKO 2/";

rawDIR= directory + "rawData/Plate b_405DAPI_488_561_640_Cycle_01/"
outputDIR = directory + "MAXProj_original/Plate a";
File.makeDirectory(outputDIR);

list = getFileList(rawDIR) 
Array.show(list);

for (i = 0; i < lengthOf(list); i++) {
	if(endsWith(list[i], ".oir")) { 
		
		open(rawDIR + list[i]);
		rawName=File.getNameWithoutExtension(list[i]);
		run("Z Project...", "projection=[Max Intensity] all");
		//saveAs("TIFF", outputDIR + rawName + "_MAXProj.tif");
		
		Stack.setDisplayMode("composite");
		Stack.setActiveChannels("1111");
		saveAs("JPEG", outputDIR + rawName + "_MAXProj1111.jpg");
		run("Close All");
	}
	
}
