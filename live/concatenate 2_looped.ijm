directory = "/Users/yezhu/Desktop/Lab file_Ye/ECad-NCad reporter line/221118 PS_bi_13/rawFile/";
folder1 = directory + "Day0-1_445CAAX_514ECAD_561NCAD_40X_Cycle/";
folder2 = directory + "Day1-2_445CAAX_514ECAD_561NCAD_40X_Cycle/";
outputDir = "/Users/yezhu/Desktop/Lab file_Ye/ECad-NCad reporter line/221118 PS_bi_13/MATLAB_live/rawData/";
list1 = getFileList(folder1);
list2 = getFileList(folder2); 
Array.show(list1);
Array.show(list2);

for (i = 1; i < lengthOf(list1); i++) {
	if(endsWith(list1[i], ".oir")) { 
		
		open(folder1 + list1[i]);
		open(folder2 + list2[i]);
		selectWindow(list1[i]);
		selectWindow(list2[i]);
		
		title1 = "445CAAX_514ECAD_561NCAD" + String.pad(i+1,2);
		
		commandString = "  title=[" + title1 + "] keep  image1=" + list1[i] + "  image2=" + list2[i];
		print(commandString);
//concatenate
		run("Concatenate...", commandString);
		saveAs("TIFF", outputDir + title1 + ".tif");
		run("Close All");
	}
}