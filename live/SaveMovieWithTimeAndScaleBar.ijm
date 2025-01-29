directory = "/Users/yezhu/Desktop/Lab file_Ye/ECad-NCad reporter line/221118 PS_bi_13/MATLAB_live/rawData/";
outputDIR = "/Users/yezhu/Desktop/Lab file_Ye/ECad-NCad reporter line/221118 PS_bi_13/video/";
list = getFileList(directory) 
Array.show(list);

for (i = 9; i < 10; i++) {
	if(endsWith(list[i], ".tif")) { 
		
		open(directory + list[i]);
		rawName=File.getNameWithoutExtension(list[i]);
		
		run("Z Project...", "projection=[Max Intensity] all");
		Stack.setChannel(1);
		setMinAndMax(29, 4095);//HSD PS11,12
		//setMinAndMax(50, 1300);//SD PS8
		//setMinAndMax(50, 3000);//SD PS9
		//setMinAndMax(50, 1300);//SD PS10
		Stack.setChannel(2);
		setMinAndMax(56, 4095);//HSD PS11,12
		//setMinAndMax(30, 288);//SD PS8
		//setMinAndMax(30, 288);//SD PS9
		Stack.setChannel(3);
		setMinAndMax(73, 800);//HSD PS11,12
		//setMinAndMax(73, 600);//SD PS8
		//setMinAndMax(73, 720);//SD PS9
		run("Label...", "format=00:00 starting=60 interval=90 x=10 y=50 font=40 text=[] range=1-48");
		run("Scale Bar...", "width=20 height=5 thickness=10 font=20 color=White background=None location=[Lower Right] horizontal bold overlay");
		
		run("Split Channels");
		saveAs("Gif", outputDIR + rawName + "_MAX-C3.gif");
        close();
        saveAs("Gif", outputDIR + rawName + "_MAX-C2.gif");
        close();
        saveAs("Gif", outputDIR + rawName + "_MAX-C1.gif");
		run("Close All");
	}	
}
