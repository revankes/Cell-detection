/*	StackCellDetection_3Ch_V4.IJM
	***********************
	Author: 		Titia
	Date Created: 		January 16, 2015

	1) detect cells in green channel, but can be addapted 
 	2) measure intensity and other params in the Green / Red / Blue channel

 	note: only takes ROIs with circularity of 0.5 and larger
 	
 	**********************************************

		Variable initiation

	**********************************************

*/
//	Arrays
var list 	= newArray();

//	Strings and numbers
var dir = "";						//	directory
var outDir = "";					//	output directory
var index = 0;						//	Index of folder (well)

var minSize = 80, maxSize = 2000;			//	min and max size of nuclei for segmentation (in pixels)
var minCirc = 0;			
var enlarge = -1;
/*
 	**********************************************

		Signal Measurements

	**********************************************
*/

macro "[M] Measure Nuclear Foldings"
{
	setup();
	setBatchMode(true);
	print("Analysis");
	print("************************************************************");
	index = 0;
	for(i=0; i<list.length; i++)
	{
			
		path = dir + list[i];
		if(endsWith(path,'.tif'))
		{
			run("Collect Garbage");
			open(path);
			id = getImageID;
			title = getTitle;
			print(title);
			prefix = substring(title,0,indexOf(title,".tif"));
			print("Image",i,":",prefix);			
						
			selectImage(id);
			print(title);
			run("Duplicate...", "title=GREENchannel.tif duplicate channels=1");
			GREENid2 = getImageID;			
			
			selectImage(id);
			run("Duplicate...", "title=REDchannel.tif duplicate channels=1");
			GREENid = getImageID;

			selectImage(GREENid);
			roiManager("reset");
			run("Despeckle","stack");
			run("Gaussian Blur...", "sigma=1 stack");
			setAutoThreshold("Huang"+" dark stack");//Huang IJ_Isodata
			run("Convert to Mask", "stack");
			run("Fill Holes", "stack");
			run("Watershed", "stack");
			run("Analyze Particles...", "size="+80+"-"+2000+" circularity="+0.5+"-1.00 show=Nothing exclude add stack");
			n = roiManager("count");
			print(n, "Nuclei");
			roiManager("deselect");
			run("Clear Results");
			run("Set Measurements...", "area mean standard min centroid center bounding shape median stack redirect=None decimal=4");
			selectImage(GREENid2);
			roiManager("Measure");
			selectWindow("Results");
			saveAs("Measurements",dir+prefix+"_Results_Green.txt");
			roiManager("Save", dir+prefix+"_Roiset_Green.zip");
			
			//measure values in other channels			
			selectImage(id);
			ns = nSlices;
			run("Duplicate...", "title=REDchannel.lsm duplicate channels=2 slices=1-"+ns);
			REDid = getImageID;
			titleRED = getTitle;
			print(titleRED);
			roiManager("deselect");
			run("Clear Results");
			roiManager("Measure");
			saveAs("Measurements",dir+prefix+"_Results_Red.txt");
			selectImage(REDid);
			close;

			selectImage(id);
			ns = nSlices;
			run("Duplicate...", "title=BLUEchannel.lsm duplicate channels=3 slices=1-"+ns);
			BLUEid = getImageID;
			titleBLUE = getTitle;
			print(titleBLUE);
			roiManager("deselect");
			run("Clear Results");
			roiManager("Measure");
			saveAs("Measurements",dir+prefix+"_Results_Blue.txt");
			selectImage(BLUEid);
			close;
		
			selectImage(id);
			close;
			selectImage(GREENid2);
			close;

		}
	}
	print("************************************************************");
	
}

/*	
 	**********************************************

		Functions

	**********************************************
*/


function getMoment()
{
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+second;
     return TimeString;
}

function setup()
{
	print("\\Clear");
	run("Close All");
	run("Clear Results");
	roiManager("reset");
	run("Collect Garbage");
	run("Colors...", "foreground=white background=black selection=yellow");	
	setOption("BlackBackground", false);
	run("Set Measurements...", "area mean shape stack redirect=None decimal=4");
	dir = getDirectory("");
	list = getFileList(dir);
	isWin=indexOf(getInfo("os.name"),"Windows")>=0;
	TimeString = getMoment();
	print(TimeString);
	print("************************************************************");
}


function segment(GREENid2)
{	
	selectImage(GREENid2);
	roiManager("reset");
	run("Despeckle");
	run("Gaussian Blur...", "sigma=1");
	setAutoThreshold("Li"+" dark");
	run("Convert to Mask");
	run("Fill Holes");
	run("Watershed");
	run("Analyze Particles...", "size="+80+"-"+2000+" pixel circularity="+0+"-1.00 show=Nothing exclude clear add");
	n = roiManager("count");
	print(n, "Nuclei");
	return n;
}