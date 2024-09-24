The three folders contain DICOM stacks that could be output by using
	functions laid out in iplfe_register3d.txt
	I could not get this text file to copy over correctly to the 
	Scanco software, but if you run the 'iplfe' command (without 
	quotes) in the Scanco software terminal, you should be able 
	to run the commands on two scans of your choice. This code
	was built to register two scans of the same bone, imaged
	at different times. 
	
The MATLAB code will chew through these DICOM stacks and ask for 
	two key pieces of information: the AP axis and medial side. 
	For the AP axis, please choose the anterior side first.
 	Each subsequent scan will default to the same AP axis as
	the previous.
	
The scan 1, 2, and difference figures can be navigated by slice #

The TV, BV, BMC, BMD variables are tables that report total volume,
	bone volume, bone mineral content, and bone mineral density in the 
	anterior, medial, lateral, and posterior quadrants, respectively.
	
Voxel resolution should be changed in the code if necessary, 
and the calibration slope and intercept should be machine-specific.

To setup a batch run, change the subjects initialization line to 
include all of the subjects you want to process.
E.g. subjects = ["MARA14-","MARA15-","MARA16-"];

If you would like to repeat the same axes for any two scans,
list them together with the first stack being the reference.

PLEASE NOTE: Each time you run this code a new output.xlsx will be
written to the code's home directory. This WILL be overwritten if 
you run the code again.

Update 09/24/24: Now requires a single center dicom nested in a folder
called {sub} {location} scan1_raw
This is required only from scan1 because scan2 is later registered to 
scan1. We use this single slice to calculate the medial side and 
angle of rotation for regional analysis.