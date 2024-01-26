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
	
The scan 1, 2, and difference figures can be navigated by slice #

The BV, BMC, BMD variables are tables that report TOTAL bone volume,
	bone mineral content, and bone mineral density in the anterior
	medial, anterior lateral, posterior medial, and posterior lateral 
	quadrants, respectively.
	
Voxel resolution should be changed in the code if necessary, 
and the calibration slope and intercept should be machine-specific.

To setup a batch run, change the subjects initialization line to 
include all of the subjects you want to process.
E.g. subjects = ["MARA14-","MARA15-","MARA16-"];

PLEASE NOTE: Each time you run this code a new output.xlsx will be
written to the code's home directory. This WILL be overwritten if 
you run the code again.

Also, the excel document does not save the numbers as numbers, but 
instead as text. You'll need to highlight the numbers and click the 
warning symbol that pops up to choose "convert to number"
