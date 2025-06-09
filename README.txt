This version of the code reads in the full grayscale AIM files to correctly calculate vBMD.
Updated and tested by Karen Troy 6/6/2025

To run this code, you first need to export the full gray AIM files from the scanner.
This is accomplished by running the IPL_TRANSFORM_AIM.COM file.

In the COM file, set the following variables before you run it:
orig_dir  	the directory where your time 1 AIM file is located
orig_name  	the name of your original AIM file (without the AIM extension)
follow_dir	the directory where your time 2 AIM file is located
follow_name	the name of your time 2 AIM file (without the AIM extension)
tmat_name	the name of the transformation matrix that was generated 
		using image 3D registration (assumed to be located in the follow_dir 
		location
trans_dcm_name	the name you want your transformed time 2 AIM to have (I recommend 
		using the same as follow_name with "_trans" appended)


This file will generate 3 sets of files:
1. Dicom files for the time 1 AIM
2. A transformed time2 AIM file
3. Dicom files for the time 2 AIM

You then need to use filezilla to transfer the dicome files into folders for analysis
Name the folders '[scan name]_aim'.  For example, for scan 524_4_scan1, your folder
should be named 524_4_scan1_aim.  And the scan2 would be 524_4_scan2_aim

Each folder should contain all of the dicom files from the respective aim.  Only scan2
will be transformed.