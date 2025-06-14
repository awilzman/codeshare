$!      
$!         _/_/_/  _/_/_/    _/        
$!          _/    _/    _/  _/           Image Processing Language  
$!         _/    _/_/_/    _/ 
$!        _/    _/        _/             (c) Andres Laib, IBT/ETH Z�rich
$!     _/_/_/  _/        _/_/_/_/                         Scanco Medical
$!        
$!
$!
$!  IPL Batch Scanco
$!  This script was written by Karen Troy June 2025
$
$!  To run, type @UE:IPL_TRANSFORM_AIM.COM (or whatever the current script name is,
$!  into the command window
$
$!  Purpose is to open up an AIM file, transform it using an existing TMAT.DAT file,
$!  and then write the transformed AIM file and also export it to DICOM.
$!  This program also writes the original AIM (not-transformed) into DICOMs

$
$!Modify these each time you run the code
$   orig_dir="DISK2:[MICROCT.IMPORT]"
$   orig_name="C0001992_4_1"  
$   follow_dir="DISK2:[MICROCT.IMPORT]"
$   follow_name="C0002216_4_2"
$   tmat_name="_30_TMAT"
$   trans_dcm_name="C0002216_4_2_trans"
$
$!Files you are reading
$   raw_orig_file:='orig_dir''orig_name'.AIM
$   orig_dcm:='orig_dir''orig_name'
$   raw_follow_file:='follow_dir''follow_name'.AIM
$   tmat_file:='follow_dir''tmat_name'.DAT
$   follow_trans_dcm:='follow_dir''trans_dcm_name'
$
$!Files you are writing
$   transformed_file:='follow_dir''trans_dcm_name'.AIM
$
$ipl_scanco_prog := $um:ipl_scanco_fe.exe
$
$ON ERROR THEN EXIT
$ipl_scanco_prog 

/todicom_from_aim
  -aim_name			"raw_orig_file
  -dicom_filename 		"orig_dcm
  -multifiles			true
  -global_slicename		true

/read
  -name             	     raw_follow
  -filename  		     "raw_follow_file

/transform
  -in				raw_follow
  -out				raw_follow_trans
  -Tmat_file_name		"tmat_file
  -img_interpol_option		1
  -el_size_mm_out		-1.000 -1.000 -1.000
  -Tmat_invert			false

/write_v020
  -name				raw_follow_trans
  -filename			"transformed_file
  -compress_type		bin
  -version_020			true

/todicom_from_aim
  -aim_name			"transformed_file
  -dicom_filename 		"follow_trans_dcm
  -multifiles			true
  -global_slicename		true


..
$ exit                                      