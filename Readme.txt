This repository contains source code for "Building a user interface with 
MATLAB GUIDE for MRI data volumes in Imiomics",
Degree Project in Computer and Information Engineering,
Anna Larsson spring 2018.

To run the UI two groups of code are needed:
* The repository code, which is made by me with the exception of two scripts.
* Confidential files - image volumes and other medical data - which was only
stored locally on the desktop where I worked on this project.

I will publish the repository code to demonstrate my work, but it will be 
separated from all medical data (the "confidential files"). It should be 
possible to use the UI for a similar dataset, but depending on file naming 
and folder structure, some changes in the repository code may be needed.

If you have access to the same dataset (the confidential files) and wish to
run the UI, first check the path definitions in the beginning of the following 
code files and change them if needed:
* create_data.m
* del2_script.m  
Then follow the instructions under "Running the UI in the current setup".

If you wish to run the UI on a new (but similar) dataset, study the tables 
of the current data in the "Confidential data" section. Either rename your 
files or change the code (especially in create_data.m, del2.m and del2_script.m)
as needed. Make sure you have non-imaging data following the current format 
(tables need to have a certain number of columns, but the number of rows can vary 
depending on the size of the dataset).
 

==========
The code files in this repository:
==========

create_data.m       Run once for every new dataset, to create data structures.
del1_script.m       Is called by del1.m
del1.m              Is called by del1.fig (when running the UI part 1)
del1.fig            Opened by GUIDE
del2_script.m       Is called by del2.m
del2.m              Is called by del2.fig (when running the UI part 2)
del2.fig            Opened by GUIDE
imiFatVol.m         Storage class used in part 1
vtk_read_header.m   Volume loader script (not made by me)
vtk_read_volume.m   Volume loader script (not made by me)


==========
Running the UI in the current setup:
==========

1)

The folder containing the confidential files (large image volume files and 
non-imaging data) must be added to the MATLAB path. This will be referred
to as the "confidential data folder". In the current setup, this is E:\imiomics.

First:
In the Current Folder window in MATLAB, find the confidential data folder 
(currently E:\imiomics). Right-click on it and select 'Add to path' and 'Selected 
Folders and Subfolders'.
(Additionally, the "Set Path" option in MATLAB lets you add a folder to the path for 
both current and future sessions.) 

Then:
In the MATLAB Current Folder window, navigate to the code repository (the folder
where this Readme file is). This should be the current folder while using the UI.

2)

The script create_data.m generates a file named uidata.mat, which is placed in the 
confidential data folder. It can be described as a compressed version of the dataset.
If uidata.mat already exists for this dataset, skip to (3).

Otherwise, run the create_data.m script in MATLAB to generate the uidata.mat file.
In the current setup, this script takes about 10 minutes.

3) 

In MATLAB, type 'guide' to open the MATLAB GUIDE window. 
In GUIDE, open 'del1.fig' (for part 1) and 'del2.fig' (for part 2).
In GUIDE, press 'Run Figure' (green play button) to open the GUI. Do this for both
part 1 and part 2.
(Do not use the play button in the m-file editor in the regular MATLAB window.)
Part 1 should take about 2 minutes to open. Part 2 takes only a few seconds.


==========
Confidential data
==========

Here follows a listing of all the "confidential files" that are found in the 
"confidential data folder" (in the current setup this is E:\imiomics) and its 
subfolders. These are image volumes in the vtk-format along with non-imaging data 
stored in the mat-format.

These image volumes come from the "POEM study" and 500NNN represents a 
"POEM ID"-number, which identifies the subject.

E:\imiomics
-----------
(folder:) \correlation maps
    (8 correlation maps, named like the following:)
    correlation_analysis_female_pval_parametric_fat_dxa_n156.vtk
    ...

(folder:) \POEM data
    (300+ image volumes, named like the following:)
    fat_500XXX_500NNN.vtk
    ...
    JacDet_500XX_500NNN.vtk
    ...

(folder:) \POEM data undeformed
    (300+ image volumes, named like the following:)
    500NNN_fat_content.vtk
    ...

IDs.mat
    Table with 2 columns:
    Column 1: POEM ID (an integer in the form 500XXX)
    Column 2: 0 or 1. 1 if this POEM ID was included in the Imiomics analysis 
        (which created the correlation maps), 0 otherwise.

POEMdxa.mat
    Table with 5 columns:
    Column 1: POEM ID (an integer in the form 500XXX)
    Column 2: 1 or 0. 1 if subject is female, 0 if male.
    Column 3: Equal to column 1 if column 4 is a number, otherwise NaN.
    Column 4: DXA fat value, a float number (roughly 10^4). Is NaN if no data 
        exists for this subject.
    Column 5: 0 if column 4 is a number, otherwise a negative value.

refFemaleID.mat
    A POEM ID number converted to the string format, like '500NNN'.
    Identifies the female reference subject (used in part 2 of the UI).


refMaleID.mat
    A POEM ID number converted to the string format, like '500NNN'.
    Identifies the male reference subject (used in part 2 of the UI).

uidata.mat
    A file that is generated by the create_data.m script (based on all the 
    above files), and contains data needed for part 1 of the UI.


