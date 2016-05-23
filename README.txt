The code requires MATLAB version 2015a or later. Unzip all files in a folder and set the MATLAB current
directory to this folder or add it to MATLAB search path.

Instructions for use
--------------------
1. Type PhantomAnalysisGUI in MATLAB Command Window.
2. Press Load Data on the GUI and select a .tiff image sequence for analysis (e.g. from the validation
   dataset).
3. Fill the boxes in the Inputs panel if not filled automatically after loading the data. Images in the 
   validation dataset have voxel dimensions 0.075 x 0.075 x 1 mm. For this dataset, central slices 7 to 
   11 are often required for dimensions and volume estimation.
4. Press Measure dimensions in the Analysis panel to measure the phantom's dimensions, or Measure volume
   to measure the volume of its central cylindrical compartment. 

Notes:
- Canny edge detection was found to be the most accurate method for these analyses and is set as the 
  default, but k-means can be used instead. Choose the method you prefer in the Analysis panel. Change 
  automatically estimated Canny parameters or choose 2 or 3 clusters for k-means in the Inputs panel.
  Then, press Measure dimensions or volume accordingly.
- Increasing the scale factor for sub-pixel analysis (applied to dimensions only) can increase the error
  rate and analysis time significantly. A factor of 5 is the default, meaning that the precision of 
  measurements in any direction is 5 times higher than the voxel size in that direction.
- The Log window lists all the measurements performed thus far. More details of analysis are shown in 
  MATLAB Command Window.
- The analysis is visualized in the Output Data figure. Select Show Overlaid (edges) or Masks only. The 
  red lines represent the detected edges and the small green cross marks the positions where dimensional 
  measurements are performed; the distance between opposing cross marks is the measured dimension. 
- Check Save images? to save images showing the analysis process and final masks in a new folder in the 
  same location where the input image is. The analysis process is shown on the middle slice of those set
  in the Inputs panel.

Validation dataset
------------------
The validation dataset was created by randomly scaling and translating actual phantom scans and then 
randomly adding (Rician distribution) or removing noise (median filters). The dataset includes:
- scans acquired using a surface (01-30) and a volume coil (31-60)
- scans in the axial (01-05, 31-35), axial 90 (06-10, 36-40), coronal (11-15, 41-45), 
  coronal 90 (16-20, 46-50), sagittal (21-25, 51-55) and sagittal 90 plane (26-30, 56-60)

All scans were used for validating the tool's accuracy of dimension measurement, and axial scans alone
for validating the tool's accuracy of volume measurement. Although the tool is designed to be able to
measure the volume in scans in the other planes as well, it might have a higher chance of error. For 
both analyses, axial 90 and sagittal scans are rotated 90 degrees clockwise.



by Xenios Milidonis, last modified 23 May 2016
email: x.milidonis@sms.ed.ac.uk