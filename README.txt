/*---------------------------------------------------------------------------
 *
 * Copyright 2016 by Kitware, Inc. All Rights Reserved. Please refer to
 *
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 *
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 *
 * Author : Chengjiang Long <chengjiang.long@kitware.com>
 *
 * Description :
 *
 *
 * Created : <2016-05-10>
 *
 *-------------------------------------------------------------------------*/

-----------------------
How to build?
-----------------------
Once the "measurement" directory  has been "git clone",
	 cd measurement/build 
	 If VTK is installed:
	 cmake ..
	 If VTK is not installed but compiled on your system, you will need to specify the path to your
	 VTK build:
	 cmake -DVTK_DIR:PATH=/home/me/vtk_build ..
	 Build the project:
	 make

Attention: I use VTK 7.0 to develop. So if your version of VTK is lower than 7.0, some functions may
not be well supported.

-----------------------
How to use?
-----------------------


1. visualize the *.vtp polydata.
./VisualizePolyData *.vtp

The usage of ./ReadPolyData is the same as ./VisualizePolyData. The difference is that
./VisualizePolyData is able to display the color information.


2. visualize the fitting plane
./VisualizeFitPlane fitting_plane.txt method
For method:
0: estimated by least square regression.
1: estimated by least median of squares.
2: estimated by RANSAC.

3. smooth the the *.vtp polydata
./SmoothPolyData input.vtp output.vtp
e.g.: ./SmoothPolyData output_contour_400_coloration.vtp smooth.vtp

4. edit the *.vtp polydata
./EditPolyData input.vtp fitting_plane.txt output.vtp method colorThreshold normalThreshold
e.g.: ./EditPolyData smooth.vtp fitting_plane.txt edit.vtp 2 100 0.35

5. extract the object with the fitting_plane
./ExtractObject input.vtp fitting_plane.txt method distThreshold
e.g.: ./ExtractObject edit.vtp fitting_plane 2 0.22
Note that footseg.vtp, footcut.vtp and footcut.ply come out after taking this step. footseg.vtp is
the extracted foot, while footcut.vtp and footcut.ply are the point clouds obtained by the cutter 
plane.

The usage of ./CutPolyData is the same to ./ExtractObject.

To better under stand the cutter coutour, I would like to suggest you to run ./ContoursFromPolyData.
It usage is as follows:
./ContoursFromPolyData input.vtp fitting_plane.txt method
e.g.: ./ContoursFromPolyData footseg.vtp fitting_plane.txt 2

6. measurement with the fitting line by MAPTK tool (maptk_fit_line -c maptk_fit_line.conf)
./Measure countour.vtp fitting_plane.txt  fitting_line.txt method
e.g.: ./Measure footcut.vtp fitting_plane.txt  fitting_line.txt 2



---------------------
Remarks
---------------------
(a) fitting_plane.txt is obtained by running MAPTK tool
maptk_fit_plane -c maptk_fit_plane.conf 

*Note that configure file mapt_fit_plane.conf must mention that landmarks.ply is the input, and fitting_plane.txt is the output.

(b) fitting_line.txt is also obtained by running MAPTK tool 
maptk_fit_line -c maptk_fit_line.conf

*Note that configure file mapt_fit_line.conf must mention that contour point cloud (e.g.: footcut.ply) is the input, and fitting_line.txt is the output.

