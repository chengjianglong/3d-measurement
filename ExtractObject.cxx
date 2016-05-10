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

#include <vtkPLYReader.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>

#include "vtkUnstructuredGrid.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include <vtkGeometryFilter.h>
#include <vtkConnectivityFilter.h>
#include <vtkCutter.h>
#include <vtkStripper.h>
#include <vtkSpline.h>
#include <vtkKochanekSpline.h>
#include <vtkSplineFilter.h>
#include <vtkPLYWriter.h>

#include <vtkPlaneSource.h>
#include <vtkDataSetWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkCenterOfMass.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkImplicitDataSet.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkPlane.h>
#include <vtkUnstructuredGrid.h>
#include <vtkBridgeDataSet.h>
#include <vtkCellData.h>
#include <vtkGenericClip.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

void WritePolyData(vtkPolyData* const polyData, const std::string& filename)
{
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polyData);
#else
    writer->SetInputData(polyData);
#endif
    writer->SetFileName(filename.c_str());
    writer->Write();
}
 
void WriteDataSet(vtkDataSet* const dataSet, const std::string& filename)
{
    vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(dataSet);
#else
    writer->SetInputData(dataSet);
#endif
    writer->SetFileName(filename.c_str());
    writer->Write();
}

int main ( int argc, char *argv[] )
{
  // Parse command line arguments
  if(argc != 5)
    {
    std::cerr << "Usage: " << argv[0]
              << " Filename(.vtp) Fitting_plane_filename(.txt) method threshold" << std::endl;
    return EXIT_FAILURE;
    }

  std::string filename = argv[1];
  std::string plane_filename = argv[2];
  int method = atoi(argv[3]);
  std::string strdist = argv[4];

  //-------------------------------------------------------------//
  // Read the parameter for the fitting plane.
  //-------------------------------------------------------------//
  std::vector< std::vector<double> > myplanes;
  std::ifstream plane_ifs(plane_filename.c_str(), std::ifstream::in);
  std::string strline;
  if(plane_ifs.is_open())
  {
      int lno = 0;
	  int vno;
	  while(plane_ifs.good())
	  {
		  lno++;
		  std::getline(plane_ifs, strline);
		  
		  if(lno == 2 || lno == 5 || lno == 8)
		  {
              //std::cout << lno << ": " << line << std::endl;
			  std::stringstream ss(strline);
			  std::string token;
			  std::vector<std::string> valstrs;
			  while(std::getline(ss, token, ' '))
			  {
				  valstrs.push_back(token);
				  //std::cout << token << std::endl;
			  }

			  std::vector<double> params;
			  for(int i=0; i<4; i++)
			  {
				  params.push_back(stod(valstrs[i]));
			  }

			  myplanes.push_back(params);
		  }
	  }
  }

  plane_ifs.close();


  double planeParams[4];
  for(int i=0; i<myplanes.size(); i++)
  {
      std::cout << "Plane normal = " << i << ": ";
	  for(int j=0; j<myplanes[i].size(); j++)
	  {
	      std::cout << " " << myplanes[i][j];
		  
		  if(i==method && j<4)
		  {
			  planeParams[j] = myplanes[i][j];
		  }
	  }
	  std::cout <<  std::endl;
  }

  //----------------------------------------------------//
  // adjust the normal direction.
  // Based on the observation, when z-value in planeParams
  // is positive, then it is the upper plane.
  //----------------------------------------------------//
  if(planeParams[2] < 0)
  {
	  planeParams[0] = -planeParams[0];
	  planeParams[1] = -planeParams[1];
	  planeParams[2] = -planeParams[2];
	  planeParams[3] = -planeParams[3];
  }


  //----------------------------------------------------------//
  // Read the input polydata.
  //----------------------------------------------------------//
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  // Compute the center of mass
  vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter = vtkSmartPointer<vtkCenterOfMass>::New();
  centerOfMassFilter->SetInputConnection(reader->GetOutputPort());
  centerOfMassFilter->SetUseScalarsAsWeights(false);
  centerOfMassFilter->Update();

  double center[3];
  centerOfMassFilter->GetCenter(center);
	  
  std::cout << "Center of mass is " << center[0] << " " << center[1] << " " << center[2] << std::endl;

  // Add ids to the points and cells
  vtkSmartPointer<vtkIdFilter> idFilter = vtkSmartPointer<vtkIdFilter>::New();
  idFilter->SetInputConnection(reader->GetOutputPort());
  idFilter->Update();
  
  // Create a plane to clip with
  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

  double dist = stod(strdist);
  plane->SetOrigin(center[0] + planeParams[0]*dist, center[1] + planeParams[1]*dist, center[2] + planeParams[2]*dist);
  plane->SetNormal(planeParams[0], planeParams[1], planeParams[2]);

  // Create a plane
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetCenter(center[0], center[1], center[2]);
  // RANSIC
  planeSource->SetNormal(planeParams[0], planeParams[1], planeParams[2]);
  
  // LMS
  //planeSource->SetNormal(0.183286, 0.0549184, -0.981524);
  
  // Regression
  //planeSource->SetNormal(-0.688477, 0.318007, 0.651822);
  
  
  //planeSource->SetOutputPointsPrecision(6);
  planeSource->Update();

  vtkPolyData* planedata = planeSource->GetOutput();

  // Covert the DataSet to a GenericDataset
  vtkSmartPointer<vtkBridgeDataSet> bridgeDataSet = vtkSmartPointer<vtkBridgeDataSet>::New();
  bridgeDataSet->SetDataSet(idFilter->GetOutput());

  vtkSmartPointer<vtkGenericClip> clipper = vtkSmartPointer<vtkGenericClip>::New();
  clipper->SetClipFunction(plane);
  clipper->SetInputData(bridgeDataSet);
  clipper->Update();

// Get the clipped cell ids
//  vtkUnstructuredGrid* clipped = clipper->GetOutput();
//  vtkIdTypeArray* originalIds = vtkIdTypeArray::SafeDownCast(clipped->GetCellData()->GetArray("vtkIdFilter_Ids"));
   
//    {
//    std::cout << "new id " << i << ", original id " << originalIds->GetValue(i) << std::endl;
//    }
  
  
  
  //-----------------------------------------------------------------//
  //  Save the results
  //-----------------------------------------------------------------//
  //
  // convert the unstructured grid into polydata
  vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  geometryFilter->SetInputData(clipper->GetOutput());
  geometryFilter->Update(); 

  /*
  vtkXMLPolyDataWriter *writer=vtkXMLPolyDataWriter::New();
  writer->SetInputConnection(geometryFilter->GetOutputPort());
  //writer->SetFileName("clipped.vtu");
  writer->SetFileName("clipper.vtp");
  writer->SetDataModeToAscii();
  writer->Write();
  writer->Delete();
  */

  // Write the file
  WritePolyData(planedata, "planedata.vtp");
  //WritePolyData(geometryFilter->GetOutput(), "clipper.vtp"); 

  //------------------------------------------------------------------//
  // Extract the color information.
  // -----------------------------------------------------------------//
  //

  // Extract the polydata
  vtkSmartPointer<vtkPolyData> output = geometryFilter->GetOutput();
  
  // Get the number of points in the polydata
  vtkIdType idNumPointsInFile = output->GetNumberOfPoints();

  std::cout << "NumPoint: " << idNumPointsInFile << std::endl;

  // access the color information
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  // Now check for point data
  vtkPointData *pd = output->GetPointData();
  if (pd)
  {
      std::cout << " contains point data with "
                << pd->GetNumberOfArrays()
                << " arrays." << std::endl;

      for (int i = 0; i < pd->GetNumberOfArrays(); i++)
      {
           std::cout << "\tArray " << i
                << " is named "
                << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL")
                << std::endl;

		   //vtkSmartPointer<vtkDataArray> arrayData =  vtkSmartPointer<vtkDataArray>::New();
		   //arrayData = pd->GetArray(pd->GetArrayName(i));
		   vtkSmartPointer<vtkDataArray> arrayData = pd->GetArray(pd->GetArrayName(i));

		   if(i==1)
		   {
			   //colors->DeepCopy(arrayData);
			   
			   
			   for(int j=0; j<arrayData->GetNumberOfTuples(); j++)
			   {
				   unsigned char color[3];
				   //for(int k=0; k<arrayData->GetElementComponentSize(); k++)
				   for(int k=0; k<3; k++)
				   {
					   color[k] = arrayData->GetComponent(j,k);
					   //std::cout <<"("<< i << " , " << j << " , " << k << "): " << arrayData->GetComponent(j,k) << std::endl;
				   }

				   //std::cout<<"color: " << (int)color[0] << " " << (int)color[1] << " " << (int)color[2] << std::endl;

				   //colors->InsertNextTypedTuple(color);
				   colors->InsertNextTuple3(color[0], color[1], color[2]);
			   }
			

		   }


       }
  }

  std::cout << "color # : " << colors->GetNumberOfTuples() << std::endl;

  output->GetPointData()->SetScalars(colors);


  // Apply vtkConnectivityFilter
  vtkSmartPointer<vtkConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
  connectivityFilter->SetInputData(output);
  connectivityFilter->SetExtractionModeToLargestRegion();
  //connectivityFilter->SetExtractionModeToAllRegions();
  //connectivityFilter->ColorRegionsOn();
  connectivityFilter->Update();

  // Apply vtkCutter
  vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
  cutter->SetCutFunction(plane);
  cutter->SetInputConnection(connectivityFilter->GetOutputPort());
  //cutter->GenerateValues(20, -0.816385, 0.816385);
  cutter->GenerateValues(1, 0, 0);
  
  // Create a mapper and actor for clipped points
  vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  //mapper->SetInputConnection(clipper->GetOutputPort());
  // mapper->SetInputData(output);
  mapper->SetInputConnection(connectivityFilter->GetOutputPort());

  vtkSmartPointer<vtkGeometryFilter> clippergeometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  clippergeometryFilter->SetInputData(connectivityFilter->GetOutput());
  clippergeometryFilter->Update(); 

  vtkXMLPolyDataWriter *writer=vtkXMLPolyDataWriter::New();
  writer->SetInputConnection(clippergeometryFilter->GetOutputPort());
  //writer->SetFileName("clipped.vtu");
  writer->SetFileName("clipper.vtp");
  writer->SetDataModeToAscii();
  writer->Write();
  writer->Delete();
  
  //-----------------------------------------------------------------//
  //  Save the largest commponent
  //-----------------------------------------------------------------//
  //
  
  vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
  stripper->SetInputConnection(cutter->GetOutputPort());

  vtkSmartPointer<vtkKochanekSpline> spline = vtkSmartPointer<vtkKochanekSpline>::New();
  spline->SetDefaultTension(0.5);

  vtkSmartPointer<vtkSplineFilter> sf = vtkSmartPointer<vtkSplineFilter>::New();
  sf->SetInputConnection(stripper->GetOutputPort());
  sf->SetSubdivideToSpecified();
  sf->SetNumberOfSubdivisions(50);
  sf->SetSpline(spline);
  sf->GetSpline()->ClosedOn();

  vtkSmartPointer<vtkGeometryFilter> cutgeometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  cutgeometryFilter->SetInputConnection(cutter->GetOutputPort());
  //cutgeometryFilter->SetInputConnection(sf->GetOutputPort());
  cutgeometryFilter->Update(); 

  vtkXMLPolyDataWriter *cutwriter=vtkXMLPolyDataWriter::New();
  cutwriter->SetInputConnection(cutgeometryFilter->GetOutputPort());
  //writer->SetFileName("clipped.vtu");
  cutwriter->SetFileName("footcut.vtp");
  cutwriter->SetDataModeToAscii();
  cutwriter->Write();
  cutwriter->Delete();
 
  vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
  plyWriter->SetFileName("footcut.ply");
  plyWriter->SetFileTypeToASCII();
  plyWriter->SetInputConnection(cutgeometryFilter->GetOutputPort());
  plyWriter->Write();

  // Create a mapper and actor for cutter
  vtkSmartPointer<vtkPolyDataMapper> cutterMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  //cutterMapper->SetInputConnection(cutter->GetOutputPort());
  cutterMapper->SetInputConnection(sf->GetOutputPort());
  cutterMapper->ScalarVisibilityOff();
  
 
  // convert the unstructured grid into polydata
  //vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  geometryFilter->SetInputConnection(connectivityFilter->GetOutputPort());
  geometryFilter->Update(); 

  vtkXMLPolyDataWriter *segwriter=vtkXMLPolyDataWriter::New();
  segwriter->SetInputConnection(geometryFilter->GetOutputPort());
  //writer->SetFileName("clipped.vtu");
  segwriter->SetFileName("footseg.vtp");
  segwriter->SetDataModeToAscii();
  segwriter->Write();
  segwriter->Delete();
 
  //-----------------------------------------------------------------//
  //  Visualize
  //-----------------------------------------------------------------//
  //
  vtkSmartPointer<vtkActor> actor =  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  //actor->GetProperty()->SetPointSize(2);

  // Create a mapper
  vtkSmartPointer<vtkPolyDataMapper> planeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  planeMapper->SetInputData(planedata);

  // Create an actor for the plane
  vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
  planeActor->SetMapper(planeMapper);

  // Create a actor for the cutter
  vtkSmartPointer<vtkActor> cutterActor = vtkSmartPointer<vtkActor>::New();
  cutterActor->SetMapper(cutterMapper);
  cutterActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
  cutterActor->GetProperty()->SetLineWidth(6);

  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  // Add the actor to the scene
  renderer->AddActor(actor);
  renderer->AddActor(planeActor);
  renderer->AddActor(cutterActor);
  renderer->SetBackground(0.3,0.6,0.3); // Background color white
 
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
  
  return EXIT_SUCCESS;
}
