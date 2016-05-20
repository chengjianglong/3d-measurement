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
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>

#include <vtkPlaneSource.h>
#include <vtkPointData.h>
#include <vtkGeometryFilter.h>

#include <vtkXMLPolyDataReader.h>
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

std::vector<std::vector<double> > ReadFitParamFile(std::string filename)
{
  std::vector< std::vector<double> > myparams;
  std::ifstream ifs(filename.c_str(), std::ifstream::in);
  std::string strline;
  if(ifs.is_open())
  {
      int lno = 0;
	  int vno;
	  while(ifs.good())
	  {
		  lno++;
		  std::getline(ifs, strline);
		  
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
			  for(int i=0; i<valstrs.size(); i++)
			  {
				  params.push_back(stod(valstrs[i]));
			  }

			  myparams.push_back(params);
		  }
	  }
  }

  ifs.close();

  return myparams;

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
  std::vector< std::vector<double> > myplanes = ReadFitParamFile(plane_filename);
  double planeNormal[3];
  for(int i=0; i<myplanes.size(); i++)
  {
      std::cout << "Plane normal = " << i << ": ";
	  for(int j=0; j<myplanes[i].size(); j++)
	  {
	      std::cout << " " << myplanes[i][j];
		  
		  if(i==method && j<3)
		  {
			  planeNormal[j] = myplanes[i][j];
		  }
	  }
	  std::cout <<  std::endl;
  }

  //----------------------------------------------------//
  // adjust the normal direction.
  // Based on the observation, when z-value in planeNormal
  // is positive, then it is the upper plane.
  //----------------------------------------------------//
  if(planeNormal[2] < 0)
  {
	  planeNormal[0] = -planeNormal[0];
	  planeNormal[1] = -planeNormal[1];
	  planeNormal[2] = -planeNormal[2];
  }

  //-----------------------------------------------------//
  // Read the polydata.
  //-----------------------------------------------------//
  // Read all the data from the file
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  //---------------------------------------------------//
  // Extract the center.
  //---------------------------------------------------//

  // Compute the center of mass
  vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter = vtkSmartPointer<vtkCenterOfMass>::New();
  centerOfMassFilter->SetInputConnection(reader->GetOutputPort());
  centerOfMassFilter->SetUseScalarsAsWeights(false);
  centerOfMassFilter->Update();

  double center[3];
  centerOfMassFilter->GetCenter(center);
	  
  std::cout << "Center of mass is " << center[0] << " " << center[1] << " " << center[2] << std::endl;


  //---------------------------------------------------//
  // Clip the polydata.
  //---------------------------------------------------//
  
  // Add ids to the points and cells
  vtkSmartPointer<vtkIdFilter> idFilter = vtkSmartPointer<vtkIdFilter>::New();
  idFilter->SetInputConnection(reader->GetOutputPort());
  idFilter->Update();
  
  double dist = stod(strdist);
  center[0] += dist*planeNormal[0];
  center[1] += dist*planeNormal[1];
  center[2] += dist*planeNormal[2];
  
  //---------------------------------------------------//
  // Create a plane source.
  //---------------------------------------------------//
  
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetCenter(center[0], center[1], center[2]);
  planeSource->SetNormal(planeNormal[0], planeNormal[1], planeNormal[2]);
  planeSource->SetXResolution(20);
  planeSource->SetYResolution(20);
  planeSource->Update();
 
  //---------------------------------------------------//
  // Create a plane to clip with
  //---------------------------------------------------//
  //
  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
  plane->SetOrigin(center[0], center[1], center[2]);
  plane->SetNormal(planeNormal[0], planeNormal[1], planeNormal[2]);


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
 

  //------------------------------------------------------------------//
  // Extract the color information.
  // -----------------------------------------------------------------//
  // Convert the unstructured grid to polydata
  vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  geometryFilter->SetInputData(clipper->GetOutput());
  geometryFilter->Update();

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

		   vtkSmartPointer<vtkDataArray> arrayData = pd->GetArray(pd->GetArrayName(i));

		   if(strcmp(pd->GetArrayName(i), "medianColoration")==0 || strcmp(pd->GetArrayName(i),"MedianColoration")==0)
		   {
			   colors->DeepCopy(arrayData);
			   
			   
			   for(int j=0; j<arrayData->GetNumberOfTuples(); j++)
			   {
				   unsigned char color[3];
			
				   //std::cout << "\t ----- size: " << arrayData->GetElementComponentSize() << std::endl;

				   //for(int k=0; k<arrayData->GetElementComponentSize(); k++)
				   for(int k=0; k<3; k++)
				   {
					   color[k] = arrayData->GetComponent(j,k);
					   //std::cout <<"("<< i << " , " << j << " , " << k << "): " << arrayData->GetComponent(j,k) << std::endl;
				   }

				   //std::cout<<"color: " << (int)color[0] << " " << (int)color[1] << " " << (int)color[2] << std::endl;

				   //colors->InsertNextTuple3(color[0], color[1], color[2]);
			   }

		   }


       }
  }

  std::cout << "color # : " << colors->GetNumberOfTuples() << std::endl;

  output->GetPointData()->SetScalars(colors);
 

  //----------------------------------------------------------------//
  //  Visualization.
  //----------------------------------------------------------------//

  // Create a mapper and actor for clipped points
  vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  //mapper->SetInputConnection(clipper->GetOutputPort());
  mapper->SetInputData(output);


  vtkSmartPointer<vtkActor> actor =  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  //actor->GetProperty()->SetPointSize(2);

    
  // Create a mapper
  vtkSmartPointer<vtkPolyDataMapper> planeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  planeMapper->SetInputData(planedata);

  // Create an actor for the sphere
  vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
  planeActor->SetMapper(planeMapper);
  
  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  // Add the actor to the scene
  renderer->AddActor(actor);
  renderer->AddActor(planeActor);
  renderer->SetBackground(.3,.6,.3); // Background color white
 
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
  
  return EXIT_SUCCESS;
}
