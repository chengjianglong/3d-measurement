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

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkMath.h>

#include <iostream>
#include <sstream>
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

int main(int argc, char *argv[])
{ 
  // Parse command line arguments
  if(argc != 4)
    {
    std::cerr << "Usage: " << argv[0]
              << " Filename(.vtp) Fitting_plane_filename(.txt) method" << std::endl;
    return EXIT_FAILURE;
    }

  std::string filename = argv[1];
  std::string plane_filename = argv[2];
  int method = atoi(argv[3]);

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


  //----------------------------------------------------//
  // Read the data form the input file.
  //----------------------------------------------------//

  vtkSmartPointer<vtkPolyData> inputPolyData;
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  inputPolyData = reader->GetOutput();
  

  vtkSmartPointer<vtkPolyDataMapper> inputMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  inputMapper->SetInput(inputPolyData);
#else
  inputMapper->SetInputData(inputPolyData);
#endif
  

  //---------------------------------------------------//
  // Create a plane to cut
  //---------------------------------------------------//
  //
  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
  plane->SetOrigin(inputPolyData->GetCenter());
  plane->SetNormal(planeNormal[0], planeNormal[1], planeNormal[2]);
  
  double minBound[3];
  minBound[0] = inputPolyData->GetBounds()[0];
  minBound[1] = inputPolyData->GetBounds()[2];
  minBound[2] = inputPolyData->GetBounds()[4];

  double maxBound[3];
  maxBound[0] = inputPolyData->GetBounds()[1];
  maxBound[1] = inputPolyData->GetBounds()[3];
  maxBound[2] = inputPolyData->GetBounds()[5];

  double center[3];
  center[0] = inputPolyData->GetCenter()[0];
  center[1] = inputPolyData->GetCenter()[1];
  center[2] = inputPolyData->GetCenter()[2];

  double distanceMin = sqrt(vtkMath::Distance2BetweenPoints(minBound, center));
  double distanceMax = sqrt(vtkMath::Distance2BetweenPoints(maxBound, center));

  std::cout << "distanceMin : " << distanceMin << std::endl;
  std::cout << "distanceMax : " << distanceMax << std::endl;

  // Create cutter
  vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
  cutter->SetCutFunction(plane);
#if VTK_MAJOR_VERSION <= 5
  cutter->SetInput(inputPolyData);
#else
  cutter->SetInputData(inputPolyData);
#endif
  cutter->GenerateValues(20, -distanceMin, distanceMax);
  vtkSmartPointer<vtkPolyDataMapper> cutterMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  cutterMapper->SetInputConnection( cutter->GetOutputPort());
  cutterMapper->ScalarVisibilityOff();



  //-------------------------------------------------------------------//
  // Visualiztion.
  //-------------------------------------------------------------------//
  
  // Create plane actor
  vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
  planeActor->GetProperty()->SetColor(1.0,1,0);
  planeActor->GetProperty()->SetLineWidth(3);
  planeActor->SetMapper(cutterMapper);
  
  // Create input actor
  vtkSmartPointer<vtkActor> inputActor = vtkSmartPointer<vtkActor>::New();
  inputActor->GetProperty()->SetColor(1.0, 0.8941, 0.7686); // bisque
  inputActor->SetMapper(inputMapper);
  
  // Create renderers and add actors of plane and cube
  vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(planeActor); //display the rectangle resulting from the cut
  renderer->AddActor(inputActor); //display the cube
  
  //Add renderer to renderwindow and render
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(600, 600);
  
  vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);
  renderer->SetBackground(.3, .6, .3);
  renderWindow->Render();

  interactor->Start();
  
  return EXIT_SUCCESS;
}
