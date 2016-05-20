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

#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkPlaneSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataSetMapper.h>
#include <vtkCenterOfMass.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

std::vector<std::vector<double> > ReadPointcloudFromPLYFile(std::string filename)
{
  std::vector< std::vector<double> > mylandmarks;
  std::ifstream ifs(filename.c_str(), std::ifstream::in);
  std::string line;
  if(ifs.is_open())
  {
      int lno = 0;
	  int vno;
	  int startlno = 1000000;
	  while(ifs.good())
	  {
		  lno++;
		  std::getline(ifs, line);
		  if(line.find("element vertex") != std::string::npos)
		  {
			  std::string vnostr = line.substr(15);
			  vno = stoi(vnostr);
			  //std::cout << "Vertex # = " << vno << std::endl;
		  }

		  if(line == "end_header")
		  {
			  startlno = lno+1;
		  }
		  
		  if(lno >= startlno && lno < vno+startlno)
		  {
              //std::cout << lno << ": " << line << std::endl;
			  std::stringstream ss(line);
			  std::string token;
			  std::vector<std::string> valstrs;
			  while(std::getline(ss, token, ' '))
			  {
				  valstrs.push_back(token);
				  //std::cout << token << std::endl;
			  }

			  std::vector<double> vertex;
			  for(int i=0; i<3; i++)
			  {
				  vertex.push_back(stod(valstrs[i]));
			  }

			  mylandmarks.push_back(vertex);
		  }
	  }
  }

  ifs.close();


  return mylandmarks;

}

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
  if(argc != 4)
    {
    std::cerr << "Usage: " << argv[0]
              << " Filename(.ply) Fitting_plane_filename(.txt) method" << std::endl;
    return EXIT_FAILURE;
    }

  std::string filename = argv[1];
  std::string plane_filename = argv[2];
  int method = atoi(argv[3]);

  //----------------------------------------------------------------------//
  // read the *.ply data 
  // ---------------------------------------------------------------------//
  std::vector< std::vector<double> > mylandmarks = ReadPointcloudFromPLYFile(filename);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();

  for(int i=0; i<mylandmarks.size(); i++)
  {
	  double p[3];
      std::cout << "i = " << i << ": ";
	  
	  for(int j=0; j<3; j++)
	  {
		  p[j] = mylandmarks[i][j];
	      std::cout << " " << mylandmarks[i][j];
	  }
	  std::cout <<  std::endl;

	  vtkIdType pid[1];
	  pid[0] = points->InsertNextPoint(p);
	  vertices->InsertNextCell(1,pid);
  }

  std::cout<< " check here 2. " << std::endl;
  vtkSmartPointer<vtkPolyData> point = vtkSmartPointer<vtkPolyData>::New();
  point->SetPoints(points);
  point->SetVerts(vertices);



  //-------------------------------------------------------------//
  // Read the parameter for the fitting plane.
  //-------------------------------------------------------------//
  std::vector< std::vector<double> > myplanes = ReadFitParamFile(plane_filename);
  double planeParams[4];
  for(int i=0; i<myplanes.size(); i++)
  {
      std::cout << "Fitting plane normal " << i << ": ";
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



  //-------------------------------------------------------------//
  // Create a plane for visualization.
  // ------------------------------------------------------------//

  // Compute the center of mass
  vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter = vtkSmartPointer<vtkCenterOfMass>::New();
  centerOfMassFilter->SetInputData(point);
  centerOfMassFilter->SetUseScalarsAsWeights(false);
  centerOfMassFilter->Update();

  double center[3];
  centerOfMassFilter->GetCenter(center);
	  
  std::cout << "Center of mass is " << center[0] << " " << center[1] << " " << center[2] << std::endl;

  // Create a plane
  double dist = -planeParams[3];
  vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetCenter(center[0] + planeParams[0]*dist, center[1] + planeParams[1]*dist, center[2] + planeParams[2]*dist);
  planeSource->SetNormal(planeParams[0], planeParams[1],  planeParams[2]);
  
  
  planeSource->SetXResolution(20);
  planeSource->SetYResolution(20);
  planeSource->Update();

  vtkPolyData* planedata = planeSource->GetOutput();


  
  //-----------------------------------------------------------//
  //  Visualization.
  //-----------------------------------------------------------//

  // Create a mapper and actor for landmarks
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(point);


  vtkSmartPointer<vtkActor> actor =  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetPointSize(3);
  actor->GetProperty()->SetColor(0.0, 0.0, 1.0);


  // Create a mapper and actor for the plane
  vtkSmartPointer<vtkPolyDataMapper> planeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  planeMapper->SetInputData(planedata);

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
  renderer->SetBackground(0.3, 0.6, 0.3); // Background color white
 
  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
  
  return EXIT_SUCCESS;
}
