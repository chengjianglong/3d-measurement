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

#include <sstream>
#include <iostream>

#include <vtkVersion.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkSphereSource.h>
#include <vtkExtractEdges.h>
#include <vtkTriangleFilter.h>

#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

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
  if(argc != 7)
    {
    std::cerr << "Usage: " << argv[0]
              << " filename(.vtp fitting_plane_filename(.txt) edit_filename(.vtp) method grayThreshold normalThreshold" << std::endl;
    return EXIT_FAILURE;
    }

  std::string filename = argv[1];
  std::string plane_filename = argv[2];
  std::string edit_filename = argv[3];
  int method = atoi(argv[4]);
  int grayThreshold = atoi(argv[5]);
  double normalThreshold = stod(std::string(argv[6]));


  //-------------------------------------------------------------//
  // Read the parameter for the fitting plane.
  //-------------------------------------------------------------//
  std::vector< std::vector<double> > myplanes = ReadFitParamFile(plane_filename);
  double normal[3];
  for(int i=0; i<myplanes.size(); i++)
  {
      std::cout << "Fitting plane  = " << i << ": ";
	  for(int j=0; j<myplanes[i].size(); j++)
	  {
	      std::cout << " " << myplanes[i][j];
		  
		  if(i==method && j<3)
		  {
			  normal[j] = myplanes[i][j];
		  }
	  }
	  std::cout <<  std::endl;
  }


  //----------------------------------------------------------//
  // Read the polydata (*.vtp)
  //----------------------------------------------------------//

  // Read all the data from the file
  vtkSmartPointer<vtkXMLPolyDataReader> reader =  vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();


  //------------------------------------------------------------------//
  // Extract the color information.
  // -----------------------------------------------------------------//
  // Extract the polydata
  vtkSmartPointer<vtkPolyData> output = reader->GetOutput();
  
  // Get the number of points in the polydata
  vtkIdType idNumPointsInFile = output->GetNumberOfPoints();

  std::cout << "NumPoint: " << idNumPointsInFile << std::endl;
  std::vector<bool> bDeleteList(idNumPointsInFile, false);

  // access the color information
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName(filename.c_str());

  vtkSmartPointer<vtkIdList> myIdList = vtkSmartPointer<vtkIdList>::New(); 
  
  // Now check for point data
  vtkPointData *pd = output->GetPointData();
  if (pd)
  {
      std::cout << "contains point data with "
                << pd->GetNumberOfArrays()
                << " arrays." << std::endl;


      for (int i = 0; i < pd->GetNumberOfArrays(); i++)
      {
           std::cout << "\tArray " << i
                << " is named "
                << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL")
                << std::endl;

		   vtkSmartPointer<vtkDataArray> arrayData = pd->GetArray(pd->GetArrayName(i));

		   //std::cout << "\t----- array_length : " << arrayData->GetNumberOfTuples() << std::endl;
		   
		   if(strcmp(pd->GetArrayName(i), "medianColoration")==0 || strcmp(pd->GetArrayName(i),"MedianColoration")==0)
		   {
			   //colors->DeepCopy(arrayData);
			   
			   
			   for(int j=0; j<arrayData->GetNumberOfTuples(); j++)
			   {
				   unsigned char color[3];
				   int totalcolor = 0;
				   //for(int k=0; k<arrayData->GetElementComponentSize(); k++)
				   for(int k=0; k<3; k++)
				   {
					   color[k] = arrayData->GetComponent(j,k);
					   totalcolor += (int)color[k];
					   //std::cout <<"("<< i << " , " << j << " , " << k << "): " << arrayData->GetComponent(j,k) << std::endl;
				   }

				   // Select the vertexes by the gray color information. 
				   if(totalcolor < grayThreshold)
				   {
					   myIdList->InsertNextId(j);
					   bDeleteList[j] = true;
				   }


				   //std::cout<<"color: " << (int)color[0] << " " << (int)color[1] << " " << (int)color[2] << std::endl;

				   colors->InsertNextTuple3(color[0], color[1], color[2]);

			   }

		   }


       }
  }

  std::cout << "color # : " << colors->GetNumberOfTuples() << std::endl;

  output->GetPointData()->SetScalars(colors);
  
  
  //----------------------------------------------------------------------/
  // Parse the vertex, line and poly information
  //----------------------------------------------------------------------/
  //
  std::cout << std::endl << "------------- Vertex, Line and Poly ---------" << std::endl;
  std::cout << "Number of Cells : " << output->GetNumberOfCells() << std::endl;
  vtkSmartPointer<vtkCellArray> vert = output->GetVerts();
  std::cout << "verts # : " << vert->GetNumberOfCells() << std::endl;
  vtkSmartPointer<vtkCellArray> line = output->GetLines();
  std::cout << "lines # : " << line->GetNumberOfCells() << std::endl;
  vtkSmartPointer<vtkCellArray> poly = output->GetPolys();
  std::cout << "polys # : " << poly->GetNumberOfCells() << std::endl;
  
  for(int i=0; i<vert->GetNumberOfCells(); i++)
  {
	   vtkSmartPointer<vtkIdList> vertIdList = vtkSmartPointer<vtkIdList>::New();
	   vert->GetNextCell(vertIdList);
	   std::cout << "i = " << i << ", vert id # = " << vertIdList->GetNumberOfIds();
	   std::cout << " : ";
	   for(int j=0; j<vertIdList->GetNumberOfIds(); j++)
	   {
		   std::cout << " " << vertIdList->GetId(j);
	   }
	   std::cout << std::endl;

  }

  for(int i=0; i<line->GetNumberOfCells(); i++)
  {
	   vtkSmartPointer<vtkIdList> lineIdList = vtkSmartPointer<vtkIdList>::New();
	   line->GetNextCell(lineIdList);
	   std::cout << "i = " << i << ", line id # = " << lineIdList->GetNumberOfIds();
	   std::cout << " : ";
	   for(int j=0; j<lineIdList->GetNumberOfIds(); j++)
	   {
		   std::cout << " " << lineIdList->GetId(j);
	   }
	   std::cout << std::endl;
  }


  //---------------------------------------------------------//
  // Extract the polydata that satisfy the requirements.
  // 2 constraints: (1) gray color, and (2) normal information.
  //---------------------------------------------------------//

  vtkSmartPointer<vtkCellArray> newPoly = vtkSmartPointer<vtkCellArray>::New();


  vtkSmartPointer<vtkDataArray> normalArray = pd->GetArray(pd->GetArrayName(0));

  for(int i=0; i<poly->GetNumberOfCells(); i++)
  {
	   vtkSmartPointer<vtkIdList> polyIdList = vtkSmartPointer<vtkIdList>::New();
	   poly->GetNextCell(polyIdList);

	   vtkIdType id1 = polyIdList->GetId(0);
	   vtkIdType id2 = polyIdList->GetId(1);
	   vtkIdType id3 = polyIdList->GetId(2);

	   double avgnormal[3];
	   avgnormal[0] = (normalArray->GetComponent(id1,0) + normalArray->GetComponent(id2,0) + normalArray->GetComponent(id3,0))/3.0;
	   avgnormal[1] = (normalArray->GetComponent(id1,1) + normalArray->GetComponent(id2,1) + normalArray->GetComponent(id3,1))/3.0;
	   avgnormal[2] = (normalArray->GetComponent(id1,2) + normalArray->GetComponent(id2,2) + normalArray->GetComponent(id3,2))/3.0;

	   double dist = sqrt(avgnormal[0]*avgnormal[0] + avgnormal[1]*avgnormal[1] + avgnormal[2]*avgnormal[2]);
	   double crossVal = fabs(normal[0]*avgnormal[0] + normal[1]*avgnormal[1] + normal[2]*avgnormal[2])/dist;
	   

	   if(bDeleteList[id1] == false && bDeleteList[id2] == false && bDeleteList[id3] == false && crossVal < normalThreshold)
	   {
		   //std::cout << "Keep " << i << "-th polygon. " << std::endl;
		   newPoly->InsertNextCell(polyIdList);
	   }
  }

  output->SetPolys(newPoly);

  

  //-----------------------------------------------------------------//
  //  Save the results
  //-----------------------------------------------------------------//
  //
  // convert the unstructured grid into polydata
  vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
  geometryFilter->SetInputData(output);
  geometryFilter->Update(); 

  vtkXMLPolyDataWriter *writer=vtkXMLPolyDataWriter::New();
  writer->SetInputConnection(geometryFilter->GetOutputPort());
  writer->SetFileName(edit_filename.c_str());
  writer->SetDataModeToAscii();
  writer->Write();
  writer->Delete();

  
  //----------------------------------------------------------------// 
  // Visualize
  //----------------------------------------------------------------//
  
  
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInputConnection(output->GetProducerPort());
#else
  mapper->SetInputData(output);
#endif

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(.3, .6, .3); // Background color green

  renderWindow->Render();
  renderWindowInteractor->Start();


  return EXIT_SUCCESS;
}
