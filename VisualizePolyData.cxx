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
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkProperty.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

int main ( int argc, char *argv[] )
{
  // Parse command line arguments
  if(argc != 2)
    {
    std::cerr << "Usage: " << argv[0]
              << " Filename(.vtp)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string filename = argv[1];

  // Read all the data from the file
  vtkSmartPointer<vtkXMLPolyDataReader> reader =  vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  //------------------------------------------------------------------//
  // Extract the color information.
  // -----------------------------------------------------------------//
  //

  // Extract the polydata
  vtkSmartPointer<vtkPolyData> output = reader->GetOutput();
  
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
  
  
  //---------------------------------------------------------------------//
  // Visualization
  //---------------------------------------------------------------------//
 
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInputConnection(output->GetProducerPort());
#else
  mapper->SetInputData(output);
#endif

  vtkSmartPointer<vtkActor> actor =  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =  vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(.3, .6, .3); // Background color green

  renderWindow->Render();
  renderWindowInteractor->Start();


  return EXIT_SUCCESS;
}
