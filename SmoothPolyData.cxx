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

#include <vtkXMLPolyDataReader.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkDelaunay2D.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkCamera.h>


int main(int argc, char *argv[])
{
    // Parse command line arguments
    if(argc != 3)
    {
    std::cerr << "Usage: " << argv[0]
              << " input_filename(.vtp)  output_filename(.vtp)" << std::endl;
    return EXIT_FAILURE;
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];

	//---------------------------------------------------//
    // Read all the data from the file
	//---------------------------------------------------//
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(input_filename.c_str());
    reader->Update();


	//---------------------------------------------------//
    // Smooth the polydata
	//---------------------------------------------------//
	//
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputConnection(reader->GetOutputPort());
    smoothFilter->SetNumberOfIterations(15);
    //smoothFilter->SetRelaxationFactor(0.1);
    smoothFilter->SetRelaxationFactor(1.0);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();

    // Update normals on newly smoothed polydata
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputConnection(smoothFilter->GetOutputPort());
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOn();
    normalGenerator->Update();

    vtkSmartPointer<vtkPolyDataMapper> inputMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    inputMapper->SetInputConnection(reader->GetOutputPort());
    vtkSmartPointer<vtkActor> inputActor = vtkSmartPointer<vtkActor>::New();
    inputActor->SetMapper(inputMapper);

    vtkSmartPointer<vtkPolyDataMapper> smoothedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    smoothedMapper->SetInputConnection(normalGenerator->GetOutputPort());
    vtkSmartPointer<vtkActor> smoothedActor = vtkSmartPointer<vtkActor>::New();
    smoothedActor->SetMapper(smoothedMapper);
  
	
    //-----------------------------------------------------------------//
    //  Save the results
    //-----------------------------------------------------------------//
    //
    // convert the unstructured grid into polydata
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputConnection(normalGenerator->GetOutputPort());
    geometryFilter->Update(); 

    vtkXMLPolyDataWriter *writer=vtkXMLPolyDataWriter::New();
    writer->SetInputConnection(geometryFilter->GetOutputPort());
    writer->SetFileName(output_filename.c_str());
    writer->SetDataModeToAscii();
    writer->Write();
    writer->Delete();



    //-----------------------------------------------------------------//
    // Visualization.
    //-----------------------------------------------------------------//

    // There will be one render window
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(600, 300);

    // And one interactor
    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    // Define viewport ranges
    // (xmin, ymin, xmax, ymax)
    double leftViewport[4] = { 0.0, 0.0, 0.5, 1.0 };
    double rightViewport[4] = { 0.5, 0.0, 1.0, 1.0 };

    // Setup both renderers
    vtkSmartPointer<vtkRenderer> leftRenderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow->AddRenderer(leftRenderer);
    leftRenderer->SetViewport(leftViewport);
    leftRenderer->SetBackground(.6, .5, .4);

    vtkSmartPointer<vtkRenderer> rightRenderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow->AddRenderer(rightRenderer);
    rightRenderer->SetViewport(rightViewport);
    rightRenderer->SetBackground(.4, .5, .6);

    // Add the input parabola to the left and the smoothed parabola to the right
    leftRenderer->AddActor(inputActor);
    rightRenderer->AddActor(smoothedActor);

    leftRenderer->ResetCamera();
    leftRenderer->GetActiveCamera()->Azimuth(130);
    leftRenderer->GetActiveCamera()->Elevation(-80);

    rightRenderer->ResetCamera();
    rightRenderer->GetActiveCamera()->Azimuth(130);
    rightRenderer->GetActiveCamera()->Elevation(-80);

    renderWindow->Render();
    interactor->Start();

    return EXIT_SUCCESS;
}
