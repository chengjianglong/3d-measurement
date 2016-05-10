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
#include <vtkXMLPolyDataWriter.h>
#include <vtkCenterOfMass.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPerspectiveTransform.h>
#include <vtkPLYWriter.h>
#include <vtkAxesActor.h>

#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkBoundingBox.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <vtkVersion.h>
#include <vtkSphereSource.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOutlineFilter.h>

int main(int argc, char *argv[])
{
  // Parse command line arguments
  if(argc != 5)
  {
    std::cerr << "Usage: " << argv[0]
              << " Filename(.vtp) Fitting_plane_filename(.txt) Fitting_line_filename(.txt) method" << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];
  std::string plane_filename = argv[2];
  std::string line_filename = argv[3];
  int method = atoi(argv[4]);

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

  //-------------------------------------------------------------//
  // Read the parameter for the fitting line.
  //-------------------------------------------------------------//
  std::vector< std::vector<double> > mylines;
  std::ifstream line_ifs(line_filename.c_str(), std::ifstream::in);
  if(line_ifs.is_open())
  {
      int lno = 0;
	  int vno;
	  while(line_ifs.good())
	  {
		  lno++;
		  std::getline(line_ifs, strline);
		  
		  if(lno == 2 || lno == 5 || lno == 8)
		  {
              //std::cout << lno << ": " << strline << std::endl;
			  std::stringstream ss(strline);
			  std::string token;
			  std::vector<std::string> valstrs;
			  while(std::getline(ss, token, ' '))
			  {
				  valstrs.push_back(token);
				  //std::cout << token << std::endl;
			  }

			  std::vector<double> params;
			  for(int i=0; i<3; i++)
			  {
				  params.push_back(stod(valstrs[i]));
			  }

			  mylines.push_back(params);
		  }
	  }
  }

  line_ifs.close();


  double lineParams[3];
  for(int i=0; i<myplanes.size(); i++)
  {
      std::cout << "Plane normal = " << i << ": ";
	  for(int j=0; j<myplanes[i].size(); j++)
	  {
	      std::cout << " " << mylines[i][j];
		  
		  if(i==method && j<3)
		  {
			  lineParams[j] = mylines[i][j];
		  }
	  }
	  std::cout <<  std::endl;
  }


  //-------------------------------------------------------------//
  // Read the point clouds from the file.
  //-------------------------------------------------------------//

  // Read all the data from the file
  vtkSmartPointer<vtkXMLPolyDataReader> reader =  vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter = vtkSmartPointer<vtkCenterOfMass>::New();
  centerOfMassFilter->SetInputConnection(reader->GetOutputPort());
  centerOfMassFilter->SetUseScalarsAsWeights(false);
  centerOfMassFilter->Update();

  double center[3];
  centerOfMassFilter->GetCenter(center);
  std::cout << "Center:  ( " << center[0] << " , " << center[1] << " , " << center[2] << " ) " << std::endl;

  double zNormal[3]; // = {-0.156764, -0.0543708, 0.986138};
  double xPara[3]; // = {-0.0564942, 0.158968, 0.0551351};

  zNormal[0] = planeNormal[0];
  zNormal[1] = planeNormal[1];
  zNormal[2] = planeNormal[2];
  xPara[0] = lineParams[0];
  xPara[1] = lineParams[1];
  xPara[2] = lineParams[2];

  double xNormal[3];
  
  xNormal[0] = -xPara[2];
  xNormal[1] = xPara[0]+xPara[1];
  xNormal[2] = xPara[0]*xPara[2];

  double xconst = sqrt(xNormal[0]*xNormal[0] + xNormal[1]*xNormal[1] + xNormal[2]*xNormal[2]);
  xNormal[0] /= xconst;
  xNormal[1] /= xconst;
  xNormal[2] /= xconst;

  std::cout << "xconst = " << xconst << std::endl;
  

  double yNormal[3];
  vtkMath::Cross(xNormal, zNormal, yNormal);

  double yconst = sqrt(yNormal[0]*yNormal[0] + yNormal[1]*yNormal[1] + yNormal[2]*yNormal[2]);
  yNormal[0] /= yconst;
  yNormal[1] /= yconst;
  yNormal[2] /= yconst;

  std::cout << "yconst = " << yconst << std::endl;
  
  
  std::cout << "X-axis:  ( " << xNormal[0] << " , " << xNormal[1] << " , " << xNormal[2] << " ) " << std::endl;
  std::cout << "Y-axis:  ( " << yNormal[0] << " , " << yNormal[1] << " , " << yNormal[2] << " ) " << std::endl;
  std::cout << "Z-axis:  ( " << zNormal[0] << " , " << zNormal[1] << " , " << zNormal[2] << " ) " << std::endl;

  //check the axis
  double checkxy = xNormal[0]*yNormal[0] + xNormal[1]*yNormal[1] + xNormal[2]*yNormal[2];
  double checkxz = xNormal[0]*zNormal[0] + xNormal[1]*zNormal[1] + xNormal[2]*zNormal[2];
  double checkyz = yNormal[0]*zNormal[0] + yNormal[1]*zNormal[1] + yNormal[2]*zNormal[2];
  std::cout << "Check the direction by CheckXY = " << checkxy << std::endl;
  std::cout << "Check the direction by CheckXZ = " << checkxz << std::endl;
  std::cout << "Check the direction by CheckYX = " << checkyz << std::endl;

  double px[3];
  px[0] = center[0]+xNormal[0];
  px[1] = center[1]+xNormal[1];
  px[2] = center[2]+xNormal[2];

  double py[3];
  py[0] = center[0]+yNormal[0];
  py[1] = center[1]+yNormal[1];
  py[2] = center[2]+yNormal[2];

  double pz[3];
  pz[0] = center[0]+zNormal[0];
  pz[1] = center[1]+zNormal[1];
  pz[2] = center[2]+zNormal[2];

  //---------------------------------------------------------------//
  // Visualize the directions for (xNormal, yNormal, zNormal).
  //---------------------------------------------------------------//

  // Create the polydata where we will store all the geometric data
  vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();

  vtkSmartPointer<vtkPoints> pts =  vtkSmartPointer<vtkPoints>::New();
  pts->InsertNextPoint(center);
  pts->InsertNextPoint(px);
  pts->InsertNextPoint(py);
  pts->InsertNextPoint(pz);

  // Add the points to the polydata container
  linesPolyData->SetPoints(pts);

  // Create the first line (between Origin and P0)
  vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
  line0->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
  line0->GetPointIds()->SetId(1, 1); // the second 1 is the index of P0 in linesPolyData's points

  vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
  line1->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
  line1->GetPointIds()->SetId(1, 2); // the second 1 is the index of P0 in linesPolyData's points

  vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
  line2->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
  line2->GetPointIds()->SetId(1, 3); // the second 1 is the index of P0 in linesPolyData's points

  // Create a vtkCellArray container and store the lines in it
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  lines->InsertNextCell(line0);
  lines->InsertNextCell(line1);
  lines->InsertNextCell(line2);

  // Add the lines to the polydata container
  linesPolyData->SetLines(lines);
  
  // Create two colors - one for each line
  unsigned char red[3] = { 255, 0, 0 };
  unsigned char green[3] = { 0, 255, 0 };
  unsigned char blue[3] = { 0, 0, 255 };

  // Create a vtkUnsignedCharArray container and store the colors in it
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
#if VTK_MAJOR_VERSION < 7
  colors->InsertNextTupleValue(red);
  colors->InsertNextTupleValue(green);
  colors->InsertNextTupleValue(blue);
#else
  colors->InsertNextTypedTuple(red);
  colors->InsertNextTypedTuple(green);
  colors->InsertNextTypedTuple(blue);
#endif
                          
  // Color the lines.
  // SetScalars() automatically associates the values in the data array passed as 
  // to the elements in the same indices of the cell data array on which it is called.
  // This means the first component (red) of the colors array
  // is matched with the first component of the cell array (line 0)
  // and the second component (green) of the colors array
  // is matched with the second component of the cell array (line 1)
  linesPolyData->GetCellData()->SetScalars(colors);
  

  // Setup the visualization pipeline
  vtkSmartPointer<vtkPolyDataMapper> linemapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
  linemapper->SetInputData(linesPolyData);

  vtkSmartPointer<vtkActor> lineactor = vtkSmartPointer<vtkActor>::New();
  lineactor->SetMapper(linemapper);
  
  //----------------------------------------------------------------//
  // Coordinate transformation to extract bounding box.
  //----------------------------------------------------------------//

  vtkSmartPointer<vtkMatrix4x4> RTMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  RTMatrix->SetElement(0, 0, xNormal[0]);
  RTMatrix->SetElement(1, 0, yNormal[0]);
  RTMatrix->SetElement(2, 0, zNormal[0]);
  RTMatrix->SetElement(3, 0, 0.0);
  RTMatrix->SetElement(0, 1, xNormal[1]);
  RTMatrix->SetElement(1, 1, yNormal[1]);
  RTMatrix->SetElement(2, 1, zNormal[1]);
  RTMatrix->SetElement(3, 1, 0.0);
  RTMatrix->SetElement(0, 2, xNormal[2]);
  RTMatrix->SetElement(1, 2, yNormal[2]);
  RTMatrix->SetElement(2, 2, zNormal[2]);
  RTMatrix->SetElement(3, 2, 0.0);

  RTMatrix->SetElement(0, 3, -xNormal[0]*center[0] - xNormal[1]*center[1] - xNormal[0]*center[2]);
  RTMatrix->SetElement(1, 3, -yNormal[0]*center[0] - yNormal[1]*center[1] - yNormal[1]*center[2]);
  RTMatrix->SetElement(2, 3, -zNormal[0]*center[0] - zNormal[1]*center[1] - zNormal[2]*center[2]);
  RTMatrix->SetElement(3, 3, 1.0);
  

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  double orient[3];
  transform->GetOrientation(orient, RTMatrix);
  std::cout << "Orientation:  " << orient[0] << " , " << orient[1] << " , " << orient[2] << std::endl;
  transform->SetMatrix(RTMatrix);
  transform->Update();


  
  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetTransform(transform);
  transformFilter->SetInputConnection(reader->GetOutputPort());
  transformFilter->Update();


  //-------------------------------------------------------//
  // Save the transformed point cloud (.ply) for verification.
  //-------------------------------------------------------//

  vtkSmartPointer<vtkPLYWriter> plywriter = vtkSmartPointer<vtkPLYWriter>::New();
  plywriter->SetInputConnection(transformFilter->GetOutputPort());
  plywriter->SetFileTypeToASCII();
  plywriter->SetFileName("transformed.ply");
  plywriter->Write();
 

  
  vtkSmartPointer<vtkPerspectiveTransform> perspectiveTransform = vtkSmartPointer<vtkPerspectiveTransform>::New();
  perspectiveTransform->SetMatrix(RTMatrix);
  perspectiveTransform->Update();
  
  
  //-------------------------------------------------------//
  //  Bounding box
  // ------------------------------------------------------//

  vtkBoundingBox boundingBox;

  std::cout << std::endl << "------------- Point ---------" << std::endl;
  //vtkSmartPointer<vtkPolyData> output =  reader->GetOutput();
  vtkSmartPointer<vtkPolyData> output =  transformFilter->GetOutput();
  std::cout << "point # : " << output->GetNumberOfPoints() << std::endl;
  for(int i=0; i<output->GetNumberOfPoints(); i++)
  {
	  double p[3];
	  output->GetPoint(i,p);
	 
	  boundingBox.AddPoint(p);

	  
	  /*
	  double p[3];
	  output->GetPoint(i,p);
	  std::cout << "P_" << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ") -> ";
	  
	  double dist = (p[0]-center[0])*zNormal[0] + (p[1]-center[1])*zNormal[1] + (p[2]-center[2])*zNormal[2];
	  std::cout << " dist = " << dist << "  ===  ";
	  //double perspectiveProjection[4];
	  double p4[4];
	  p4[0] = p[0];
	  p4[1] = p[1];
	  p4[2] = p[2];
	  p4[3] = 1.0;
	  double *perspectiveProjection = RTMatrix->MultiplyDoublePoint(p4);
	  //perspectiveTransform->TransformPoint(p, perspectiveProjection);
	  std::cout << "(" << perspectiveProjection[0] << " " << perspectiveProjection[1] << " "  << perspectiveProjection[2]  << ")" << std::endl;
	  */
  }

  //-------------------------------------------------------//
  // Extract the bounding box.
  //-------------------------------------------------------//

  double bounds[6]; 
  boundingBox.GetBounds(bounds);
  std::cout << "Bounds: " << bounds[0] << ", " << bounds[1] << ", " << bounds[2] << ", "
	        << bounds[3] << ", " << bounds[4] << ", " << bounds[5] << std::endl;

  double width = bounds[1] - bounds[0];
  double length = bounds[3] - bounds[2];
  if(width > length)
  {
	  double tmp = width;
	  width = length;
	  length = tmp;
  }

  std::cout << "Measurement (foot length) : " << length << std::endl;
  std::cout << "Measurement (foot width ) : " << width << std::endl;



  //-------------------------------------------------------//
  // Visualization.
  //-------------------------------------------------------//

  // Setup the text and add it to the renderer
  vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
  textActor->SetInput ( "Hello world" );
  textActor->SetPosition2 (400, 400 );
  textActor->GetTextProperty()->SetFontSize ( 24 );
  textActor->GetTextProperty()->SetColor ( 1.0, 0.0, 0.0 );

  

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  //mapper->SetInputConnection(reader->GetOutputPort());
  mapper->SetInputConnection(transformFilter->GetOutputPort());
  vtkSmartPointer<vtkActor> actor =  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  //actor->SetUserTransform(transform);

  // Create the outline
  vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();
  //outline->SetInputConnection(reader->GetOutputPort());
  outline->SetInputConnection(transformFilter->GetOutputPort());
  
  vtkSmartPointer<vtkPolyDataMapper> outlineMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
  outlineMapper->SetInputConnection(outline->GetOutputPort());
  vtkSmartPointer<vtkActor> outlineActor =  vtkSmartPointer<vtkActor>::New();
  outlineActor->SetMapper(outlineMapper);
  outlineActor->GetProperty()->SetColor(0,0,0);
  //outlineActor->SetUserTransform(transform);
  
  // Setup the window
  vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);



  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
  axes->AxisLabelsOff();
  axes->SetTotalLength(1,1,1);

  // Add the actors to the scene
  renderer->AddActor(axes);
  renderer->AddActor(actor);
  //renderer->AddActor(lineactor);
  renderer->AddActor(outlineActor);
  renderer->AddActor2D ( textActor );
  renderer->SetBackground(0.3, 0.6, 0.3); // Background color white

  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
