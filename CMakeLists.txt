
cmake_minimum_required(VERSION 2.8)

ADD_DEFINITIONS(
	-std=c++11 # Or -std=c++0x
	# Other flags
	)

PROJECT(ReadPolyData)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(ReadPolyData MACOSX_BUNDLE ReadPolyData)
add_executable(VisualizePolyData MACOSX_BUNDLE VisualizePolyData)
add_executable(SmoothPolyData MACOSX_BUNDLE SmoothPolyData)
add_executable(ContoursFromPolyData MACOSX_BUNDLE ContoursFromPolyData)
add_executable(EditPolyData MACOSX_BUNDLE EditPolyData)
add_executable(CutPolyData MACOSX_BUNDLE CutPolyData)
add_executable(VisualizeFitPlane MACOSX_BUNDLE VisualizeFitPlane)
add_executable(ExtractObject MACOSX_BUNDLE ExtractObject)
add_executable(Measure MACOSX_BUNDLE Measure)

if(VTK_LIBRARIES)
  target_link_libraries(ReadPolyData ${VTK_LIBRARIES})
else()
  target_link_libraries(ReadPolyData vtkHybrid vtkWidgets)
endif()

if(VTK_LIBRARIES)
  target_link_libraries(VisualizePolyData ${VTK_LIBRARIES})
else()
  target_link_libraries(VisualizePolyData vtkHybrid vtkWidgets)
endif()

if(VTK_LIBRARIES)
  target_link_libraries(SmoothPolyData ${VTK_LIBRARIES})
else()
  target_link_libraries(SmoothPolyData vtkHybrid vtkWidgets)
endif()

if(VTK_LIBRARIES)
  target_link_libraries(ContoursFromPolyData ${VTK_LIBRARIES})
else()
  target_link_libraries(ContoursFromPolyData vtkHybrid vtkWidgets)
endif()

if(VTK_LIBRARIES)
  target_link_libraries(EditPolyData ${VTK_LIBRARIES})
else()
  target_link_libraries(EditPolyData vtkHybrid vtkWidgets)
endif()

if(VTK_LIBRARIES)
	target_link_libraries(CutPolyData ${VTK_LIBRARIES})
else()
	target_link_libraries(CutPolyData vtkHybrid vtkWidgets)
endif()

if(VTK_LIBRARIES)
	target_link_libraries(VisualizeFitPlane ${VTK_LIBRARIES})
else()
	target_link_libraries(VisualizeFitPlane vtkHybrid vtkWidgets)
endif()

if(VTK_LIBRARIES)
	target_link_libraries(ExtractObject ${VTK_LIBRARIES})
else()
	target_link_libraries(ExtractObject vtkHybrid vtkWidgets)
endif()

if(VTK_LIBRARIES)
	target_link_libraries(Measure ${VTK_LIBRARIES})
else()
	target_link_libraries(Measure vtkHybrid vtkWidgets)
endif()
