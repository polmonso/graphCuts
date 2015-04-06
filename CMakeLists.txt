cmake_minimum_required(VERSION 2.8)
 
project(GraphCuts)
 
set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()
 
include_directories( BEFORE ~/code/tclap-1.2.1/include )

set(MAXFLOW_SRC )
set(MAXFLOW_HDR )

set(MAXFLOW_INCLUDE_DIRS maxflow)

# we need this so headers show up on qt creator
set(GRAPHCUTS_HEADERS graphCutsAdapter.h graphCutsAdapter.hxx VerbosityConstant.h visualization.h)

# Add this if we wan to provide a library of the itk wrapper + maxflow.
# The itk wrapper is header only, so no need for now
#add_library(graphCuts ${GRAPHCUTS_HEADERS})
#set_target_properties(graphCuts PROPERTIES LINKER_LANGUAGE CXX)

add_executable(graphCutter graphCutter.cpp ${GRAPHCUTS_HEADERS})
target_link_libraries(graphCutter
                      #graphCuts
                      ${Glue} ${VTK_LIBRARIES} ${ITK_LIBRARIES})
