cmake_minimum_required(VERSION 2.8)
 
project(GraphCuts)
 
set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

if(${VTK_FOUND})
add_definitions(-DVTK_FOUND)
endif()

set(TCLAP_INCLUDE_DIR ~/code/tclap-1.2.1/include/ CACHE PATH "Path to tclap/ where CmdLine.h is")

include_directories( BEFORE ${TCLAP_INCLUDE_DIR} )

set(MAXFLOW_SRC maxflow/maxflow.cpp)

set(MAXFLOW_INCLUDE_DIRS maxflow)

# we need this so headers show up on qt creator
set(GRAPHCUTS_HEADERS graphCutsAdapter.h
                      graphCutsAdapter.hxx
                      verbosityConstant.h)

# Add this if we wan to provide a library of the itk wrapper + maxflow.
# The itk wrapper is header only, so no need for now
#add_library(graphCuts ${GRAPHCUTS_HEADERS})
#set_target_properties(graphCuts PROPERTIES LINKER_LANGUAGE CXX)

add_executable(graphCutter graphCutter.cpp ${MAXFLOW_SRC} ${GRAPHCUTS_HEADERS})
target_link_libraries(graphCutter
                      #graphCuts
                      ${Glue} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
