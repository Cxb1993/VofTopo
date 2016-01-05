cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

if (ParaView_SOURCE_DIR)
  include_directories(${VTK_INCLUDE_DIRS})
else (ParaView_SOURCE_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif (ParaView_SOURCE_DIR)

find_package(CUDA REQUIRED)
include_directories(${CUDA_SDK_ROOT_DIR}/common/inc)
include_directories(${CUDA_INCLUDE_DIRS})

add_library(vofTopology vofTopology.cxx)
add_library(marchingCubes_cpu marchingCubes_cpu.cxx)

add_paraview_plugin(VofTopo "1.0"
  SERVER_MANAGER_XML VofTopo.xml
  SERVER_MANAGER_SOURCES vtkVofTopo.cxx)

target_link_libraries(VofTopo PRIVATE vofTopology marchingCubes_cpu)
