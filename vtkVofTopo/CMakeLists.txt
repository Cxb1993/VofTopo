cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11 -fdiagnostics-color")

if (ParaView_SOURCE_DIR)
  include_directories(${VTK_INCLUDE_DIRS})
else (ParaView_SOURCE_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif (ParaView_SOURCE_DIR)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

FIND_PATH(INTERVOL_DIR intervol_app_common.h)
INCLUDE_DIRECTORIES(${INTERVOL_DIR})

add_library(vofTopology vofTopology.cxx)
add_library(marchingCubes_cpu marchingCubes_cpu.cxx)
add_library(correspondence correspondence.cxx)

add_paraview_plugin(VofTopo "1.0"
  SERVER_MANAGER_XML VofTopo.xml
  SERVER_MANAGER_SOURCES vtkVofTopo.cxx)

target_link_libraries(VofTopo PRIVATE
  vofTopology
  marchingCubes_cpu
  correspondence)
