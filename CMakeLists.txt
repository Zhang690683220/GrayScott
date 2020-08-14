cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 14)


project(gray-scott C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

find_package(MPI REQUIRED)

find_package(DataSpaces REQUIRED)


option(USE_TIMERS "Use profiling timers")
if(USE_TIMERS)
  message(STATUS "Enabling profiling timers")
  add_definitions(-DENABLE_TIMERS)
endif()

# We are not using the C++ API of MPI, this will stop the compiler look for it
add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

add_executable(gray-scott
  simulation/main.cpp
  simulation/gray-scott.cpp
  simulation/settings.cpp
  simulation/writer.cpp
)
target_link_libraries(gray-scott DataSpaces::DataSpaces MPI::MPI_C)

install(FILES
          ${CMAKE_MODULE_PATH}/FindDataSpaces.cmake
          DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake/Modules )

  # Need to install/export Dataspaces find module for downstream dependencies
  configure_file(cmake/Modules/FindDataSpaces.cmake Modules/FindDataSpaces.cmake COPYONLY)

option(VTK "Build VTK apps")
if (VTK_ROOT)
  set(VTK ON)
endif(VTK_ROOT)
  
if (VTK)
  message(STATUS "Configuring VTK apps")
  
  find_package(VTK COMPONENTS
    vtkCommonCore
    vtkFiltersCore
    vtkFiltersGeneral
    vtkIOImage
    vtkIOXML
  )

  if(VTK_FOUND)
    include(${VTK_USE_FILE})
    add_executable(isosurface_ds
      analysis/isosurface_ds.cpp
      analysis/analysis.cpp
      analysis/reader.cpp)
    target_link_libraries(isosurface_ds DataSpaces::DataSpaces ${VTK_LIBRARIES}
      MPI::MPI_C)
  endif(VTK_FOUND)


  endif(VTK)



