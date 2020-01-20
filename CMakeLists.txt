cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 14)


project(gray-scott C CXX)

set (DSINCLUDE /home/bz186/dataspaces/build_tcp/include)
set (DSLINKDIR /home/bz186/dataspaces/build_tcp/lib)
set (DSLINKLIB -ldspaces -ldscommon -ldart   -lm  -lrt)

set (ZFPINCLUDE /home/bz186/zfp/include)
set (ZFPLINKDIR /home/bz186/zfp/lib)
set (ZFPLINKLIB -lzfp -lm)

find_package(MPI REQUIRED)

option(USE_COMPRESSION "Use compression mode")
if(USE_COMPRESSION)
  message(STATUS "Enabling compression")
  add_definitions(-DENABLE_COMPRESSION)
endif()

option(USE_TIMERS "Use profiling timers")
if(USE_TIMERS)
  message(STATUS "Enabling profiling timers")
  add_definitions(-DENABLE_TIMERS)  
endif()





add_executable(gray-scott
  simulation/main.cpp
  simulation/gray-scott.cpp
  simulation/settings.cpp
)
target_include_directories(gray-scott PUBLIC ${gray-scott}/simulation ${DSINCLUDE} ${ZFPINCLUDE})
target_link_directories(gray-scott PUBLIC ${DSLINKDIR} ${ZFPLINKDIR})
target_link_libraries(gray-scott MPI::MPI_C ${DSLINKLIB} ${ZFPLINKLIB})


