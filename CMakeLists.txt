cmake_minimum_required(VERSION 2.6)
PROJECT(qmc_master)
INCLUDE(CheckFunctionExists)

find_package(MPI)
Find_Package(SDL)

SET(CMAKE_CXX_FLAGS_RELEASE "-march=corei7 -O3 -mmmx -finline-functions -fopenmp -Wconversion")
SET(CMAKE_C_FLAGS_RELEASE "-march=corei7 -O3 -mmmx -finline-functions -fopenmp -Wconversion")
SET(CMAKE_CXX_FLAGS_DEBUG "-march=corei7 -gstabs+ -Wall -Wabi -Wconversion")
SET(CMAKE_C_FLAGS_DEBUG "-march=corei7 -O3 -mmmx -finline-functions -Wconversion -gstabs+")
SET(CMAKE_CXX_FLAGS_PROFILED "-O2 -march=corei7 -pg -Wconversion")
SET(CMAKE_C_FLAGS_PROFILED "-O3 -march=corei7 -pg -Wconversion")

include_directories(common/ SDL_QMC_Master/ MPI_QMC_Master)


Set(MPI_Master_src ./MPI_QMC_Master/qmc_master.cpp ./MPI_QMC_Master/registry.cpp)

Set(SDL_Master_src ./SDL_QMC_Master/qmc_master.cpp ./SDL_QMC_Master/stripped_SDL_Net.c ./SDL_QMC_Master/registry.cpp)

ADD_EXECUTABLE(MPI_QMC_Master ${MPI_Master_src})

ADD_EXECUTABLE(SDL_QMC_Master ${SDL_Master_src})

SET(MPI_Master_Defines "-DUSE_MPI -I. -DHAVE_CONFIG_H")

#let's build the proper include file string from MPI_INCLUDE_PATH
foreach(elem ${MPI_INCLUDE_PATH})
set(var "-I${elem}")
set(tempflags "${tempflags} ${var}")
endforeach(elem)

set(MPI_files_flags "-DUSE_MPI ${tempflags}")

#let's retrieve the mpi library string...

set(templinkflags "")
foreach(elem ${MPI_EXTRA_LIBRARY})
get_filename_component(fln ${elem} NAME_WE)#strip path
string(REPLACE "lib" "" fln ${fln}) # strip the lib prefix of the libraries
set(templinkflags "${templinkflags} -l${fln}")
endforeach(elem)

IF(CMAKE_COMPILER_IS_GNUCXX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

set_source_files_properties(${MPI_Master_src} PROPERTIES COMPILE_FLAGS "${MPI_files_flags}")

set_target_properties(MPI_QMC_Master PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/MPI_QMC_Master COMPILE_FLAGS "${MPI_Master_Defines}" LINK_FLAGS "${MPI_LINK_FLAGS} ${templinkflags}")

set_target_properties(SDL_QMC_Master PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SDL_QMC_Master COMPILE_FLAGS " ${COMPILE_FLAGS} -DUSE_SDL -I. -DHAVE_CONFIG_H" )

target_link_libraries(SDL_QMC_Master ${SDL_LIBRARY})
target_link_libraries(MPI_QMC_Master ${MPI_LIBRARIES}) 

CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)

CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)

CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)