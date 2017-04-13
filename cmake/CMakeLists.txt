cmake_minimum_required(VERSION 3.6)
project(MOM5 C Fortran)

find_package(MPI)

if ( NOT MPI_FOUND )
  set(MPI_FORTRAN_COMPILER mpifort)
  set(MPI_C_COMPILER mpicc)
  set(CMAKE_Fortran_COMPILER mpifort)
  set(CMAKE_C_COMPILER mpicc)
endif ( NOT MPI_FOUND )

set(CMAKE_Fortran_COMPILER mpifort)
set(CMAKE_C_COMPILER mpicc)

execute_process(COMMAND nf-config --fflags OUTPUT_VARIABLE NETCDF_FFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND nf-config --flibs OUTPUT_VARIABLE NETCDF_FLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

message("Using Fortran: ${CMAKE_Fortran_COMPILER} and C: ${CMAKE_C_COMPILER}")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(FFLAGS "-fno-alias -stack_temps -safe-cray-ptr -ftz -i_dynamic -assume byterecl -i4 -r8 -nowarn -sox -traceback")
  set(FFLAGS_RELEASE "${FFLAGS} -O2 -fp-model source")
  set(FFLAGS_DEBUG "${FFLAGS} -g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(FFLAGS "-fcray-pointer -fdefault-real-8 -fdefault-double-8 -Waliasing -ffree-line-length-none -fno-range-check")
  set(FFLAGS_RELEASE "${FFLAGS} -O2 -fbounds-check")
  set(FFLAGS_DEBUG "${FFLAGS} -O0 -g -W -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${FFLAGS_DEBUG} ${NETCDF_FFLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${FFLAGS_RELEASE} ${NETCDF_FFLAGS}")

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -sox -traceback -O2 -debug minimal")
elseif(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
endif()

message("Using Fortran flags: ${CMAKE_Fortran_FLAGS}")
message("NetCDF libs: ${NETCDF_FLIBS}")

add_definitions(-Duse_libMPI -Duse_netCDF -Duse_LARGEFILE -DSPMD -D__IFC)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/../src")

add_subdirectory(FMS)
add_subdirectory(ocean)
add_subdirectory(MOM5_solo)
add_subdirectory(ice_sis)
add_subdirectory(MOM5_SIS)
