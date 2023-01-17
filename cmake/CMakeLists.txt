cmake_minimum_required(VERSION 3.6)
project(MOM5 C Fortran)

set(CMAKE_VERBOSE_MAKEFILE ON)

# make sure that the default is RELWITHDEBINFO
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELWITHDEBINFO CACHE STRING
      "Choose the type of build, options are: Debug Release Relwithdebinfo."
      FORCE)
endif()
message("Build type: " ${CMAKE_BUILD_TYPE})

# If your system uses a bespoke MPI compiler wrapper then define this on the
# command line using the MPIFORT_EXE and MPICC_EXE environment variables
if ((DEFINED ENV{MPIFORT_EXE}) AND (DEFINED ENV{MPICC_EXE}))
  message("mpifort executable defined: " $ENV{MPIFORT_EXE})
  message("Will assume system MPI implementation is sound. Remove MPIFORT_EXE from environment to automatically configure MPI")
  set(CMAKE_Fortran_COMPILER $ENV{MPIFORT_EXE})
  message("mpicc executable found: " $ENV{MPICC_EXE})
  message("Will assume system MPI implementation is sound. Remove MPICC_EXE from environment to automatically configure MPI")
  set(CMAKE_C_COMPILER $ENV{MPICC_EXE})
else ()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND (ENV{FC} STREQUAL "gfortran" AND ENV{CC} STREQUAL "gcc"))
    # Use compiler wrappers if FC and CC are default values 
    SET( CMAKE_Fortran_COMPILER mpif90 )
    SET( CMAKE_C_COMPILER mpicc )
  endif()
endif()
find_package(MPI REQUIRED COMPONENTS Fortran)

message("Using Fortran: ${MPI_Fortran_COMPILER} ${CMAKE_Fortran_COMPILER} and C: ${MPI_C_COMPILER} ${CMAKE_C_COMPILER}")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS                "-fno-alias -stack-temps -safe-cray-ptr -ftz -shared-intel -assume byterecl -i4 -r8 -nowarn -sox -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE        "-O2 -fp-model source")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")
  set(CMAKE_Fortran_FLAGS_DEBUG          "-g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv") 
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  include(CheckFortranCompilerFlag)
  set(CMAKE_Fortran_FLAGS                "-fcray-pointer -fdefault-real-8 -fdefault-double-8 -Waliasing -ffree-line-length-none -fno-range-check -fbacktrace")
  check_fortran_compiler_flag("-fallow-invalid-boz" _boz_flag)
  check_fortran_compiler_flag("-fallow-argument-mismatch" _argmis_flag)
  if(_boz_flag)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-invalid-boz" )
  endif()
  if(_argmis_flag)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch" )
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE        "-O2 -fbounds-check")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")
  set(CMAKE_Fortran_FLAGS_DEBUG          "-O0 -g -W -fbounds-check -ffpe-trap=invalid,zero,overflow")
else ()
  message ("Unknown Fortran compiler default flags only...")
endif()

execute_process(COMMAND nf-config --includedir OUTPUT_VARIABLE NETCDF_INCDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND nf-config --flibs OUTPUT_VARIABLE NETCDF_FLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${NETCDF_INCDIR}")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${NETCDF_FLIBS}")

message("Linker flags: ${CMAKE_EXE_LINKER_FLAGS}")

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS                "-sox -traceback") 
  set(CMAKE_C_FLAGS_DEBUG          "-O0 -debug")
  set(CMAKE_C_FLAGS_RELEASE        "-O2 -debug minimal")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELEASE} -g")
  add_compile_definitions(__IFC)
elseif(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS                "-frecord-gcc-switches")
  set(CMAKE_C_FLAGS_DEBUG          "-O0 -Werror -Wuninitialized -Wno-stringop-overflow")
  set(CMAKE_C_FLAGS_RELEASE        "-O2")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELEASE} -g")
else ()
  message ("Unknown C compiler default flags only...")
endif()
execute_process(COMMAND nc-config --cflags OUTPUT_VARIABLE NETCDF_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${NETCDF_CFLAGS}")

# Transfer command line setting of CPP defs to CMake 
set(CPPDEFS 
  use_netCDF
  use_netCDF3
  use_netCDF4
  use_libMPI
  SPMD
  ENABLE
  ODA
  USE_OCEAN_BGC
  LAND_BND_TRACERS 
  OVERLOAD_C8 
  OVERLOAD_C4 
  OVERLOAD_R4
)
foreach(VAR ${CPPDEFS})
  if(${VAR})
    message("Setting CPP def ${VAR}")
    add_compile_definitions(${VAR})
  endif()
endforeach()

string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
message("Using Fortran flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}}")

set(CMAKE_DEBUG_TARGET_PROPERTIES
            # INCLUDE_DIRECTORIES
            # COMPILE_DEFINITIONS
            # POSITION_INDEPENDENT_CODE
            # CONTAINER_SIZE_REQUIRED
            # LIB_VERSION
          )

set(SRC_DIR "${CMAKE_SOURCE_DIR}/../src")

add_subdirectory(version)
add_subdirectory(FMS)
add_subdirectory(ocean)
add_subdirectory(atmos_bg)
add_subdirectory(atmos_ebm)
add_subdirectory(atmos_fv)
add_subdirectory(atmos_null)
add_subdirectory(atmos_phys)
add_subdirectory(land_null)
add_subdirectory(land_lad)
add_subdirectory(land_lad2)
add_subdirectory(ice_sis)
add_subdirectory(MOM5_solo)
add_subdirectory(MOM5_SIS)
add_subdirectory(MOM5_CM2M)
add_subdirectory(MOM5_ESM2M)
add_subdirectory(MOM5_ICCM)
add_subdirectory(MOM5_EBM)
