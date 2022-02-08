#!/bin/csh

# Test compile script for cmake

set generator  = "Unix Makefiles"
set build_type = "relwithdebinfo"
set target     = "MOM5_SIS"
set use_netcdf4 = 1
set help       = 0 
set FC=mpifort
set CC=mpicc

set argv = (`getopt -u -o h -l generator: -l build_type: -l target: -l fcomp: -l ccomp: -l help --  $*`)

while ("$argv[1]" != "--")
    switch ($argv[1])
        case --generator:
                set generator = $argv[2]; shift argv; breaksw
        case --build_type:
                set build_type = $argv[2]; shift argv; breaksw
        case --target:
                set target = $argv[2]; shift argv; breaksw
        case --netcdf3:
                set use_netcdf4 = 0; breaksw
        case --help:
                set help = 1;  breaksw
        case -h:
                set help = 1;  breaksw
        case --fcomp:
                set FC = $argv[2]; shift argv; breaksw
        case --ccomp:
                set CC = $argv[2]; shift argv; breaksw
    endsw
    shift argv
end
shift argv
if ( $help ) then
    echo "The optional arguments are:"
    echo "--generator  followed by the CMake generator, one of the following (default is 'Unix Makefiles'):"
    echo "             'Unix Makefiles'    : Unix makefiles"
    echo "             Ninja               : fast and lightweight build system (cmake version 3.8.0 required)"
    echo "--build_type followed by the CMake build type, one of the following (default is relwithdebinfo):"
    echo "             relwithdebinfo      : release build flags with some debugging information"
    echo "             release             : release build"
    echo "             debug               : debug build, no optimisation"
    echo "--target     followed by the type of the model, one of the following (default is MOM5_SIS):"
    echo "             MOM5_solo  : solo ocean model"
    echo "             MOM5_SIS   : ocean-seaice model"
    echo "             MOM5_CM2M  : ocean-seaice-land-atmosphere coupled climate model"
    echo "             MOM5_ESM2M : ocean-seaice-land-atmosphere coupled climate model with biogeochemistry, EarthSystemModel"
    echo "             MOM5_ICCM  : ocean-seaice-land-atmosphere coupled model"
    echo "             MOM5_EBM   : ocean-seaice-land-atmosphere coupled model with energy balance atmosphere"
    echo "--netcdf3    force netCDF3, default is to use netCDF4"
    echo "--fcomp      Path to FORTRAN compiler"
    echo "--ccomp      Path to C compiler"
    echo ""
    echo "The compiler options are provided in the case where a site-specific compiler wrappers should be used to"
    echo "configure builds correctly, currently defaults to mpifort and mpicc"
    echo
    exit 1
endif

set root          = $cwd:h                            # The directory you created when you checkout
set exec_dir      = $root/exec                        # source code directory
set build_dir     = $exec_dir/$build_type             # source code directory

# Update version hash
source ./update_version.csh

echo cmake -E make_directory $build_dir
cmake -E make_directory $build_dir
echo cmake -G "$generator" -DCMAKE_BUILD_TYPE=$build_type -DNETCDF4=use_netcdf4 -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_C_COMPILER=$CC -S$root/cmake -B$build_dir
cmake -G "$generator" -DCMAKE_BUILD_TYPE=$build_type -DNETCDF4=use_netcdf4 -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_C_COMPILER=$CC -S$root/cmake -B$build_dir
echo cmake --build $build_dir --target $target --config "$generator" -j
cmake --build $build_dir --target $target --config "$generator" -j
