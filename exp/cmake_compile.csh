#!/bin/csh

# Test compile script for cmake

set generator  = "Unix Makefiles"
set build_type = "relwithdebinfo"
set target     = "MOM5_SIS"
set help       = 0 

set argv = (`getopt -u -o h -l generator: -l build_type: -l target: -l help --  $*`)

while ("$argv[1]" != "--")
    switch ($argv[1])
        case --generator:
                set generator = $argv[2]; shift argv; breaksw
        case --build_type:
                set build_type = $argv[2]; shift argv; breaksw
        case --target:
                set target = $argv[2]; shift argv; breaksw
        case --help:
                set help = 1;  breaksw
        case -h:
                set help = 1;  breaksw
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
    echo
    exit 1
endif

set root          = $cwd:h                            # The directory you created when you checkout
set exec_dir      = $root/exec                        # source code directory
set build_dir     = $exec_dir/$build_type             # source code directory

cmake -E make_directory $build_dir
(cd $build_dir && cmake -G "$generator" -DCMAKE_BUILD_TYPE=$build_type $root/cmake && cmake --build . --target $target --config "$generator" -- -j10)
