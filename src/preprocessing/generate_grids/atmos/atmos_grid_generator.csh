#!/bin/tcsh 
#=======================================================================
#      preprocessing : atmos_grid_generator
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscritp can be used to generate bgrid or spectral horizontal grid
#  for atmosphere model or land model. It creates Makefile, compile if necessary
#  and saves output. The generated data sets are atmos_grid.nc (default value of
#  namelist varible output_file of atmos_grid_generator_nml) and are saved in
#  the directory $output_dir/data. The standard output is stored at
#  $output_dir/ascii/.
#=======================================================================
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform     = "ncrc.intel"                       # A unique identifier for your platform
#
  set root         = $cwd:h:h:h:h                       # The directory that contains src/ and bin/
  set npes         = 1
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set name         = "atmos_grid"                           # name of the grid file will be generated
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include                      # fms include directory  
  set atmtool      = $root/src/atmos_spectral/tools         # directory contains atmos_spectral tool code.
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/atmos_grid_generator.exe # executable created after compilation
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/atmos_grid_generator.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI -DOVERLOAD_C8" )        # list of cpp #defines to be passed to the source files

# list the source code
  set CORE         = " $tooldir/{atmos_grid.f90,atmos_grid_generator.f90} "
  set CORE         = " $CORE $atmtool"
  set UTILITIES    = "$root/src/shared/{axis_utils,constants,fms,mpp,fft,platform,memutils,mosaic}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir "
  set srclist      =  ( $CORE $UTILITIES )

   
# compile the model code and create executable
  if( ! -d $executable:h ) mkdir $executable:h
  cd $executable:h
  $mkmf -t $mkmfTemplate -p $executable:t -c "$cppDefs" $srclist 

  make $executable:t

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir

  cd $workdir

# get executable  
  cp $executable $executable:t

# --- set up namelist
    cat >input.nml <<!
     &atmos_grid_generator_nml
       output_file = '$name.nc'  /    
    &atmos_grid_nml
       grid_type = 'bgrid'
       num_lon=144, num_lat=90   /
!

#  run the executable
  $mpirunCommand $npes ./$executable:t >fms.out
  cat fms.out
       
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo

