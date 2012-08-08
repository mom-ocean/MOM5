#!/bin/tcsh 
#=======================================================================
#      idealized_ic : generate idealized initial condition
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscript can be used to generate tracer initial conditions for mom4.
#  It creates Makefile, compiles if necessary and  saves output.The generated
#  data sets are saved in the directory $output_dir/data. A grid_spec file
#  is needed for this program. 
#  ******Important notice: to obtain desired initial condition, you need to
#  specify initial_nml which is listed in this runscript.
#=======================================================================
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform     = "ncrc.intel"                           # A unique identifier for your platform
  set npes         = 1                                      # number of processors
#
# input data file and destination grid
#
  set grid_spec_file  = $cwd:h:h/generate_grids/ocean/workdir/ocean_grid.nc  # destination grid
#
  set root         = $cwd:h:h:h:h                         # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set name         = "idealized_ic"                         # name of the tool
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include                      # fms include directory
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/idealized_ic.exe            # executable created after compilation
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/idealized_ic.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )        # list of cpp #defines to be passed to the source files

# list the source code
  set CORE      = $tooldir/{idealized_ic.f90,idealized_ic_driver.f90}
  set UTILITIES = "$sharedir/{fms,mpp,platform,constants,memutils}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir " 
  set srclist   = ( $CORE $UTILITIES )

# compile the code and create executable
  if( ! -d $executable:h ) mkdir $executable:h
  cd $executable:h
  $mkmf -t $mkmfTemplate -p $executable:t -c "$cppDefs" $srclist /usr/local/include

  make $executable:t

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir

  cd $workdir

# get executable  
  cp $executable $executable:t

# --- set up namelist

    cat >input.nml <<!
    &idealized_ic_nml
       temp_type = 'zonal_levitus_temp_ic',
       salt_type = 'salinity_profile_ic'
       constant_temp_value = 15.1 ,
       constant_salt_value = 32. 
       passive_tracer_file = 'ocean_age.res.nc'
       grid_file           = '$grid_spec_file' /
    &fms_nml
       domains_stack_size = 395850 /       
!

#  run the executable
  $mpirunCommand $npes $executable:t >fms.out
  
  cat fms.out
       
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo

