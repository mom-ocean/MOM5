#!/bin/tcsh 
#=======================================================================
#      grid_transer 
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscritp can be used to convert the old name convention grid specification
#  netcdf file to new name convention grid netcdf file. The user need to provide the
#  input grid file $old_grid and the program will generate the output grid file
#  $new_grid. The output file $new_grid contains both the new name convention grid
#  infofmation and all the fields in $old_grid.
#
#=======================================================================
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform     = "ncrc.intel"                                # A unique identifier for your platform
  set npes         = 1                                      # number of processors

#
  set root         = $cwd:h:h:h:h                       # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options
# grid file to be converted
  set old_grid     = "$FMS_ARCHIVE/mom4/mom4p1/nalanda/mom4_om3_core/preprocessing/grid_spec_v7.nc"
                    # grid file you want to transfer to new format.

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 
  set name         = "grid_transfer"                    # name of the grid file will be generated
  set tooldir      = $cwd                               # directory of the tool
  set sharedir     = $root/src/shared                   # directory of the shared code.
  set includedir   = $sharedir/include                      # fms include directory
  set workdir      = $tooldir/workdir                   # where the tool is run and output is produced
  set mkmf         = $root/bin/mkmf                     # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )   # list of cpp #defines to be passed to the source files
  set executable   = $tooldir/exec/grid_transfer.exe    # executable created after compilation


# list the source code
  
  set CORE         = "$tooldir/{grid_transfer.F90}"
  set UTILITIES    = "$sharedir/{axis_utils,constants,fms,mpp,platform,memutils}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir " 
  set srclist      =  ( $CORE $UTILITIES )                       # list of the source code


# compile the code and create executable
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
    &grid_transfer_nml
      old_grid = '$old_grid'
      new_grid = '$name.nc'
       /
!

#  run the executable
  $mpirunCommand $npes ./$executable:t >fms.out
  
  cat fms.out
       
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo


