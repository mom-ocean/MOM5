#!/bin/tcsh 
#=======================================================================
#      river_regrid : remap river network data.
#         Contact : Zhi Liang   email : Zhi.Liang
#  This program can remap river network data from spherical grid onto another
#  spherical grid. 
#
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform     = "ncrc.intel"                                    # A unique identifier for your platform
  set npes         = 1
#
  set root         = $cwd:h:h:h                         # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 
# needed data set and field name of the runoff data
# NOTE: Users must change the following according to their need. These particular values are here for testing purposes only.
#
 set river_input_file = $FMS_ARCHIVE/home/z1l/fms/src/preprocessing/river_regrid/input/orig_runoff.rivers_init.nc
 set land_grid_file   = $FMS_ARCHIVE/mom4/mom4p1/riga/preprocess_tools/grid_spec.partial_areas.1deg.nc

#############################################################################
# Users need not change anything below this line except the namelists values.
##############################################################################
  set name         = "river_regrid"                         # name of the data file will be generated
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include                          # fms include directory
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/river_regrid.exe            # executable created after compilation
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )                     # list of cpp #defines to be passed to the source files

# list the source code
  set CORE      = "$tooldir/{river_regrid.f90}"
  set UTILITIES = "$sharedir/{axis_utils,constants,fms,mpp,platform,memutils,horiz_interp,mosaic}"
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
     &river_regrid_nml
       river_input_file = '$river_input_file'
       land_grid_file   = '$land_grid_file' /
!


#  run the executable
     $mpirunCommand $npes ./$executable:t >fms.out

     cat fms.out
  
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo
