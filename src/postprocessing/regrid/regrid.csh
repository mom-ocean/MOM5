#!/bin/tcsh 
#=======================================================================
#      regrid_2d : interp 2-D data
#         Contact : Zhi Liang   email : z1l
#
#  This runscritp can be used to regrid data $src_data on any grid $src_grid to logically
#  rectangular grid described by grid descriptor file $dst_grid. 
#
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform     = "ncrc.intel"                               # A unique identifier for your platform
  set npes         = 1                                  # number of processors
#
  set root         = $cwd:h:h:h                         # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 
# input data file and destination grid
  set src_data       = $FMS_ARCHIVE/z1l/CM2Q-d2_1PctTo4x_j1/ocean.000101-000512.eta_t.nc
  set grid_spec_file = $FMS_ARCHIVE/z1l/postprocessing/regrid/input/om3ocn_1degatm.nc

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set name         = "regrid_eta_t"                         # name of the data file will be generated
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/regrid.exe               # executable created after compilation
  set mkmf         = $root/bin/mkmf                     # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )        # list of cpp #defines to be passed to the source files

# list the source code
  set CORE      = "$tooldir/{regrid.F90}"
  set UTILITIES = "$sharedir/{axis_utils,constants,fms,mpp,horiz_interp,platform,memutils,mosaic}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir " 
  set srclist   = ( $CORE $UTILITIES )

# compile the code and create executable
  if( ! -d $executable:h ) mkdir $executable:h
  cd $executable:h
  $mkmf -t $mkmfTemplate -p $executable:t -c "$cppDefs" $srclist /usr/local/include 

  make $executable:t

  if($status != 0) then
     unset echo
     echo " Error in compilation "
     exit 1
  endif

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir

  cd $workdir

# get executable  
  cp $executable $executable:t

# --- set up namelist

    cat >input.nml <<!
    &regrid_nml
       src_data       = '$src_data',
       grid_spec_file = '$grid_spec_file',
       dst_data       = '$name.nc',
       num_flds       = 1
       fld_name       = 'eta_t'
       fld_pos        =  'T'
       vector_fld     = .false.
       debug          = .false. /
!

#  run the executable
  $mpirunCommand $npes ./$executable:t >fms.out 
  cat fms.out
       
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out

  unset echo
