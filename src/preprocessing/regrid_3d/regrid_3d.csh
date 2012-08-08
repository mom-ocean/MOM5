#!/bin/tcsh 
#=======================================================================
#      regrid_3d : interp 3-D data
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscritp can be used to regrid 3-d lat-lon gridded data to
#  logically rectangular grid described by grid descriptor file $dest_file.
#  Applies only to scalar fields.
#
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
  set root         = $cwd:h:h:h                         # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template
#
# input data file and destination grid
#
  set src_file        = $FMS_ARCHIVE/mom4/input_data/levitus_ewg.nc                     # source data file
  set dest_grid       = $cwd:h/generate_grids/ocean/workdir/ocean_grid.nc # destination grid
 
#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set name         = "regrid_3d"                            # name of the data file will be generated
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include                      # fms include directory
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/regrid_3d.exe            # executable created after compilation
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/regrid_3d.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )        # list of cpp #defines to be passed to the source files


# list the source code
  set CORE      = "$tooldir/{regrid_3d.f90}"
  set UTILITIES = "$sharedir/{axis_utils,constants,fms,mpp,horiz_interp,platform,memutils,mosaic}"
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
    &regrid_3d_nml
       src_file = '$src_file'
       numfields = 2
       src_field_name = 'temp', 'salt'
       dest_grid = '$dest_grid'
       dest_file = '$name.nc'
       ntimes_saved = 1,
       timelevel_saved = 1,
       num_nbrs = 4 /
    &fms_nml
       domains_stack_size = 395850 /
!

#  run the executable
        $mpirunCommand $npes ./$executable:t >fms.out

  cat fms.out
       
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo

