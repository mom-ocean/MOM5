#!/bin/tcsh 
#=======================================================================
#      regrid : interp 2-D/3-D logical rectangular grid to logical rectangular grid.
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
  set platform     = "ncrc.intel"      # A unique identifier for your platform
  set npes         = 1        # number of processors
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
  set src_data     = $FMS_ARCHIVE/mom4/mom4p1/riga/preprocess_tools/ocean.000101-000512.uv.nc
  set src_grid     = $FMS_ARCHIVE/mom4/mom4p1/riga/preprocess_tools/grid_spec_v5.nc
  set dst_grid     = $FMS_ARCHIVE/mom4/mom4p1/riga/preprocess_tools/grid_spec_v7.nc

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################

  set name         = "regrid"                               # name of the data file will be generated
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/regrid.exe               # executable created after compilation
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/regrid.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
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
      src_data = '$src_data',
      src_grid = '$src_grid',
      dst_grid = '$dst_grid',
      dst_data = '$name.nc',
      num_flds = 2
      fld_name = 'u','v'
      fld_pos  =  'C','C'
      vector_fld = .true., .true.
      use_source_vertical_grid = .false.
      apply_mask = .true.
      debug      = .false. /
!

#  run the executable
  $mpirunCommand $npes ./$executable:t >fms.out
  
  cat fms.out
       
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo
