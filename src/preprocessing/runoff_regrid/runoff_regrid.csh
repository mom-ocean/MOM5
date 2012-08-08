#!/bin/tcsh 
#=======================================================================
#      regrid_2d : interp 2-D data
#         Contact : Zhi Liang   email : z1l
#
#  This runscritp can be remap runoff data from spherical grid onto any grid
#  ( spherical or tripolar ) using conservative scheme. The runscript will
#  run executable $create_grid to generate atmos_grid based on the $src_data.
#  Then generate exchange grid grid_spec.nc by running $make_xgrids.
#  $make_xgrids need file atmos_grid.nc and file $dst_grid as command line
#  input. Finally remap the runoff data onto $dst_grid and write out output
#  file by executing $runoff_regrid.
#
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform      = "ncrc.intel"                                # A unique identifier for your platform

#  set ocean_only_dst_grid                  # comment this if the $dst_grid contains exchange grid information
                                            # uncomment this if the $dst_grid only contains ocean grid information.
			                    # $make_xgrids only accept ocean only grid as the -o option.
 
#
  set root         = $cwd:h:h:h                         # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 
  set make_xgrids  = $root/bin/make_xgrids_$platform
# needed data set and field name of the runoff data
# NOTE: Users must change the following according to their need. These particular values are here for testing purposes only.
#
  set dst_grid      = $FMS_ARCHIVE/mom4/mom4p1/memphis/mom4_om1p5/preprocessing/om1p5_grid_spec.nc
  set src_data      = $FMS_ARCHIVE/mom4/mom4p1/omsk/om3_core1/INPUT/RUNOFF.nc
  set src_fld_name  = 'RUNOFF'

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set name          = "runoff_regrid"                    # name of the data file will be generated
  set tooldir       = $cwd                               # directory of the tool
  set workdir       = $tooldir/workdir                   # where the tool is run and output is produced
  set create_grid   = $tooldir/exec/create_grid.exe      # executable created after compilation create_grid.F90
  set runoff_regrid = $tooldir/exec/runoff_regrid.exe    # executable created after compilation runoff_regrid.F90
  if ( $?DEBUG ) then
    set create_grid   = $tooldir/debug/debug_create_grid.exe 
    set runoff_regrid = $tooldir/debug/debug_runoff_regrid.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )                     # list of cpp #defines to be passed to the source files
 
# compile create_grid.F90 to create the executable.
  if( ! -d $create_grid:h ) mkdir $create_grid:h
  cd $create_grid:h
  $mkmf -t $mkmfTemplate -p $create_grid:t -c "$cppDefs"  $tooldir/create_grid.F90 /usr/local/include
  make $create_grid:t

# if make_xgrids does not exist, compare the source code
  if ( ! -f $make_xgrids ) then
     cc -O2 -o $make_xgrids $tooldir:h/generate_grids/make_xgrids/make_xgrids.c -L/usr/local/lib -lnetcdf -lm
  endif

# compile runoff_field.F90 to create the executable.
  if( ! -d $runoff_regrid:h ) mkdir $runoff_regrid:h
  cd $runoff_regrid:h
  $mkmf -t $mkmfTemplate -p $runoff_regrid:t -c "$cppDefs"  $tooldir/runoff_regrid.F90 /usr/local/include
  make $runoff_regrid:t

# compile runoff_regrid.F90 to create the executable.

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir

  cd $workdir

# get executable  
  cp $create_grid $create_grid:t
  cp $runoff_regrid $runoff_regrid:t
  cp $make_xgrids make_xgrids

# --- set up namelist

    cat >input.nml <<!
    &create_grid_nml
      src_data     = '$src_data',
      src_fld_name = '$src_fld_name'
      /
    &runoff_regrid_nml
      src_fld_name = '$src_fld_name',
      src_data     = '$src_data',
      dst_data     = '$name.nc'
       /
!

# first generate the atmosphere grid to be used to generate exchange grid

  if ( $?DEBUG ) then
     totalview $create_grid:t >fms.out
  else
     $create_grid:t >fms.out
  endif

  if ( $status != 0 ) then
    unset echo
    echo "ERROR: create_grid failed"
    exit 1
  endif

# then generate the exchange grid between dst_grid and the atmos_grid just generated
  if( $?ocean_only_dst_grid ) then
     make_xgrids -a atmos_grid.nc -l atmos_grid.nc -o $dst_grid >> fms.out
  else

# if $dst_grid contains exchange grid and atmos, land grid information, extract
# ocean grid from $dst_grid.
     ncea -x -v AREA_ATMxOCN,DI_ATMxOCN,DJ_ATMxOCN,I_ATM_ATMxOCN,J_ATM_ATMxOCN,I_OCN_ATMxOCN,J_OCN_ATMxOCN,AREA_ATMxLND,DI_ATMxLND,DJ_ATMxLND,I_ATM_ATMxLND,J_ATM_ATMxLND,I_LND_ATMxLND,J_LND_ATMxLND,AREA_LNDxOCN,DI_LNDxOCN,DJ_LNDxOCN,I_LND_LNDxOCN,J_LND_LNDxOCN,I_OCN_LNDxOCN,J_OCN_LNDxOCN,xba,yba,xta,yta,AREA_ATM,xbl,ybl,xtl,ytl,AREA_LND,AREA_LND_CELL,xto,yto,AREA_OCN $dst_grid $dst_grid:t
     make_xgrids -a atmos_grid.nc -l atmos_grid.nc -o $dst_grid:t >> fms.out
     rm -f $dst_grid:t
  endif

  if ( $status != 0 ) then
    unset echo
    echo "ERROR: make_xgrids failed"
    exit 1
  endif

# remove atm_grid.nc
  rm -f atmos_grid.nc

# Finally remap the runoff data onto destination grid
     $runoff_regrid:t >> fms.out

  if ( $status != 0 ) then
    unset echo
    echo "ERROR: runoff_regrid failed"
    exit 1
  endif

  cat fms.out

# remove the exchange grid
  rm -f grid_spec.nc
  
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out

  unset echo
