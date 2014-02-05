#!/bin/csh 
#=======================================================================
#      preprocessing : make_xgrids.csh
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscritp can be used to generate exchange grids between ocean, atmosphere and land grid.
#  The user need to create atmos_grid.nc, land_grid.nc and ocean_grid.nc before using make_xgrids.
#  Normally land_grid.nc is the same as atmos_grid.nc. 
#=======================================================================
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################
  set echo
  set name         = "grid_spec"                                  # name of the grid file will be generated
  set platform     = "ncrc.intel"                                    # A unique identifier for your platform
  set npes         = 1                                      # number of processors
#
  set root         = $cwd:h:h:h:h                       # The directory that contains src/ and bin/
  set tooldir      = $cwd                                         # directory of the tool 
  set workdir      = $tooldir/workdir                             # where the tool is run and output is produced
  set executable   = $tooldir/exec/make_xgrids                    # executable created after compilation
  set xgrids_code  = $tooldir/make_xgrids.c                       # source code of make_xgrids
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

# model grid data
  set atmos_grid   = $tooldir:h/atmos/workdir/atmos_grid.nc       # atmos grid
  set land_grid    = $tooldir:h/atmos/workdir/atmos_grid.nc       # land grid
  set ocean_grid   = $tooldir:h/ocean/workdir/ocean_grid.nc       # ocean grid

#--create the executable -----------------------------------------------------------
  if( ! -d $executable:h ) mkdir $executable:h
  cd $executable:h
  cc -O -o $executable:t $xgrids_code -I$netcdf3_inc_dir -L$netcdf3_lib_dir -lnetcdf -lm

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir

  cd $workdir

# get executable  
  cp $executable $executable:t

#  run the executable
  ./$executable:t -o $ocean_grid -a $atmos_grid -l $land_grid >>fms.out   
  cat fms.out
  
# renmae output file
  mv fms.out $name.fms.out
  mv grid_spec.nc $name.nc

  unset echo 
