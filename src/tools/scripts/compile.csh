#!/bin/csh 
set echo

source /usr/local/Modules/default/init/tcsh
module purge
module load intel_compilers/11.1.073 mpich2/1.2.1p1 netcdf/4.0.1
alias make make HDF5_HOME=/usr/local/hdf5-1.8.5-patch1_optimized NETCDF_HOME=/usr/local/netcdf-4.0.1_optimized SITE=pan -f fre-nctools.mk


set bindir   = $HOME/bin/tools_siena_201205
set tooldir  = $PWD:h

cd $tooldir

foreach tool (check_mask fregrid make_coupler_mosaic make_hgrid make_regional_mosaic make_solo_mosaic make_topog make_vgrid remap_land river_regrid runoff_regrid transfer_to_mosaic_grid)
  cd $tool
  rm -f *.o
  make
  if ( $status != 0 ) then
    unset echo
    echo "ERROR: make failed for $tool"
    exit 1
  endif
  cp  $tool $bindir
  set parallel_tool = ${tool}_parallel
  if( -e ${parallel_tool} ) then
     cp  ${parallel_tool} $bindir
  endif
  cd ..
end
