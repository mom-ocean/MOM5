#!/bin/csh
#contact: Zhi Liang (Zhi.Liang)

set echo
set release  = siena_201205
set bindir   = /home/z1l/bin/tools_$release
set basedir  = /archive/z1l/tools/test_$release
set npes     = 10
set inputdir              = /archive/z1l/tools/input
set om3_grid_dir          = $basedir/om3_grid
set dome_grid_dir         = $basedir/dome_grid
set bowl_grid_dir         = $basedir/bowl_grid
set channel_grid_dir      = $basedir/channel_grid
set torus_grid_dir        = $basedir/torus_grid
set gyre_grid_dir         = $basedir/gyre_grid
set coupled_grid_dir      = $basedir/coupled_grid
set C48_to_N45_dir        = $basedir/C48_to_N45
set cmip5_regrid_dir      = $basedir/cmip5_regrid
set regrid_extrap_dir     = $basedir/regrid_with_extrap
set runoffdir             = $basedir/runoff_regrid
set riverregriddir        = $basedir/river_regrid
set checkmask_baltic1_dir = $basedir/checkmask_baltic1
set checkmask_coupled_dir = $basedir/checkmask_coupled
set remaplanddir          = $basedir/remap_land
set regional_regrid_dir   = $basedir/regional_regrid

source /usr/local/Modules/default/init/tcsh
module purge
module load intel_compilers/11.1.073 mpich2/1.2.1p1 netcdf/4.0.1
module load nccmp

#**********************************************************************************
#
#    Test 1: om3 kind of ocean_grid
#
#
#**********************************************************************************
if( ! -d $om3_grid_dir ) mkdir -p $om3_grid_dir
cd $om3_grid_dir
#create hgrid
$bindir/make_hgrid --grid_type tripolar_grid --nxbnd 2 --nybnd 7 --xbnd -280,80 --ybnd -82,-30,-10,0,10,30,90 --dlon 1.0,1.0 --dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 --grid_name ocean_hgrid --center c_cell
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for om3_grid"
   exit 1
endif

#create vgrid
$bindir/make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for om3_grid"
   exit 1
endif

#create solo mosaic
$bindir/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc --periodx 360
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for om3_grid"
   exit 1
endif

#create topography data
$bindir/make_topog --mosaic ocean_mosaic.nc --topog_type realistic --topog_file /archive/fms/mom4/input_data/OCCAM_p5degree.nc --topog_field TOPO --scale_factor -1 --vgrid ocean_vgrid.nc --output topog.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for om3_grid"
   exit 1
endif

#test the reproducing ability of make_topog
mpirun -np 10 $bindir/make_topog_parallel --mosaic ocean_mosaic.nc --topog_type realistic --topog_file /archive/fms/mom4/input_data/OCCAM_p5degree.nc --topog_field TOPO --scale_factor -1 --vgrid ocean_vgrid.nc --output topog.10pe.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_topog_parallel for om3_grid"
   exit 1
endif

nccmp -md topog.nc topog.10pe.nc
if( $status != 0 ) then
    echo "ERROR: make_topog could not reproduce from processor count for om3_grid"
    exit 1
endif
  
#***********************************************************************************
#
#    Test 2: ocean_grid for dome experiment
#
#**********************************************************************************
if( ! -d $dome_grid_dir ) mkdir -p $dome_grid_dir
cd $dome_grid_dir
#create hgrid
$bindir/make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,25 --ybnd 63,70 --dlon 0.5,0.5 --dlat 0.5,0.5 --grid_name ocean_hgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for dome_grid"
   exit 1
endif

#create vgrid
$bindir/make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for dome_grid"
   exit 1
endif

#create solo mosaic
$bindir/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for dome_grid"
   exit 1
endif

#create topography data
$bindir/make_topog --mosaic ocean_mosaic.nc --topog_type dome --dome_slope 0.01 --dome_embayment_south 69 --dome_embayment_west 19.25 --dome_embayment_east 21.25 --dome_bottom 5500
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for dome_grid"
   exit 1
endif

#***********************************************************************************
#
#    Test 3: ocean_grid for bowl experiment
#
#**********************************************************************************
if( ! -d $bowl_grid_dir ) mkdir -p $bowl_grid_dir
cd $bowl_grid_dir
#create hgrid
$bindir/make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,20 --ybnd 60,72 --dlon 0.5,0.5 --dlat 0.5,0.5 --grid_name ocean_hgrid
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for bowl_grid"
   exit 1
endif

#create vgrid
$bindir/make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for bowl_grid"
   exit 1
endif

#create solo mosaic
$bindir/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for bowl_grid"
   exit 1
endif

#create topography data
$bindir/make_topog --mosaic ocean_mosaic.nc --topog_type bowl --bottom_depth 5500
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for bowl_grid"
   exit 1
endif

#***********************************************************************************
#
#    Test 4: ocean_grid for box channel experiment  
#
#**********************************************************************************
if( ! -d $channel_grid_dir ) mkdir -p $channel_grid_dir
cd $channel_grid_dir
#create hgrid
$bindir/make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,10 --ybnd -40,-10 --dlon 0.5,0.5 --dlat 0.5,0.5 --grid_name ocean_hgrid
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for channel_grid"
   exit 1
endif

#create vgrid
$bindir/make_vgrid --nbnds 2 --bnds 0,1200 --dbnds 40,40 --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for channel_grid"
   exit 1
endif

#create solo mosaic
$bindir/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for channel_grid"
   exit 1
endif

#create topography data
$bindir/make_topog --mosaic ocean_mosaic.nc --topog_type box_channel --jwest_south 5 --jwest_north 15 --jeast_south 30 --jeast_north 45 --bottom_depth 1200
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for channel_grid"
   exit 1
endif

#***********************************************************************************
#
#    Test 5: ocean_grid for torus experiment  
#
#**********************************************************************************
if( ! -d $torus_grid_dir ) mkdir -p $torus_grid_dir
cd $torus_grid_dir
#create hgrid
$bindir/make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,20 --ybnd 10,30 --dlon 1,1 --dlat 1,1 --grid_name ocean_hgrid
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for torus_grid"
   exit 1
endif

#create vgrid
$bindir/make_vgrid --nbnds 2 --bnds 0.,200 --dbnds 10,10 --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for torus_grid"
   exit 1
endif

#create solo mosaic
$bindir/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc --periodx 20 --periody 20
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for torus_grid"
   exit 1
endif

#create topography data
$bindir/make_topog --mosaic ocean_mosaic.nc --topog_type rectangular_basin --bottom_depth 200
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for torus_grid"
   exit 1
endif

#***********************************************************************************
#
#    Test 6: ocean_grid for gyre experiment   
#
#**********************************************************************************
if( ! -d $gyre_grid_dir ) mkdir -p $gyre_grid_dir
cd $gyre_grid_dir
#create hgrid
$bindir/make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,10 --ybnd 15,35 --dlon 2,2 --dlat 2,2 --grid_name ocean_hgrid
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for gyre_grid"
   exit 1
endif

#create vgrid
$bindir/make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for gyre_grid"
   exit 1
endif

#create solo mosaic
$bindir/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for gyre_grid"
   exit 1
endif

#create topography data
$bindir/make_topog --mosaic ocean_mosaic.nc --topog_type rectangular_basin --bottom_depth 5500
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for gyre_grid"
   exit 1
endif

#***********************************************************************************
#
#    Test 7: grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid.
#
#**********************************************************************************
if( ! -d $coupled_grid_dir ) mkdir -p $coupled_grid_dir
cd $coupled_grid_dir
#create ocean_hgrid 
$bindir/make_hgrid --grid_type tripolar_grid --nxbnd 2 --nybnd 7 --xbnd -280,80 --ybnd -82,-30,-10,0,10,30,90 --dlon 1.0,1.0 --dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 --grid_name ocean_hgrid --center c_cell
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for coupled_grid ocean"
   exit 1
endif

#create ocean_vgrid
$bindir/make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for coupled_grid ocean"
   exit 1
endif

#create ocean solo mosaic
$bindir/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc --periodx 360
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for coupled_grid ocean"
   exit 1
endif

#create ocean topography data
$bindir/make_topog --mosaic ocean_mosaic.nc --topog_type realistic --topog_file /archive/fms/mom4/input_data/OCCAM_p5degree.nc --topog_field TOPO --scale_factor -1 --vgrid ocean_vgrid.nc --output topog.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for coupled_grid ocean"
   exit 1
endif

#Create C48 grid for atmos/land
$bindir/make_hgrid --grid_type gnomonic_ed --nlon 96 --grid_name C48_grid
if( $status != 0 ) then
    echo "ERROR:  failed at make_hgrid for C48_grid"
    exit 1
endif
#create C48 solo mosaic for atmos/land
$bindir/make_solo_mosaic --num_tiles 6 --dir ./ --mosaic C48_mosaic --tile_file C48_grid.tile1.nc,C48_grid.tile2.nc,C48_grid.tile3.nc,C48_grid.tile4.nc,C48_grid.tile5.nc,C48_grid.tile6.nc
if( $status != 0 ) then
    echo "ERROR: failed at make_solo_mosaic for C48_grid"
    exit 1
endif

#make the coupler_mosaic
$bindir/make_coupler_mosaic --atmos_mosaic C48_mosaic.nc --ocean_mosaic ocean_mosaic.nc --ocean_topog  topog.nc --mosaic_name grid_spec
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_coupler_mosaic for coupled_grid"
    exit 1
endif

#check reproducing ability between processor count for make_coupler_mosaic
if( ! -d multiple_pe ) mkdir multiple_pe
cd multiple_pe
mpirun -np 10 $bindir/make_coupler_mosaic_parallel --atmos_mosaic ../C48_mosaic.nc --ocean_mosaic ../ocean_mosaic.nc --ocean_topog  ../topog.nc --mosaic_name grid_spec
  
if( $status != 0 ) then
    echo "ERROR: failed at make_coupler_mosaic on 10 pe for coupled_grid"
    exit 1
endif
foreach file ( `ls -1 *.nc` )
   nccmp -md $file ../$file
   if( $status != 0 ) then
      echo "ERROR: make_coupler_mosaic could not reproduce for $file"
      exit 1
   endif
endif

#***********************************************************************************
#
#    Test 8: remap data from C48 to regular lat-lon grid
#
#**********************************************************************************
if( ! -d $C48_to_N45_dir ) mkdir -p $C48_to_N45_dir
cd $C48_to_N45_dir
$bindir/fregrid --input_mosaic $inputdir/C48/C48_mosaic.nc --input_dir $inputdir --input_file 19800101.atmos_daily --scalar_field zsurf,temp,t_surf --nlon 144 --nlat 90 --interp_method conserve_order2 --output_file 19800101.atmos_daily.nc --check_conserve --remap_file C48_to_N45_remap.nc

if( $status != 0 ) then
    unset echo
    echo "ERROR: run failed for fregrid from C48 to N45 using conservative_order2"
    exit 1
endif

#running on multiple processors
mpirun -np $npes $bindir/fregrid_parallel --input_mosaic $inputdir/C48/C48_mosaic.nc --input_dir $inputdir --input_file 19800101.atmos_daily --scalar_field zsurf,temp,t_surf --nlon 144 --nlat 90 --interp_method conserve_order2 --output_file 19800101.atmos_daily.${npes}pe.nc --check_conserve --remap_file C48_to_N45_remap.nc

if( $status != 0 ) then
    unset echo
    echo "ERROR: failed for fregrid from C48 to N45 using conservative_order2 on $npes pe"
    exit 1
endif

nccmp -md 19800101.atmos_daily.nc 19800101.atmos_daily.${npes}pe.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: fregrid could not reproduce between 1 and $npes processors"
    exit 1
endif

#test reading from remap file
mpirun -np $npes $bindir/fregrid_parallel --input_mosaic $inputdir/C48/C48_mosaic.nc --input_dir $inputdir --input_file 19800101.atmos_daily --scalar_field zsurf,temp,t_surf --nlon 144 --nlat 90 --interp_method conserve_order2 --output_file 19800101.atmos_daily.${npes}pe.read.nc --check_conserve --remap_file C48_to_N45_remap.nc

if( $status != 0 ) then
    unset echo
    echo "ERROR: failed for fregrid from C48 to N45 using conservative_order2 on $npes pe and read remap file"
    exit 1
endif

nccmp -md 19800101.atmos_daily.${npes}pe.read.nc 19800101.atmos_daily.${npes}pe.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: fregrid could not reproduce when reading from remap file"
    exit 1
endif

#***********************************************************************************
#
#    Test 9: remap cmip5 ocean data onto regular lat-lon grid (use GFDL-CM3 data as example)
#
#**********************************************************************************
if( ! -d $cmip5_regrid_dir ) mkdir -p $cmip5_regrid_dir
cd $cmip5_regrid_dir
# create GFDL-CM3_grid.nc
$bindir/make_hgrid --grid_type from_file --my_grid_file $inputdir/GFDL-CM3_sic_OImon.nc --grid_name GFDL-CM3_grid
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_hgrid for cmip5 ocean data"
    exit 1
endif

#. create GFDL-CM3_mosaic.nc
$bindir/make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name GFDL-CM3_mosaic --tile_file GFDL-CM3_grid.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_solo_mosaic for cmip5 ocean data"
    exit 1
endif

#. create remap file 
$bindir/fregrid --input_mosaic GFDL-CM3_mosaic.nc --nlon 180 --nlat 90 --remap_file remap_file.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at creating remap_file for cmip5 ocean data"
    exit 1
endif

#. remap data
$bindir/fregrid --input_mosaic GFDL-CM3_mosaic.nc --nlon 180 --nlat 90 --remap_file remap_file.nc --input_dir $inputdir --input_file GFDL-CM3_sic_OImon.nc --output_file sic_OImon.out.nc --scalar_field sic
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at regrid for cmip5 ocean data"
    exit 1
endif

#***********************************************************************************
#
#    Test 10: remap data onto cm2m ocean grid with extrapolation and vertical interpolation.
#
#**********************************************************************************
#create levitus_grid and levitus_mosaic
if( ! -d $regrid_extrap_dir ) mkdir -p $regrid_extrap_dir
cd $regrid_extrap_dir
$bindir/make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,360 --ybnd -90,90 --nlon 720 --nlat 360 --grid_name levitus_grid
if( $status != 0 ) then
    echo "ERROR: failed at make_hgrid for levitus grid"
    exit 1 
endif 
$bindir/make_solo_mosaic --num_tiles 1 --dir . --mosaic_name levitus_mosaic --tile_file levitus_grid.nc --periodx 360
if( $status != 0 ) then
    echo "ERROR: failed at make_solo_mosaic for levitus grid"
    exit 1 
endif 

$bindir/fregrid --input_mosaic levitus_mosaic.nc --input_dir $inputdir --input_file WOA09_ann_theta.nc --scalar_field POTENTIAL_TEMP --output_file WOA09_ann_theta_cm2g_extrap.nc --output_mosaic $inputdir/cm2m_grid/ocean_mosaic.nc --extrapolate --dst_vgrid $inputdir/cm2m_grid/ocean_vgrid.nc
if( $status != 0 ) then
    echo "ERROR: failed at regrid data onto cm2m ocean grid with extrapolation and vertical interpolation."
    exit 1 
endif 


#***********************************************************************************
#
#    Test 11: remap runoff data from regular lat-lon grid onto cm2m grid
#
#**********************************************************************************
if( ! -d $runoffdir ) mkdir -p $runoffdir
cd $runoffdir
$bindir/runoff_regrid --input_file $inputdir/runoff.daitren.iaf.nc --input_fld_name runoff --output_mosaic $inputdir/cm2m_grid/ocean_mosaic.nc --output_topog $inputdir/cm2m_grid/topog.nc --output_file runoff.cm2m.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at runoff_regrid from lat-lon grid onto cm2m ocean grid"
    exit 1
endif


#***********************************************************************************
#
#    Test 12: river_regrid: remap data from lat-lon onto C48 grid.
#
#**********************************************************************************
if( ! -d $riverregriddir ) mkdir -p $riverregriddir
cd $riverregriddir

$bindir/river_regrid --mosaic $inputdir/C48/mosaic.nc --river_src $inputdir/z1l_river_output_M45_tripolar_aug24.nc --output river_data_C48 --land_thresh 0.000001
if( $status != 0 ) then 
    unset echo 
    echo "ERROR: failed for river_regrid to regrid river data onto C48 mosaic grid."
    exit 1 
endif


#***********************************************************************************
#
#    Test 13: check_mask: create mask_table to mask out some all-land domain
#                         to save processor usage
# 
#             for a sea-ice model, baltic1 experiment
#**********************************************************************************
if( ! -d $checkmask_baltic1_dir ) mkdir -p $checkmask_baltic1_dir
cd $checkmask_baltic1_dir
$bindir/check_mask --grid_file $inputdir/baltic1_grid_spec.nc --min_pe 60 --max_pe 80 

if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at check_mask for baltic1 experiment"
    exit 1
endif

#***********************************************************************************
#
#    Test 14: check_mask: create mask_table to mask out some all-land domain
#                         to save processor usage
# 
#             for a coupled model with cm2m grid.
#**********************************************************************************
if( ! -d $checkmask_coupled_dir ) mkdir -p $checkmask_coupled_dir
cd $checkmask_coupled_dir
$bindir/check_mask --grid_file $inputdir/cm2m_grid/ocean_mosaic.nc --ocean_topog $inputdir/cm2m_grid/topog.nc --min_pe 100 --max_pe 200
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at check_mask for a coupled model grid"
    exit 1
endif
#***********************************************************************************
#
#    Test 15: remap_land: remap land restart file.
#
#**********************************************************************************
if( ! -d $remaplanddir ) mkdir -p $remaplanddir
cd $remaplanddir
$bindir/remap_land --src_mosaic $inputdir/C48/C48_mosaic.nc --dst_mosaic $inputdir/C180/C180_mosaic.nc --src_restart $inputdir/land_src_restart/soil.res --dst_cold_restart $inputdir/land_dst_cold_restart/soil.res --dst_restart soil.res
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at remap_land from C48 to C180"
    exit 1
endif

if( ! -d multiple_pe ) mkdir multiple_pe
cd multiple_pe
mpirun -n $npes $bindir/remap_land_parallel --src_mosaic $inputdir/C48/C48_mosaic.nc --dst_mosaic $inputdir/C180/C180_mosaic.nc --src_restart $inputdir/land_src_restart/soil.res --dst_cold_restart $inputdir/land_dst_cold_restart/soil.res --dst_restart soil.res
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at remap_land from C48 to C180 on parallel job"
    exit 1
endif
foreach file (`ls -1 soil*.nc`)
   nccmp -md $file ../$file
   if( $status != 0 ) then
      echo "ERROR: remap_land could not reproduce  between processor count for $file"
      exit 1
   endif


#***********************************************************************************
#
#    Test 16: make_regional_mosaic: create mosaic file for regional output
#
#**********************************************************************************
if( ! -d $regional_regrid_dir ) mkdir -p $regional_regrid_dir
cd $regional_regrid_dir

$bindir/make_regional_mosaic --global_mosaic $inputdir/C2560/C2560_mosaic.nc --regional_file $inputdir/rregionatmos_4xdaily.tile1.nc
$bindir/fregrid --input_mosaic regional_mosaic.nc --input_dir $inputdir --input_file rregionatmos_4xdaily.tile1.nc --lonBegin 0  --lonEnd 5 --latBegin -20 --latEnd 20 --nlon 60 --nlat 480 --scalar_field temp
