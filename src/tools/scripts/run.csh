#!/bin/csh
#contact: Zhi Liang (Zhi.Liang)
#PBS -N test_tools
#PBS -l nodes=1:ppn=6
#PBS -l walltime=3:00:00
#PBS -o $HOME/tools/tikal/tools/scripts
#PBS -j oe
#PBS -m abe
#PBS -r y
#PBS -q batch

#set echo
set release  = tikal
set platform = `perl -T -e "use Net::Domain(hostdomain); print hostdomain"`
if( $platform == "princeton.rdhpcs.noaa.gov" ) then
   set indir    = /archive/fms/tools/input
   set outdir   = /work/$USER/tools/test_$release
   set basedir  = $TMPDIR/tools
   module load mpich2/1.2.1p1
   alias aprun mpirun
else if( $platform == "ncrc.gov" ) then
   set indir    = /lustre/f1/pdata/fms/tools/input
   set outdir   = /lustre/f1/unswept/$USER/tools/test_$release
   set basedir  = /lustre/f1/$USER/work/tools/test_$release
   alias aprun aprun
endif

set npes = 4
set npes2 = 6

set dome_grid         = dome_grid
set bowl_grid         = bowl_grid
set channel_grid      = channel_grid
set torus_grid        = torus_grid
set gyre_grid         = gyre_grid
set coupled_grid      = coupled_grid
set coupled_nest_grid = coupled_nest_grid
set C48_to_N45        = C48_to_N45
set cmip5_regrid      = cmip5_regrid
set regrid_extrap     = regrid_extrap
set runoff_regrid     = runoff_regrid
set river_regrid      = river_regrid
set checkmask_baltic1 = checkmask_baltic1
set checkmask_MOM6    = checkmask_MOM6
set remap_land        = remap_land
set regional_regrid   = regional_regrid
set make_quick_mosaic = make_quick_mosaic
set mppncscatter      = mppncscatter
set transfer_mosaic   = transfer_to_mosaic_grid
set CM21_to_CM25      = CM2.1_to_CM2.5

module load fre-nctools/$release
module load nccmp

    
#***********************************************************************************
#
#    Test 1: ocean_grid for dome experiment
#
#**********************************************************************************
set workdir = $basedir/$dome_grid
set output  = $outdir/$dome_grid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#create hgrid
make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,25 --ybnd 63,70 --dlon 0.5,0.5 --dlat 0.5,0.5 --grid_name ocean_hgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for dome_grid"
   exit 1
endif

#create vgrid
make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for dome_grid"
   exit 1
endif

#create solo mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for dome_grid"
   exit 1
endif

#create topography data
make_topog --mosaic ocean_mosaic.nc --topog_type dome --dome_slope 0.01 --dome_embayment_south 69 --dome_embayment_west 19.25 --dome_embayment_east 21.25 --dome_bottom 5500
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for dome_grid"
   exit 1
endif

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 1: ocean_grid for dome experiment                       *"
echo "*                                                                         *"
echo "***************************************************************************"

#copy data to the output directory
mv *.nc $output 

#***********************************************************************************
#
#    Test 2: ocean_grid for bowl experiment
#
#**********************************************************************************
set workdir = $basedir/$bowl_grid
set output  = $outdir/$bowl_grid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#create hgrid
make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,20 --ybnd 60,72 --dlon 0.5,0.5 --dlat 0.5,0.5 --grid_name ocean_hgrid
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for bowl_grid"
   exit 1
endif

#create vgrid
make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for bowl_grid"
   exit 1
endif

#create solo mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for bowl_grid"
   exit 1
endif

#create topography data
make_topog --mosaic ocean_mosaic.nc --topog_type bowl --bottom_depth 5500
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for bowl_grid"
   exit 1
endif

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 2: ocean_grid for bowl experiment                       *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 3: ocean_grid for box channel experiment  
#
#**********************************************************************************
set workdir = $basedir/$channel_grid
set output  = $outdir/$channel_grid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#create hgrid
make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,10 --ybnd -40,-10 --dlon 0.5,0.5 --dlat 0.5,0.5 --grid_name ocean_hgrid
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for channel_grid"
   exit 1
endif

#create vgrid
make_vgrid --nbnds 2 --bnds 0,1200 --dbnds 40,40 --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for channel_grid"
   exit 1
endif

#create solo mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for channel_grid"
   exit 1
endif

#create topography data
make_topog --mosaic ocean_mosaic.nc --topog_type box_channel --jwest_south 5 --jwest_north 15 --jeast_south 30 --jeast_north 45 --bottom_depth 1200
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for channel_grid"
   exit 1
endif

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 3: ocean_grid for box channel experiment                *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 4: ocean_grid for torus experiment  
#
#**********************************************************************************
set workdir = $basedir/$torus_grid
set output  = $outdir/$torus_grid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#create hgrid
make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,20 --ybnd 10,30 --dlon 1,1 --dlat 1,1 --grid_name ocean_hgrid
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for torus_grid"
   exit 1
endif

#create vgrid
make_vgrid --nbnds 2 --bnds 0.,200 --dbnds 10,10 --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for torus_grid"
   exit 1
endif

#create solo mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc --periodx 20 --periody 20
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for torus_grid"
   exit 1
endif

#create topography data
make_topog --mosaic ocean_mosaic.nc --topog_type rectangular_basin --bottom_depth 200
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for torus_grid"
   exit 1
endif

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 4: ocean_grid for torus experiment                      *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 5: ocean_grid for gyre experiment   
#
#**********************************************************************************
set workdir = $basedir/$gyre_grid
set output  = $outdir/$gyre_grid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#create hgrid
make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,10 --ybnd 15,35 --dlon 2,2 --dlat 2,2 --grid_name ocean_hgrid
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for gyre_grid"
   exit 1
endif

#create vgrid
make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for gyre_grid"
   exit 1
endif

#create solo mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for gyre_grid"
   exit 1
endif

#create topography data
make_topog --mosaic ocean_mosaic.nc --topog_type rectangular_basin --bottom_depth 5500
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for gyre_grid"
   exit 1
endif

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 5: ocean_grid for gyre experiment                       *"
echo "*                                                                         *"
echo "***************************************************************************"


#***********************************************************************************
#
#    Test 6: grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid.
#
#**********************************************************************************
set workdir = $basedir/$coupled_grid
set output  = $outdir/$coupled_grid
set input   = $indir
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#create ocean_hgrid 
make_hgrid --grid_type tripolar_grid --nxbnd 2 --nybnd 7 --xbnd -280,80 --ybnd -82,-30,-10,0,10,30,90 --dlon 1.0,1.0 --dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 --grid_name ocean_hgrid --center c_cell
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for coupled_grid ocean"
   exit 1
endif

#create ocean_vgrid
make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for coupled_grid ocean"
   exit 1
endif

#create ocean solo mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc --periodx 360
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for coupled_grid ocean"
   exit 1
endif

#create ocean topography data
make_topog --mosaic ocean_mosaic.nc --topog_type realistic --topog_file $input/OCCAM_p5degree.nc --topog_field TOPO --scale_factor -1 --vgrid ocean_vgrid.nc --output topog.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for coupled_grid ocean"
   exit 1
endif

if( ! -d parallel ) mkdir parallel
cd parallel

aprun -n $npes make_topog_parallel --mosaic ../ocean_mosaic.nc --topog_type realistic --topog_file $input/OCCAM_p5degree.nc --topog_field TOPO --scale_factor -1 --vgrid ../ocean_vgrid.nc --output topog.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_topog_parallel for coupled_grid ocean"
   exit 1
endif

nccmp -md topog.nc ../topog.nc
if( $status != 0 ) then
   echo "ERROR: make_topog did not reproduce from processor count for coupled_grid"
   exit 1
endif

cd ..

#Create C48 grid for atmos/land
make_hgrid --grid_type gnomonic_ed --nlon 96 --grid_name C48_grid
if( $status != 0 ) then
    echo "ERROR:  failed at make_hgrid for C48_grid"
    exit 1
endif
#create C48 solo mosaic for atmos/land
make_solo_mosaic --num_tiles 6 --dir ./ --mosaic C48_mosaic --tile_file C48_grid.tile1.nc,C48_grid.tile2.nc,C48_grid.tile3.nc,C48_grid.tile4.nc,C48_grid.tile5.nc,C48_grid.tile6.nc
if( $status != 0 ) then
    echo "ERROR: failed at make_solo_mosaic for C48_grid"
    exit 1
endif

#make the coupler_mosaic
make_coupler_mosaic --atmos_mosaic C48_mosaic.nc --ocean_mosaic ocean_mosaic.nc --ocean_topog  topog.nc --mosaic_name grid_spec --check --area_ratio_thresh 1.e-10 --check 
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_coupler_mosaic for coupled_grid"
    exit 1
endif

#check reproducing ability between processor count for make_coupler_mosaic
cd parallel
aprun -n $npes make_coupler_mosaic_parallel --atmos_mosaic ../C48_mosaic.nc --ocean_mosaic ../ocean_mosaic.nc --ocean_topog  ../topog.nc --mosaic_name grid_spec  --area_ratio_thresh 1.e-10 
  
if( $status != 0 ) then
    echo "ERROR: failed at make_coupler_mosaic_parallel for coupled_grid"
    exit 1
endif
foreach file ( `ls -1 *.nc` )
   nccmp -md $file ../$file
   if( $status != 0 ) then
      echo "ERROR: make_coupler_mosaic could not reproduce for $file for different processor count"
      exit 1
   endif
endif

cd ..
#copy data to the output directory
mv *.nc $output 


echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 6: grid for coupled model                               *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 7: grid for coupled nest model (land are C48 and ocean is 1 degree tripolar grid, atmosphere is C48 with nested region
#
#**********************************************************************************
set workdir = $basedir/$coupled_nest_grid
set output  = $outdir/$coupled_nest_grid
set input   = $indir
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#create ocean_hgrid 
make_hgrid --grid_type tripolar_grid --nxbnd 2 --nybnd 7 --xbnd -280,80 --ybnd -82,-30,-10,0,10,30,90 --dlon 1.0,1.0 --dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 --grid_name ocean_hgrid --center c_cell
if( $status != 0 ) then
   echo "ERROR: failed at make_hgrid for coupled_nest_grid ocean"
   exit 1
endif

#create ocean_vgrid
make_vgrid --nbnds 3 --bnds 0.,220.,5500. --dbnds 10.,10.,367.14286 --center c_cell --grid_name ocean_vgrid 
if( $status != 0 ) then
   echo "ERROR: failed at make_vgrid for coupled_nest_grid ocean"
   exit 1
endif

#create ocean solo mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc --periodx 360
if( $status != 0 ) then
   echo "ERROR: failed at make_solo_mosaic for coupled_nest_grid ocean"
   exit 1
endif

#create ocean topography data
make_topog --mosaic ocean_mosaic.nc --topog_type realistic --topog_file $input/OCCAM_p5degree.nc \
           --topog_field TOPO --scale_factor -1 --vgrid ocean_vgrid.nc --output topog.nc
if( $status != 0 ) then
   echo "ERROR: failed at make_topog for coupled_nest_grid ocean"
   exit 1
endif

#Create C48 grid with atmos nested grid.
make_hgrid --grid_type gnomonic_ed --nlon 96 --grid_name atmos_grid --do_schmidt \
           --target_lat 48.15 --target_lon -100.15 --halo 3 --stretch_factor 3 \
	   --great_circle_algorithm --nest_grid --refine_ratio 3 --parent_tile 4 \
	   --istart_nest 21 --iend_nest 60 --jstart_nest 11 --jend 70
if( $status != 0 ) then
    echo "ERROR:  failed at make_hgrid for nested atmos_grid"
    exit 1
endif
#create C48 solo mosaic for atmos
make_solo_mosaic --num_tiles 7 --dir ./ --mosaic atmos_mosaic --tile_file atmos_grid.tile1.nc,atmos_grid.tile2.nc,atmos_grid.tile3.nc,atmos_grid.tile4.nc,atmos_grid.tile5.nc,atmos_grid.tile6.nc,atmos_grid.tile7.nc
if( $status != 0 ) then
    echo "ERROR: failed at make_solo_mosaic for nested atmosphere grid"
    exit 1
endif

#Create C144 grid for land
make_hgrid --grid_type gnomonic_ed --nlon 288 --grid_name land_grid --do_schmidt \
           --target_lat 48.15 --target_lon -100.15 --halo 3 --stretch_factor 3   \
	   --great_circle_algorithm --nest_grid --refine_ratio 3 --parent_tile 0
if( $status != 0 ) then
    echo "ERROR:  failed at make_hgrid for nested C48_grid"
    exit 1
endif
#create C144 solo mosaic for land
make_solo_mosaic --num_tiles 6 --dir ./ --mosaic land_mosaic \
  --tile_file land_grid.tile1.nc,land_grid.tile2.nc,land_grid.tile3.nc,land_grid.tile4.nc,land_grid.tile5.nc,land_grid.tile6.nc

if( $status != 0 ) then
    echo "ERROR: failed at make_solo_mosaic for land"
    exit 1
endif

#make the coupler_mosaic
aprun -n $npes2 make_coupler_mosaic_parallel --atmos_mosaic atmos_mosaic.nc --land_mosaic land_mosaic.nc \
          --ocean_mosaic ocean_mosaic.nc --ocean_topog  topog.nc --interp_order 1 --mosaic_name grid_spec --check
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_coupler_mosaic for coupled_nest_grid"
    exit 1
endif

#check reproducing ability between processor count for make_coupler_mosaic
if( ! -d parallel ) mkdir parallel
cd parallel
aprun -n $npes make_coupler_mosaic_parallel --atmos_mosaic ../atmos_mosaic.nc --land_mosaic ../land_mosaic.nc \
         --ocean_mosaic ../ocean_mosaic.nc --ocean_topog  ../topog.nc --interp_order 1 --mosaic_name grid_spec

if( $status != 0 ) then
    echo "ERROR: failed at make_coupler_mosaic on 10 pe for coupled_nest_grid"
    exit 1
endif
foreach file ( `ls -1 *.nc` )
   nccmp -md $file ../$file
   if( $status != 0 ) then
      echo "ERROR: make_coupler_mosaic could not reproduce for $file from processor count for nested grid"
      exit 1
   endif
endif

cd ..

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 7: grid for coupled model with nested grid              *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 8: remap data from C48 to regular lat-lon grid
#
#**********************************************************************************
set workdir = $basedir/$C48_to_N45
set output  = $outdir/$C48_to_N45
set input   = $indir/$C48_to_N45
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

fregrid --input_mosaic $input/C48_mosaic.nc --input_dir $input --input_file 19800101.atmos_daily \
        --scalar_field zsurf,temp,t_surf --nlon 144 --nlat 90 --interp_method conserve_order2 \
	--output_file 19800101.atmos_daily.nc --check_conserve --remap_file C48_to_N45_remap.nc

if( $status != 0 ) then
    unset echo
    echo "ERROR: run failed for fregrid from C48 to N45 using conservative_order2"
    exit 1
endif

if( ! -d parallel ) mkdir parallel
cd parallel

#running on multiple processors
aprun -n $npes fregrid_parallel --input_mosaic $input/C48_mosaic.nc --input_dir $input --input_file 19800101.atmos_daily \
                                --scalar_field zsurf,temp,t_surf --nlon 144 --nlat 90 --interp_method conserve_order2    \
				--output_file 19800101.atmos_daily.nc --check_conserve

if( $status != 0 ) then
    unset echo
    echo "ERROR: failed for fregrid from C48 to N45 using conservative_order2 on $npes pe"
    exit 1
endif

nccmp -md 19800101.atmos_daily.nc ../19800101.atmos_daily.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: fregrid could not reproduce between 1 and $npes processors"
    exit 1
endif

#test reading from remap file
aprun -n $npes fregrid_parallel --input_mosaic $input/C48_mosaic.nc --input_dir $input --input_file 19800101.atmos_daily \
                                --scalar_field zsurf,temp,t_surf --nlon 144 --nlat 90 --interp_method conserve_order2 \
				--output_file 19800101.atmos_daily.read.nc --check_conserve --remap_file ../C48_to_N45_remap.nc

if( $status != 0 ) then
    unset echo
    echo "ERROR: failed for fregrid from C48 to N45 using conservative_order2 on $npes pe and read remap file"
    exit 1
endif

nccmp -md 19800101.atmos_daily.read.nc ../19800101.atmos_daily.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: fregrid could not reproduce when reading from remap file"
    exit 1
endif

cd ..
#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 8: remap data from C48 to lat-lon                       *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 9: remap cmip5 ocean data onto regular lat-lon grid (use GFDL-CM3 data as example)
#
#**********************************************************************************
set workdir = $basedir/$cmip5_regrid
set output  = $outdir/$cmip5_regrid
set input   = $indir/$cmip5_regrid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

# create GFDL-CM3_grid.nc
make_hgrid --grid_type from_file --my_grid_file $input/GFDL-CM3_sic_OImon.nc --grid_name GFDL-CM3_grid
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_hgrid for GFDL-CM3 cmip5 ocean data"
    exit 1
endif

#. create GFDL-CM3_mosaic.nc
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name GFDL-CM3_mosaic --tile_file GFDL-CM3_grid.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_solo_mosaic for GFDL-CM3 cmip5 ocean data"
    exit 1
endif

#. create remap file 
fregrid --input_mosaic GFDL-CM3_mosaic.nc --nlon 180 --nlat 90 --remap_file remap_file.GFDL-CM3.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at creating remap_file for GFDL-CM3 cmip5 ocean data"
    exit 1
endif

#. remap data
fregrid --input_mosaic GFDL-CM3_mosaic.nc --nlon 180 --nlat 90 --remap_file remap_file.GFDL-CM3.nc --input_dir $input \
        --input_file GFDL-CM3_sic_OImon.nc --output_file GFDL-sic_OImon.out.nc --scalar_field sic --check_conserve
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at regrid for GFDL-CM3 cmip5 ocean data"
    exit 1
endif

#Create MPI-ESM-LR_grid.nc
make_hgrid --grid_type from_file --my_grid_file $input/MIROC-ESM_sic_OImon.nc --grid_name MIROC-ESM_grid
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_hgrid for MIROC-ESM  cmip5 ocean data"
    exit 1
endif


#create MPI-ESM-LR_mosaic.nc
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name MIROC-ESM_mosaic --tile_file MIROC-ESM_grid.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_solo_mosaic for MIROC-ESM cmip5 ocean data"
    exit 1
endif

#create remap file
fregrid --input_mosaic MIROC-ESM_mosaic.nc --nlon 180 --nlat 90 --remap_file remap_file.MIROC-ESM.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_solo_mosaic for MIROC-ESM cmip5 ocean data"
    exit 1
endif

#remap data
fregrid --input_mosaic MIROC-ESM_mosaic.nc --nlon 180 --nlat 90 --remap_file remap_file.MIROC-ESM.nc --input_dir $input \
        --input_file MIROC-ESM_sic_OImon.nc --output_file MIROC-ESM-sic_OImon.out.nc --scalar_field sic --check_conserve
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at make_solo_mosaic for MIROC-ESM cmip5 ocean data"
    exit 1
endif

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 9: remap cmip5 ocean data                               *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 10: remap data onto cm2m ocean grid with extrapolation and vertical interpolation.
#
#**********************************************************************************
set workdir = $basedir/$regrid_extrap
set output  = $outdir/$regrid_extrap
set input   = $indir/$regrid_extrap
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,360 --ybnd -90,90 --nlon 720 --nlat 360 --grid_name levitus_grid
if( $status != 0 ) then
    echo "ERROR: failed at make_hgrid for levitus grid"
    exit 1 
endif 
make_solo_mosaic --num_tiles 1 --dir . --mosaic_name levitus_mosaic --tile_file levitus_grid.nc --periodx 360
if( $status != 0 ) then
    echo "ERROR: failed at make_solo_mosaic for levitus grid"
    exit 1 
endif 

fregrid --input_mosaic levitus_mosaic.nc --input_dir $input --input_file WOA09_ann_theta.nc \
        --scalar_field POTENTIAL_TEMP --output_file WOA09_ann_theta_cm2g_extrap.nc --output_mosaic \
	 $input/ocean_mosaic.nc --extrapolate --dst_vgrid $input/ocean_vgrid.nc --check_conserve
if( $status != 0 ) then
    echo "ERROR: failed at regrid data onto cm2m ocean grid with extrapolation and vertical interpolation."
    exit 1 
endif 

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 10: remap data onto cm2m grid with                      *"
echo "*   extrpolation and vertical interpolation                               *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 11: remap runoff data from regular lat-lon grid onto cm2m grid
#
#**********************************************************************************
set workdir = $basedir/$runoff_regrid
set output  = $outdir/$runoff_regrid
set input   = $indir/$runoff_regrid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

runoff_regrid --input_file $input/runoff.daitren.iaf.nc --input_fld_name runoff --output_mosaic $input/ocean_mosaic.nc \
              --output_topog $input/topog.nc --output_file runoff.cm2m.nc
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at runoff_regrid from lat-lon grid onto cm2m ocean grid"
    exit 1
endif

#copy data to the output directory
mv *.nc $output

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 11: remap runoff data onto cm2m grid                    *"
echo "*                                                                         *"
echo "***************************************************************************"


#***********************************************************************************
#
#    Test 12: river_regrid: remap data from lat-lon onto C48 grid.
#
#**********************************************************************************
set workdir = $basedir/$river_regrid
set output  = $outdir/$river_regrid
set input   = $indir/$river_regrid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

river_regrid --mosaic $input/grid_spec.nc --river_src $input/z1l_river_output_M45_tripolar_aug24.nc \
             --output river_data_C48 --land_thresh 0.000001
if( $status != 0 ) then 
    unset echo 
    echo "ERROR: failed for river_regrid to regrid river data onto C48 mosaic grid."
    exit 1 
endif

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 12: remap river_regrid onto C48 grid                    *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 13: check_mask: create mask_table to mask out some all-land domain
#                         to save processor usage
# 
#             for a sea-ice model, baltic1 experiment
#**********************************************************************************
set workdir = $basedir/$checkmask_baltic1
set output  = $outdir/$checkmask_baltic1
set input   = $indir/$checkmask_baltic1
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

check_mask --grid_file $input/baltic1_grid_spec.nc --min_pe 60 --max_pe 80 

if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at check_mask for baltic1 experiment"
    exit 1
endif

#copy data to the output directory
mv mask_table* $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 13:  check_mask for baltic1 experiment                  *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 14: check_mask: create mask_table to mask out some all-land domain
#                         to save processor usage
# 
#             for a coupled model with cm2m grid.
#**********************************************************************************
set workdir = $basedir/$checkmask_MOM6
set output  = $outdir/$checkmask_MOM6
set input   = $indir/$checkmask_MOM6
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

check_mask --grid_file $input/ocean_mosaic.nc --ocean_topog $input/topog.nc --layout 45,72
if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at check_mask for a coupled model grid"
    exit 1
endif

mv mask_table* $output

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 14: check_mask for MOM6 grid                            *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 15: remap_land: remap land restart file.
#
#**********************************************************************************
set workdir = $basedir/$remap_land
set output  = $outdir/$remap_land
set input   = $indir/$remap_land
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

# Create remapping file and remap land.res file.
remap_land --file_type land --src_mosaic $input/C48_mosaic/C48_mosaic.nc \
           --dst_mosaic $input/C192_mosaic/C192_mosaic.nc                \
	   --src_restart $input/src_restart/land.res  \
	   --dst_restart land.res    \
	   --dst_cold_restart $input/dst_cold_restart/land.res  \
	   --remap_file remap_file_C48_to_C192 --print_memory

if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at remap_land for land"
    exit 1
endif


# Remap soil, snow, cana, glac, lake, vegn1, vegn2
foreach type ( soil snow cana glac lake vegn1 vegn2 )
 if( $type == "vegn1" || $type == "vegn2" ) then
   set filetype = vegn
 else
   set filetype = $type
 endif
 remap_land --file_type $filetype --src_mosaic $input/C48_mosaic/C48_mosaic.nc \
            --dst_mosaic $input/C192_mosaic/C192_mosaic.nc \
	    --src_restart $input/src_restart/$type.res \
	    --dst_restart $type.res  \
	    --dst_cold_restart $input/dst_cold_restart/$type.res \
	    --land_src_restart $input/src_restart/land.res \
	    --land_cold_restart $input/dst_cold_restart/land.res \
	    --remap_file remap_file_C48_to_C192 --print_memory
   if( $status != 0 ) then
      unset echo
      echo "ERROR: failed at remap_land for $type"
      exit 1
   endif
end

#Remap static_vegn will added in the future.
if( ! -d parallel ) mkdir parallel
cd parallel
aprun -n $npes remap_land_parallel --file_type land --src_mosaic $input/C48_mosaic/C48_mosaic.nc \
                                   --dst_mosaic $input/C192_mosaic/C192_mosaic.nc \
				   --src_restart $input/src_restart/land.res \
				   --dst_restart land.res \
				   --dst_cold_restart $input/dst_cold_restart/land.res \
				   --remap_file remap_file_C48_to_C192 --print_memory

if( $status != 0 ) then
    unset echo
    echo "ERROR: failed at remap_land_parallel for land"
    exit 1
endif

foreach file (`ls -1 *.nc`)
   nccmp -md $file ../$file
   if( $status != 0 ) then
      echo "ERROR: remap_land could not reproduce  between processor count for $file"
      exit 1
   endif
end
cd ..

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 15: remap land restart from C48 to C192                 *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 16: make_regional_mosaic: create mosaic file for regional output and do remapping
#             for regional output.
#
#**********************************************************************************
set workdir = $basedir/$regional_regrid
set output  = $outdir/$regional_regrid
set input   = $indir/$regional_regrid
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#Create tile grids for four regions.
make_regional_mosaic --global_mosaic $input/C192-rot-120.-0.-r1_mosaic.nc --regional_file $input/rregionatmos_4xdaily_eq_avg.tile2.nc
make_regional_mosaic --global_mosaic $input/C192-rot-120.-0.-r1_mosaic.nc --regional_file $input/rregionatmos_4xdaily_eq_avg.tile3.nc
make_regional_mosaic --global_mosaic $input/C192-rot-120.-0.-r1_mosaic.nc --regional_file $input/rregionatmos_4xdaily_eq_avg.tile5.nc
make_regional_mosaic --global_mosaic $input/C192-rot-120.-0.-r1_mosaic.nc --regional_file $input/rregionatmos_4xdaily_eq_avg.tile6.nc

#Create regional mosaic
make_solo_mosaic --num_tiles 4 --dir ./ --mosaic regional_mosaic --tile_file regional_grid.tile2.nc,regional_grid.tile3.nc,regional_grid.tile5.nc,regional_grid.tile6.nc

#Create regular lat-lon grid (100:160, -15:15, size is 360x180)
make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2 --xbnd 0,360 --ybnd -2,2 --nlon 720 --nlat 10 --grid_name latlon_grid

#Create lat-lon mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name latlon_mosaic --tile_file latlon_grid.nc

#Remap data
fregrid --input_mosaic regional_mosaic.nc --output_mosaic latlon_mosaic.nc --input_dir $input \
        --input_file rregionatmos_4xdaily_eq_avg --output_file atmos_4xdaily_latlon.nc --scalar_field ps --check_conserve


#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 16: using fregrid for regional output                   *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 17: mppncscatter first use mppncscatter to split the file
#             and then use mppnccombine to combine the file. The output should reproduce
#             original file
#
#**********************************************************************************
set workdir = $basedir/$mppncscatter
set output  = $outdir/$mppncscatter
set input   = $indir/$mppncscatter
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

foreach tile ( 1 2 3 4 5 6 )
   set file = fv_core.res.tile$tile.nc
   mppncscatter -i 2 -j 3 -x 2 -y 12  $input/$file
   mppnccombine -64 $file $file.???? 
   nccmp -md $file $input/$file
end

set file = ocean_temp_salt.res.nc
mppncscatter  -i 7 -j 7 -x 21 -y 14  $input/$file
mppnccombine -64 $file $file.???? 
nccmp -md $file $input/$file

set file = ice_model.res.nc
mppncscatter  -i 7 -j 7 -x 21 -y 14  $input/$file
mppnccombine -64 $file $file.???? 
nccmp -md $file $input/$file

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 17: mppncscatter                                        *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 18: make_quick_mosaic
#
#***********************************************************************************
set workdir = $basedir/$make_quick_mosaic
set output  = $outdir/$make_quick_mosaic
set input   = $indir/$make_quick_mosaic
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

make_quick_mosaic --input_mosaic $input/ocean_mosaic.nc --ocean_topog $input/ocean_topog.nc

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 18: make_quick_mosaic                                   *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 19: transfer to mosaic
#
#**********************************************************************************
set workdir = $basedir/$transfer_mosaic
set output  = $outdir/$transfer_mosaic
set input   = $indir/$transfer_mosaic
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

transfer_to_mosaic_grid --input_file  $input/grid_spec_v7_transfer.nc

#copy data to the output directory
mv *.nc $output 

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 19: transfer_to_mosaic_grid                             *"
echo "*                                                                         *"
echo "***************************************************************************"

#***********************************************************************************
#
#    Test 20: remap ocean restart file from CM2.1 to CM2.5
#
#**********************************************************************************
set workdir = $basedir/$CM21_to_CM25
set output  = $outdir/$CM21_to_CM25
set input   = $indir/$CM21_to_CM25
if( ! -d $workdir ) mkdir -p $workdir
if( ! -d $output  ) mkdir -p $output
cd $workdir

#Create a OM3-like regular lat-lon grid.
make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 7 --xbnd -280,80 --ybnd -82,-30,-10,0,10,30,90 \
           --dlon 1.0,1.0 --dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 --grid_name latlon_grid --center c_cell

#Create the lat-lon solo mosaic
make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name latlon_mosaic --tile_file latlon_grid.nc --periodx 360

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
aprun -n $npes fregrid_parallel --input_mosaic $input/CM2.1_mosaic.nc --input_dir $input \
                                --input_file ocean_temp_salt.res.nc --scalar_field temp,salt \
				--output_file ocean_temp_salt.res.latlon.nc --output_mosaic latlon_mosaic.nc --check_conserve

#Remap data from regular lat-lon grid onto CM2.5 ocean grid with extrapolation.
aprun -n $npes fregrid_parallel --input_mosaic latlon_mosaic.nc --input_dir ./ \
                                --input_file ocean_temp_salt.res.latlon.nc --scalar_field temp,salt \
				--output_file ocean_temp_salt.res.CM2.5.nc --output_mosaic $input/CM2.5_mosaic.nc --extrapolate  --check_conserve

#copy data to the output directory
mv *.nc $output

echo "***************************************************************************"
echo "*                                                                         *"
echo "*   complete Test 20: remap restart file from CM2.1 grid to CM2.5 grid    *"
echo "*                                                                         *"
echo "***************************************************************************"
