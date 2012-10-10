# Build the land_lad2 library

##LAND Component
# The list of source files that should be compiled for this component.
set pathnames_land_lad2  = $code_dir/path_names_land_lad2        # path to file containing list of source paths

cat > $pathnames_land_lad2 <<EOF_land_lad2
land_lad2/canopy_air/cana_tile.F90  
land_lad2/canopy_air/canopy_air.F90  
land_lad2/land_constants.F90  
land_lad2/land_data.F90     
land_lad2/land_model.F90    
land_lad2/land_tile.F90     
land_lad2/glacier/glac_tile.F90   
land_lad2/glacier/glacier.F90  
land_lad2/lake/lake.F90       
land_lad2/lake/lake_tile.F90   
land_lad2/river/river.F90       
land_lad2/river/river_physics.F90  
land_lad2/river/river_type.F90      
land_lad2/shared/nf_utils/getput.inc     
land_lad2/shared/nf_utils/getput_compressed.inc 
land_lad2/shared/nf_utils/nf_utils.F90
land_lad2/shared/nf_utils/nfc.F90 
land_lad2/shared/nf_utils/nfu.F90
land_lad2/shared/debug.inc  
land_lad2/shared/land_debug.F90  
land_lad2/shared/land_io.F90
land_lad2/shared/land_numerics.F90
land_lad2/shared/land_tile_diag.F90
land_lad2/shared/land_tile_diag_buff.F90
land_lad2/shared/land_tile_diag_sel.F90
land_lad2/shared/land_tile_io.F90
land_lad2/shared/land_utils.F90
land_lad2/shared/sphum.F90
land_lad2/snow/snow.F90
land_lad2/snow/snow_tile.F90
land_lad2/soil/soil.F90
land_lad2/soil/soil_tile.F90
land_lad2/soil/uptake.F90
land_lad2/topo_rough/topo_rough.F90
land_lad2/transitions/transitions.F90
land_lad2/vegetation/read_remap_cohort_data.inc
land_lad2/vegetation/vegetation.F90
land_lad2/vegetation/vegn_cohort.F90
land_lad2/vegetation/vegn_cohort_io.F90
land_lad2/vegetation/vegn_cohort_io.inc
land_lad2/vegetation/vegn_data.F90
land_lad2/vegetation/vegn_disturbance.F90
land_lad2/vegetation/vegn_dynamics.F90
land_lad2/vegetation/vegn_harvesting.F90
land_lad2/vegetation/vegn_photosynthesis.F90
land_lad2/vegetation/vegn_radiation.F90
land_lad2/vegetation/vegn_static_override.F90
land_lad2/vegetation/vegn_tile.F90

EOF_land_lad2

set srcList = ( )

# setup directory structure
  if ( ! -d $executable:h )             mkdir -p $executable:h
  if ( ! -d $executable:h:h/lib_land_lad2 )    mkdir -p $executable:h:h/lib_land_lad2

# compile libs
set makeFile      = Make_lib_land_lad2
cd $executable:h:h/lib_land_lad2

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_land_lad2.a -c "-DUSE_LOG_DIAG_FIELD_INFO $cppDefs"  -o "-I$executable:h:h/lib_FMS" $srcList $pathnames_land_lad2 $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_land_lad2.a"
    exit 1
endif

