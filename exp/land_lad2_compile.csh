# Build the land_lad2 library

set srcList = ( land_lad2/vegetation land_lad2/transitions land_lad2/topo_rough land_lad2/soil land_lad2/snow land_lad2/shared land_lad2/shared/nf_utils land_lad2/river land_lad2/lake land_lad2/glacier land_lad2/canopy_air land_lad2 )

set lib_name = "lib_land_lad2"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "-DUSE_LOG_DIAG_FIELD_INFO $cppDefs"  -o "-I$executable:h:h/lib_FMS" $srcList $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif

