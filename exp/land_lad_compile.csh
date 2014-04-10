# Build the land_lad library

set srcList = ( land_param land_lad/vegetation land_lad/soil land_lad )

set lib_name = "lib_land_lad"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs"  -o "-I$executable:h:h/lib_FMS" $srcList $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
