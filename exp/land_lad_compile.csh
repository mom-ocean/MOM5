# Build the land_lad library

set srcList = ( land_param land_lad/vegetation land_lad/soil land_lad )

# setup directory structure
mkdir -p $executable:h:h/lib_land_lad

# compile libs
cd $executable:h:h/lib_land_lad
$mkmf_lib -p lib_land_lad.a -c "$cppDefs"  -o "-I$executable:h:h/lib_FMS" $srcList $lib_include_dirs

make

if( $status ) then
    echo "Make failed to create  lib_land_lad.a"
    exit 1
endif
