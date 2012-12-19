# Build the atmos_BG library

set srcList = ( atmos_bgrid/tools atmos_bgrid/model atmos_bgrid/driver/coupled atmos_coupled )

set lib_name = "lib_atmos_bg"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_atmos_phys" $srcList $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
