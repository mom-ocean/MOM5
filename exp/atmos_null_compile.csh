# Build the null atmos library

set srcList = ( atmos_null atmos_param/diag_integral atmos_param/monin_obukhov )

set lib_name = "lib_atmos_null"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS" $srcList $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
