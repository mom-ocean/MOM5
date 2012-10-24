# Build the Ice library
# The list of source files that should be compiled for this component.

set srcList = ( ice_sis ice_param )

set lib_name = "lib_ice"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS" $srcList $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
