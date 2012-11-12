# Build the Ocean library

set srcList = ( mom5/ocean_bgc mom5/ocean_core mom5/ocean_diag mom5/ocean_blobs mom5/ocean_param/neutral mom5/ocean_param/sources mom5/ocean_param/lateral mom5/ocean_param/vertical mom5/ocean_param/gotm-4.0/include mom5/ocean_param/gotm-4.0/turbulence mom5/ocean_param/gotm-4.0/util mom5/ocean_tracers mom5/ocean_wave ocean_shared/generic_tracers )

set lib_name = "lib_ocean"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS" $srcList $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
