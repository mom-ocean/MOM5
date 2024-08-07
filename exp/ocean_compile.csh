# Build the Ocean library

set srcList = ( mom5/ocean_core mom5/ocean_diag mom5/ocean_wave mom5/ocean_blobs mom5/ocean_param/neutral mom5/ocean_param/sources mom5/ocean_param/lateral mom5/ocean_param/vertical mom5/ocean_param/gotm-4.0/include mom5/ocean_param/gotm-4.0/turbulence mom5/ocean_param/gotm-4.0/util mom5/ocean_tracers )

set lib_name = "lib_ocean"

if( $type == ACCESS-OM || $type == ACCESS-CM || $type == ACCESS-OM-BGC || $type == ACCESS-ESM) then
    set srcList = ( $srcList mom5/ocean_access )
    if( $type ==  ACCESS-OM-BGC || $type == ACCESS-ESM) then
        set srcList = ( $srcList mom5/ocean_csiro_bgc )
    else if ( $type ==  ACCESS-OM ) then
        set srcList = ( $srcList mom5/ocean_bgc access/generic_tracers/generic_tracers access/generic_tracers/mocsy/src )
    endif
    mkdir -p $executable:h:h/$type/$lib_name
    cd $executable:h:h/$type/$lib_name
else
    set srcList = ( $srcList mom5/ocean_bgc ocean_shared/generic_tracers )
    mkdir -p $executable:h:h/$lib_name
    cd $executable:h:h/$lib_name
endif

$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "$includes" $srcList $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
