# Build the Ocean and Ice Components library
# The list of source files that should be compiled for this component.

set srcList = ( mom5/ocean_bgc mom5/ocean_core mom5/ocean_diag mom5/ocean_blobs mom5/ocean_param/neutral mom5/ocean_param/sources mom5/ocean_param/lateral mom5/ocean_param/vertical mom5/ocean_param/gotm-4.0/include mom5/ocean_param/gotm-4.0/turbulence mom5/ocean_param/gotm-4.0/util mom5/ocean_tracers mom5/ocean_wave ocean_shared/generic_tracers  ice_sis ice_param )


# setup directory structure
  if ( ! -d $executable:h )              mkdir -p $executable:h
  if ( ! -d $executable:h:h/lib_ocean )  mkdir -p $executable:h:h/lib_ocean

# compile libs
set makeFile      = Make_lib_ocean
cd $executable:h:h/lib_ocean

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_ocean.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS" $srcList $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_ocean.a"
    exit 1
endif
