# Build the land_lad library

set pathnames_land_lad  = $code_dir/path_names_land_lad        # path to file containing list of source paths
##LAND Component
# The list of source files that should be compiled for this component.
cat > $pathnames_land_lad <<EOF_land_lad
land_lad/land_model.F90
land_lad/land_types.F90
land_lad/numerics.F90
land_lad/soil/land_properties.F90
land_lad/soil/rivers.F90
land_lad/soil/soil.F90
land_lad/vegetation/vegetation.F90
land_param/climap_albedo.F90
EOF_land_lad

# setup directory structure
mkdir -p $executable:h:h/lib_land_lad

# compile libs
cd $executable:h:h/lib_land_lad
$mkmf_lib -p lib_land_lad.a -c "$cppDefs"  -o "-I$executable:h:h/lib_FMS" $pathnames_land_lad $lib_include_dirs

make

if( $status ) then
    echo "Make failed to create  lib_land_lad.a"
    exit 1
endif
