# Build the shared FMS component library
# The list of source files that should be compiled for this component.

set pathnames_shared  = $code_dir/path_names_shared        # path to file containing list of source paths

$root/bin/list_paths $root/src/shared
mv path_names $pathnames_shared

set srcList = ( shared/tridiagonal shared/tracer_manager shared/topography shared/time_manager shared/time_interp shared/station_data shared/sat_vapor_pres shared/random_numbers shared/platform shared/oda_tools shared/diag_manager shared/data_override shared/constants shared/column_diagnostics shared/axis_utils shared/astronomy shared/amip_interp shared/fft shared/field_manager shared/horiz_interp shared/include shared/memutils shared/mosaic shared/coupler shared/test shared/version)

set lib_name = "lib_FMS"

# Set up the version string, needed for FMS build.
# The version string is the git hash of the commit used to build the code.
# The version of an executable can be found with the following command:
# readelf -p .rodata <executable> | grep -A 1 'version.F90'
setenv GIT_CONFIG_NOGLOBAL 'yes'

set old_hash=`grep 'MOM_VERSION' ../src/shared/version/version.F90 | cut -d '"' -f 2 | cut -d '=' -f 2`
set new_hash=`git rev-parse HEAD`

if ( $old_hash != $new_hash ) then
    sed -e "s/{MOM_VERSION}/$new_hash/g" $code_dir/shared/version/version.F90.template > $code_dir/shared/version/version.F90
endif

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" $srcList $pathnames_shared $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
