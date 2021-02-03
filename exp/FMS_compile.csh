# Build the shared FMS component library
# The list of source files that should be compiled for this component.

set pathnames_shared  = $code_dir/path_names_shared        # path to file containing list of source paths

$root/bin/list_paths $root/src/shared
mv path_names $pathnames_shared

set srcList = ( shared/tridiagonal shared/tracer_manager shared/topography shared/time_manager shared/time_interp shared/station_data shared/sat_vapor_pres shared/random_numbers shared/platform shared/oda_tools shared/diag_manager shared/data_override shared/constants shared/column_diagnostics shared/axis_utils shared/astronomy shared/amip_interp shared/fft shared/field_manager shared/horiz_interp shared/include shared/memutils shared/mosaic shared/coupler shared/test )

set lib_name = "lib_FMS"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" $srcList $pathnames_shared $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
