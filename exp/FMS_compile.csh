# Build the shared FMS component library
# The list of source files that should be compiled for this component.

set pathnames_shared  = $code_dir/path_names_shared        # path to file containing list of source paths

cat > $pathnames_shared <<EOF_shared
shared/drifters/cloud_interpolator.F90
shared/drifters/drifters.F90
shared/drifters/drifters_comm.F90
shared/drifters/drifters_compute_k.h
shared/drifters/drifters_core.F90
shared/drifters/drifters_input.F90
shared/drifters/drifters_io.F90
shared/drifters/drifters_push.h
shared/drifters/drifters_set_field.h
shared/drifters/fms_switches.h
shared/drifters/quicksort.F90
shared/exchange/stock_constants.F90
shared/exchange/xgrid.F90
shared/fms/fms.F90
shared/fms/fms_io.F90
shared/fms/read_data_2d.inc
shared/fms/read_data_3d.inc
shared/fms/read_data_4d.inc
shared/fms/write_data.inc
shared/mpp/mpp.F90
shared/mpp/mpp_data.F90
shared/mpp/mpp_domains.F90
shared/mpp/mpp_io.F90
shared/mpp/mpp_memutils.F90
shared/mpp/mpp_parameter.F90
shared/mpp/mpp_pset.F90
shared/mpp/mpp_utilities.F90
shared/mpp/nsclock.c
shared/mpp/threadloc.c
shared/mpp/include/mpp_chksum.h
shared/mpp/include/mpp_chksum_int.h
shared/mpp/include/mpp_chksum_scalar.h
shared/mpp/include/mpp_comm.inc
shared/mpp/include/mpp_comm_mpi.inc
shared/mpp/include/mpp_comm_nocomm.inc
shared/mpp/include/mpp_comm_sma.inc
shared/mpp/include/mpp_data_mpi.inc
shared/mpp/include/mpp_data_nocomm.inc
shared/mpp/include/mpp_data_sma.inc
shared/mpp/include/mpp_do_get_boundary.h
shared/mpp/include/mpp_do_global_field.h
shared/mpp/include/mpp_do_redistribute.h
shared/mpp/include/mpp_do_update.h
shared/mpp/include/mpp_do_updateV.h
shared/mpp/include/mpp_do_updateV_ad.h
shared/mpp/include/mpp_do_update_ad.h
shared/mpp/include/mpp_domains_comm.inc
shared/mpp/include/mpp_domains_define.inc
shared/mpp/include/mpp_domains_misc.inc
shared/mpp/include/mpp_domains_reduce.inc
shared/mpp/include/mpp_domains_util.inc
shared/mpp/include/mpp_error_a_a.h
shared/mpp/include/mpp_error_a_s.h
shared/mpp/include/mpp_error_s_a.h
shared/mpp/include/mpp_error_s_s.h
shared/mpp/include/mpp_get_boundary.h
shared/mpp/include/mpp_global_field.h
shared/mpp/include/mpp_global_reduce.h
shared/mpp/include/mpp_global_sum.h
shared/mpp/include/mpp_global_sum_ad.h
shared/mpp/include/mpp_global_sum_tl.h
shared/mpp/include/mpp_io_connect.inc
shared/mpp/include/mpp_io_misc.inc
shared/mpp/include/mpp_io_read.inc
shared/mpp/include/mpp_io_util.inc
shared/mpp/include/mpp_io_write.inc
shared/mpp/include/mpp_read_2Ddecomp.h
shared/mpp/include/mpp_reduce_mpi.h
shared/mpp/include/mpp_reduce_nocomm.h
shared/mpp/include/mpp_reduce_sma.h
shared/mpp/include/mpp_sum.inc
shared/mpp/include/mpp_sum_mpi.h
shared/mpp/include/mpp_sum_nocomm.h
shared/mpp/include/mpp_sum_sma.h
shared/mpp/include/mpp_transmit.inc
shared/mpp/include/mpp_transmit_mpi.h
shared/mpp/include/mpp_transmit_nocomm.h
shared/mpp/include/mpp_transmit_sma.h
shared/mpp/include/mpp_update_domains2D.h
shared/mpp/include/mpp_update_domains2D_ad.h
shared/mpp/include/mpp_util.inc
shared/mpp/include/mpp_util_mpi.inc
shared/mpp/include/mpp_util_nocomm.inc
shared/mpp/include/mpp_util_sma.inc
shared/mpp/include/mpp_write.h
shared/mpp/include/mpp_write_2Ddecomp.h
shared/mpp/include/system_clock.h
shared/test/dump_field.F90

EOF_shared

set srcList = ( shared/tridiagonal shared/tracer_manager shared/topography shared/time_manager shared/time_interp shared/station_data shared/sat_vapor_pres shared/random_numbers shared/platform shared/oda_tools shared/diag_manager shared/data_override shared/constants shared/column_diagnostics shared/axis_utils shared/astronomy shared/amip_interp shared/fft shared/field_manager shared/horiz_interp shared/include shared/memutils shared/mosaic shared/coupler shared/test)

set lib_name = "lib_FMS"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" $srcList $pathnames_shared $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
