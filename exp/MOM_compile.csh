#!/bin/csh -f
# Minimal compile script for fully coupled model CM2M experiments

set echo
set platform      = ncrc2.intel   # A unique identifier for your platfo
                                  # This corresponds to the mkmf templates in $root/bin dir.
set type          = MOM_solo      # Type of the experiment
set help = 0

set argv = (`getopt -u -o h -l type: -l platform:  -l help  --  $*`)
while ("$argv[1]" != "--")
    switch ($argv[1])
        case --type:
                set type = $argv[2]; shift argv; breaksw
        case --platform:
                set platform = $argv[2]; shift argv; breaksw
        case --help:
                set help = 1;  breaksw
        case -h:
                set help = 1;  breaksw
    endsw
    shift argv
end
shift argv
if ( $help ) then
    echo "The optional arguments are:"
    echo "--type       followed by the type of the experiment, currently one of the following:"
    echo "             MOM_solo : solo ocean model"
    echo "             MOM_SIS  : ocean-seaice model"
    echo "             CM2M     : ocean-seaice-land-atmosphere coupled climate model"
    echo "             ESM2M    : ocean-seaice-land-atmosphere coupled climate model with biogeochemistry, EarthSystemModel"
    echo "             ICCM     : ocean-seaice-land-atmosphere coupled model"
    echo
    echo "--platform   followed by the platform name that has a corresponfing environ file in the ../bin dir, default is ncrc.intel"
    echo
    echo
    exit 0
endif

#
# User does not need to change anything below!
#
set root          = $cwd:h                            # The directory you created when you checkout
set code_dir      = $root/src                         # source code directory
set executable    = $root/exec/$platform/$type/fms_$type.x      # executable created after compilation
set mppnccombine  = $root/bin/mppnccombine.$platform  # path to executable mppnccombine
set mkmfTemplate  = $root/bin/mkmf.template.$platform # path to template for your platform
set mkmf          = $root/bin/mkmf                    # path to executable mkmf
set cppDefs  = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS" )
#On Altrix systems you may include "-Duse_shared_pointers -Duse_SGI_GSM" in cppDefs for perfomance.
#These are included in the GFDL configuration of the model.
  
set static        = 0              # 1 if you want static memory allocation, 0 for dynamic
if($static) then
  set executable = $root/exec/$platform/${type}_static/fms_$type.x 
  set cppDefs = "$cppDefs -DMOM_STATIC_ARRAYS -DNI_=360 -DNJ_=200 -DNK_=50 -DNI_LOCAL_=60 -DNJ_LOCAL_=50"
endif

#
# Users must ensure the correct environment file exists for their platform.
#
source $root/bin/environs.$platform  # environment variables and loadable modules

# setup directory structure
  if ( ! -d $executable:h )    mkdir -p $executable:h

#
# compile mppnccombine.c, needed only if $npes > 1
  if ( ! -f $mppnccombine ) then
    cc -O -o $mppnccombine -I/usr/local/include -L/usr/local/lib $code_dir/postprocessing/mppnccombine/mppnccombine.c -lnetcdf
  endif


# Build the shared FMS component library
# The list of source files that should be compiled for this component.

set pathnames_shared  = $code_dir/path_names_shared        # path to file containing list of source paths

cat > $pathnames_shared <<EOF_shared
shared/amip_interp/amip_interp.F90
shared/astronomy/astronomy.F90
shared/axis_utils/axis_utils.F90
shared/column_diagnostics/column_diagnostics.F90
shared/constants/constants.F90
shared/coupler/atmos_ocean_fluxes.F90
shared/coupler/coupler_types.F90
shared/coupler/ensemble_manager.F90
shared/data_override/data_override.F90
shared/diag_manager/diag_axis.F90
shared/diag_manager/diag_data.F90
shared/diag_manager/diag_grid.F90
shared/diag_manager/diag_manager.F90
shared/diag_manager/diag_output.F90
shared/diag_manager/diag_util.F90
shared/diag_manager/diag_table.F90
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
shared/fft/fft.F90
shared/fft/fft99.F90
shared/field_manager/field_manager.F90
shared/field_manager/fm_util.F90
shared/field_manager/parse.inc
shared/fms/fms.F90
shared/fms/fms_io.F90
shared/fms/read_data_2d.inc
shared/fms/read_data_3d.inc
shared/fms/read_data_4d.inc
shared/fms/write_data.inc
shared/horiz_interp/horiz_interp.F90
shared/horiz_interp/horiz_interp_bicubic.F90
shared/horiz_interp/horiz_interp_bilinear.F90
shared/horiz_interp/horiz_interp_conserve.F90
shared/horiz_interp/horiz_interp_spherical.F90
shared/horiz_interp/horiz_interp_type.F90
shared/include/fms_platform.h
shared/memutils/memuse.c
shared/memutils/memutils.F90
shared/mosaic/constant.h
shared/mosaic/create_xgrid.c
shared/mosaic/create_xgrid.h
shared/mosaic/gradient.F90
shared/mosaic/gradient_c2l.c
shared/mosaic/gradient_c2l.h
shared/mosaic/grid.F90
shared/mosaic/interp.c
shared/mosaic/interp.h
shared/mosaic/mosaic.F90
shared/mosaic/mosaic_util.c
shared/mosaic/mosaic_util.h
shared/mosaic/read_mosaic.c
shared/mosaic/read_mosaic.h
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
shared/oda_tools/oda_core.F90
shared/oda_tools/oda_types.F90
shared/oda_tools/write_ocean_data.F90
shared/oda_tools/xbt_drop_rate_adjust.f90
shared/platform/platform.F90
shared/random_numbers/MersenneTwister.F90
shared/random_numbers/random_numbers.F90
shared/sat_vapor_pres/sat_vapor_pres.F90
shared/sat_vapor_pres/sat_vapor_pres_k.F90
shared/station_data/station_data.F90
shared/time_interp/time_interp.F90
shared/time_interp/time_interp_external.F90
shared/time_manager/get_cal_time.F90
shared/time_manager/time_manager.F90
shared/topography/gaussian_topog.F90
shared/topography/topography.F90
shared/tracer_manager/tracer_manager.F90
shared/tridiagonal/tridiagonal.F90

atmos_param/diag_integral/diag_integral.F90
atmos_param/monin_obukhov/monin_obukhov.F90
atmos_param/monin_obukhov/monin_obukhov_interfaces.h
atmos_param/monin_obukhov/monin_obukhov_kernel.F90

EOF_shared

# setup directory structure
set executable_h = $executable:h

if ( ! -d $executable_h )             mkdir -p $executable_h
if ( ! -d $executable_h:h/lib_FMS )    mkdir -p $executable_h:h/lib_FMS

# compile libs
set makeFile      = Make_lib_FMS
cd $executable_h:h/lib_FMS

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_FMS.a -c "$cppDefs" $pathnames_shared $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_FMS.a"
    exit
endif

# Build the Ocean and Ice Components library
# The list of source files that should be compiled for this component.

set srcList = ( mom5/ocean_bgc mom5/ocean_core mom5/ocean_diag mom5/ocean_blobs mom5/ocean_param/neutral mom5/ocean_param/sources mom5/ocean_param/lateral mom5/ocean_param/vertical mom5/ocean_param/gotm-4.0/include mom5/ocean_param/gotm-4.0/turbulence mom5/ocean_param/gotm-4.0/util mom5/ocean_tracers mom5/ocean_wave ocean_shared/generic_tracers  ice_sis ice_param )


# setup directory structure
  if ( ! -d $executable_h )              mkdir -p $executable_h
  if ( ! -d $executable_h:h/lib_ocean )  mkdir -p $executable_h:h/lib_ocean

# compile libs
set makeFile      = Make_lib_ocean
cd $executable_h:h/lib_ocean

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_ocean.a -c "$cppDefs" -o "-I$executable_h:h/lib_FMS" $srcList $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_ocean.a"
    exit
endif

if( $type == MOM_solo ) then
#   Build the MOM_solo executable   
    set srcList = ( mom5/drivers )
    set makeFile      = Make_MOM_solo
    if ( ! -d $executable_h ) mkdir -p $executable_h
    cd $executable_h
    $mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p $executable:t -o "-I$executable_h:h/lib_FMS -I$executable_h:h/lib_ocean" -l "$executable_h:h/lib_ocean/lib_ocean.a $executable_h:h/lib_FMS/lib_FMS.a"  $srcList

    make -f $makeFile
    if( $status ) then
	echo "Make failed to create the MOM_solo executable  "
	exit
    endif    
    exit
endif

if( $type == MOM_SIS ) then
#   Build the MOM_SIS executable   
    set srcList = ( atmos_null land_null coupler )
    set makeFile      = Make_MOM_SIS
    if ( ! -d $executable_h ) mkdir -p $executable_h
    cd $executable_h
    $mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p $executable:t -o "-I$executable_h:h/lib_FMS -I$executable_h:h/lib_ocean" -l "$executable_h:h/lib_ocean/lib_ocean.a $executable_h:h/lib_FMS/lib_FMS.a"  $srcList

    make -f $makeFile
    if( $status ) then
	echo "Make failed to create the MOM_SIS executable  "
	exit
    endif    
    exit
endif


if( $type == CM2M | $type == ESM2M | $type == ICCM ) then
# Build the am3 library

set pathnames_atmos_param  = $code_dir/path_names_atmos_param        # path to file containing list of source paths

cat > $pathnames_atmos_param <<EOF_atmos_param
atmos_param/betts_miller/betts_miller.F90                                                                
atmos_param/betts_miller/bm_massflux.F90                                                                      
atmos_param/betts_miller/bm_omp.F90                                   
atmos_param/cg_drag/cg_drag.F90                                  
atmos_param/cloud_generator/betaDistribution.F90                 
atmos_param/cloud_generator/cloud_generator.F90                   
atmos_param/cloud_obs/cloud_obs.F90                              
atmos_param/cloud_rad/cloud_rad.F90                              
atmos_param/cloud_zonal/cloud_zonal.F90                          
atmos_param/clouds/clouds.F90                                     
atmos_param/clubb/CLUBB_driver_SCM.F90
atmos_param/clubb/MG_microp_3D.F90
atmos_param/cosp/MISR_simulator/MISR_simulator.F
atmos_param/cosp/cosp.F90
atmos_param/cosp/cosp_constants.F90
atmos_param/cosp/cosp_defs.H
atmos_param/cosp/cosp_diagnostics.F90
atmos_param/cosp/cosp_driver.F90
atmos_param/cosp/cosp_io.F90
atmos_param/cosp/cosp_isccp_simulator.F90
atmos_param/cosp/cosp_lidar.F90
atmos_param/cosp/cosp_misr_simulator.F90
atmos_param/cosp/cosp_modis_simulator.F90
atmos_param/cosp/cosp_rttov_simulator.F90
atmos_param/cosp/cosp_simulator.F90
atmos_param/cosp/cosp_stats.F90
atmos_param/cosp/cosp_types.F90
atmos_param/cosp/cosp_utils.F90
atmos_param/cosp/MODIS_simulator/modis_simulator.F90
atmos_param/cosp/actsim/lidar_simulator.F90
atmos_param/cosp/actsim/lmd_ipsl_stats.F90
atmos_param/cosp/icarus-scops-3.7/congvec.H
atmos_param/cosp/icarus-scops-3.7/congvec.inc
atmos_param/cosp/icarus-scops-3.7/icarus.F90
atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F
atmos_param/cosp/icarus-scops-3.7/scops.F
atmos_param/cosp/llnl/cosp_radar.F90
atmos_param/cosp/llnl/llnl_stats.F90
atmos_param/cosp/llnl/pf_to_mr.F
atmos_param/cosp/llnl/prec_scops.F
atmos_param/cosp/quickbeam/array_lib.F90
atmos_param/cosp/quickbeam/atmos_lib.F90
atmos_param/cosp/quickbeam/dsd.F90
atmos_param/cosp/quickbeam/format_input.F90
atmos_param/cosp/quickbeam/gases.F90
atmos_param/cosp/quickbeam/load_hydrometeor_classes.F90
atmos_param/cosp/quickbeam/load_mie_table.F90
atmos_param/cosp/quickbeam/math_lib.F90
atmos_param/cosp/quickbeam/mrgrnk.F90
atmos_param/cosp/quickbeam/optics_lib.F90
atmos_param/cosp/quickbeam/radar_simulator.F90
atmos_param/cosp/quickbeam/radar_simulator_types.F90
atmos_param/cosp/quickbeam/zeff.F90
atmos_param/cu_mo_trans/cu_mo_trans.F90                          
atmos_param/damping_driver/damping_driver.F90                       
atmos_param/diag_cloud/diag_cloud.F90                            
atmos_param/diag_cloud_rad/diag_cloud_rad.F90                                            
atmos_param/diffusivity/diffusivity.F90                           
atmos_param/donner_deep/cumulus_closure_k.F90                     
atmos_param/donner_deep/donner_cape_k.F90                        
atmos_param/donner_deep/donner_cloud_model_k.F90                 
atmos_param/donner_deep/donner_deep.F90                           
atmos_param/donner_deep/donner_deep_k.F90                         
atmos_param/donner_deep/donner_deep_miz.F90                       
atmos_param/donner_deep/donner_lite_k.F90                        
atmos_param/donner_deep/donner_lscloud_k.F90                     
atmos_param/donner_deep/donner_meso_k.F90                         
atmos_param/donner_deep/donner_nml.h                             
atmos_param/donner_deep/donner_rad_k.F90                         
atmos_param/donner_deep/donner_types.F90                     
atmos_param/donner_deep/donner_types.h                           
atmos_param/donner_deep/donner_utilities_k.F90                    
atmos_param/donner_deep/fms_donner.F90                         
atmos_param/donner_deep/nonfms_donner.F90                        
atmos_param/donner_deep/wet_deposition_0D.F90                    
atmos_param/dry_adj/dry_adj.F90                                  
atmos_param/edt/edt.F90                                          
atmos_param/entrain/entrain.F90                                   
atmos_param/fsrad/co2_data.F90                                    
atmos_param/fsrad/co2int.F90                                                                                             
atmos_param/fsrad/fs_profile.F90                                 
atmos_param/fsrad/fsrad.F90                                      
atmos_param/fsrad/hconst.F90                                     
atmos_param/fsrad/longwave.F90                                   
atmos_param/fsrad/mcm_lw.F90                                     
atmos_param/fsrad/mcm_sw_driver.F90                              
atmos_param/fsrad/mcm_swnew.F90                                   
atmos_param/fsrad/mcm_swtbls.F90                                 
atmos_param/fsrad/rad_diag.F90                                   
atmos_param/fsrad/rdparm.F90
atmos_param/fsrad/shortwave.F90
atmos_param/grey_radiation/grey_radiation.F90
atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
atmos_param/lscale_cond/lscale_cond.F90
atmos_param/mg_drag/mg_drag.F90
atmos_param/moist_conv/moist_conv.F90
atmos_param/moist_processes/detr_ice_num.F90
atmos_param/moist_processes/moist_processes.F90
atmos_param/moist_processes/moist_processes_utils.F90
atmos_param/moist_processes/moistproc_kernels.F90
atmos_param/my25_turb/my25_turb.F90
atmos_param/physics_driver/physics_driver.F90
atmos_param/radiation_driver/radiation_driver.F90
atmos_param/ras/ras.F90
atmos_param/rh_clouds/rh_clouds.F90
atmos_param/sea_esf_rad/aerosol.F90
atmos_param/sea_esf_rad/aerosolrad_package.F90
atmos_param/sea_esf_rad/bulkphys_rad.F90
atmos_param/sea_esf_rad/cloud_spec.F90
atmos_param/sea_esf_rad/cloudrad_diagnostics.F90
atmos_param/sea_esf_rad/cloudrad_package.F90
atmos_param/sea_esf_rad/diag_clouds_W.F90
atmos_param/sea_esf_rad/donner_deep_clouds_W.F90
atmos_param/sea_esf_rad/esfsw_driver.F90
atmos_param/sea_esf_rad/esfsw_parameters.F90
atmos_param/sea_esf_rad/gas_tf.F90
atmos_param/sea_esf_rad/isccp_clouds.F90
atmos_param/sea_esf_rad/lhsw_driver.F90
atmos_param/sea_esf_rad/longwave_clouds.F90
atmos_param/sea_esf_rad/longwave_driver.F90
atmos_param/sea_esf_rad/longwave_fluxes.F90
atmos_param/sea_esf_rad/longwave_params.F90
atmos_param/sea_esf_rad/longwave_tables.F90
atmos_param/sea_esf_rad/lw_gases_stdtf.F90
atmos_param/sea_esf_rad/mgrp_prscr_clds.F90
atmos_param/sea_esf_rad/microphys_cloud.F90
atmos_param/sea_esf_rad/microphys_rad.F90
atmos_param/sea_esf_rad/optical_path.F90
atmos_param/sea_esf_rad/original_fms_rad.F90
atmos_param/sea_esf_rad/ozone.F90
atmos_param/sea_esf_rad/rad_output_file.F90
atmos_param/sea_esf_rad/rad_utilities.F90
atmos_param/sea_esf_rad/radiation_diag.F90
atmos_param/sea_esf_rad/radiative_gases.F90
atmos_param/sea_esf_rad/rh_based_clouds.F90
atmos_param/sea_esf_rad/sea_esf_rad.F90
atmos_param/sea_esf_rad/sealw99.F90
atmos_param/sea_esf_rad/shortwave_driver.F90
atmos_param/sea_esf_rad/specified_clouds_W.F90
atmos_param/sea_esf_rad/standalone_clouds.F90
atmos_param/sea_esf_rad/strat_clouds_W.F90
atmos_param/sea_esf_rad/uw_clouds_W.F90
atmos_param/sea_esf_rad/zetac_clouds_W.F90
atmos_param/shallow_conv/shallow_conv.F90
atmos_param/shallow_cu/conv_closures.F90
atmos_param/shallow_cu/conv_plumes.F90
atmos_param/shallow_cu/conv_plumes_k.F90
atmos_param/shallow_cu/conv_utilities.F90
atmos_param/shallow_cu/conv_utilities_k.F90
atmos_param/shallow_cu/deep_conv.F90
atmos_param/shallow_cu/uw_conv.F90
atmos_param/stable_bl_turb/stable_bl_turb.F90
atmos_param/strat_cloud/aerosol_cloud.F90
atmos_param/strat_cloud/check_nan.F90
atmos_param/strat_cloud/cldwat2m_micro.F90
atmos_param/strat_cloud/gamma_mg.F90
atmos_param/strat_cloud/mg_const.F90
atmos_param/strat_cloud/microphysics.F90
atmos_param/strat_cloud/morrison_gettelman_microp.F90
atmos_param/strat_cloud/nc_cond.F90
atmos_param/strat_cloud/polysvp.F90
atmos_param/strat_cloud/rotstayn_klein_mp.F90
atmos_param/strat_cloud/simple_pdf.F90
atmos_param/strat_cloud/strat_cloud.F90
atmos_param/strat_cloud/strat_cloud_legacy.F90
atmos_param/strat_cloud/strat_cloud_utilities.F90
atmos_param/strat_cloud/strat_netcdf.F90
atmos_param/strat_cloud/strat_nml.h
atmos_param/topo_drag/topo_drag.F90
atmos_param/vert_diff/vert_diff.F90
atmos_param/vert_diff_driver/vert_diff_driver.F90
atmos_param/vert_turb_driver/vert_turb_driver.F90
atmos_shared/atmos_nudge/atmos_nudge.F90
atmos_shared/interpolator/interpolator.F90
atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90
atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90
atmos_shared/tracer_driver/atmos_age_tracer.F90
atmos_shared/tracer_driver/atmos_carbon_aerosol.F90
atmos_shared/tracer_driver/atmos_ch3i.F90
atmos_shared/tracer_driver/atmos_co2.F90
atmos_shared/tracer_driver/atmos_convection_tracer.F90
atmos_shared/tracer_driver/atmos_dust.F90
atmos_shared/tracer_driver/atmos_radon.F90
atmos_shared/tracer_driver/atmos_sea_salt.F90
atmos_shared/tracer_driver/atmos_soa.F90
atmos_shared/tracer_driver/atmos_sulfate.F90
atmos_shared/tracer_driver/atmos_sulfur_hex.F90
atmos_shared/tracer_driver/atmos_tracer_driver.F90
atmos_shared/tracer_driver/atmos_tracer_utilities.F90
atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90
atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
atmos_shared/tracer_driver/tropchem/m_tracname.F90
atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90
atmos_shared/tracer_driver/tropchem/mo_chemdr.F90
atmos_shared/tracer_driver/tropchem/mo_chemini.F90
atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90
atmos_shared/tracer_driver/tropchem/mo_fastjx.F90
atmos_shared/tracer_driver/tropchem/mo_fphoto.F90
atmos_shared/tracer_driver/tropchem/mo_hook.F90
atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90
atmos_shared/tracer_driver/tropchem/mo_jpl.F90
atmos_shared/tracer_driver/tropchem/mo_photo.F90
atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90
atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90
atmos_shared/tracer_driver/tropchem/mo_setinv.F90
atmos_shared/tracer_driver/tropchem/mo_setsox.F90
atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90
atmos_shared/tracer_driver/tropchem/moz.mat.F90
atmos_shared/tracer_driver/tropchem/moz.mods.F90
atmos_shared/tracer_driver/tropchem/moz.subs.F90
atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
atmos_shared/tracer_driver/tropchem/tropchem_driver.F90
atmos_shared/vert_advection/vert_advection.F90

EOF_atmos_param

# setup directory structure
  if ( ! -d $executable_h )             mkdir -p $executable_h
  if ( ! -d $executable_h:h/lib_atmos_phys )    mkdir -p $executable_h:h/lib_atmos_phys

# compile libs
set makeFile      = Make_lib_atmos_phys
cd $executable_h:h/lib_atmos_phys

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_atmos_phys.a -c "$cppDefs" -o "-I$executable_h:h/lib_FMS" $pathnames_atmos_param $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_atmos_phys.a"
    exit
endif

endif


if( $type == CM2M | $type == ESM2M ) then
# Build the atmos_FV library

set pathnames_atmos_fv  = $code_dir/path_names_atmos_fv        # path to file containing list of source paths

cat > $pathnames_atmos_fv <<EOF_atmos_fv
atmos_coupled/atmos_model.F90
atmos_fv_dynamics/driver/coupled/atmosphere.F90
atmos_fv_dynamics/driver/coupled/fv_physics.F90
atmos_fv_dynamics/model/dyn_core.F90
atmos_fv_dynamics/model/ecmfft.F90
atmos_fv_dynamics/model/fill_module.F90
atmos_fv_dynamics/model/fv_arrays.F90
atmos_fv_dynamics/model/fv_arrays.h
atmos_fv_dynamics/model/fv_dynamics.F90
atmos_fv_dynamics/model/fv_pack.F90
atmos_fv_dynamics/model/fv_point.inc
atmos_fv_dynamics/model/mapz_module.F90
atmos_fv_dynamics/model/pft_module.F90
atmos_fv_dynamics/model/shr_kind_mod.F90
atmos_fv_dynamics/model/sw_core.F90
atmos_fv_dynamics/model/tp_core.F90
atmos_fv_dynamics/model/tracer_2d.F90
atmos_fv_dynamics/model/update_fv_phys.F90
atmos_fv_dynamics/tools/age_of_air.F90
atmos_fv_dynamics/tools/fv_diagnostics.F90
atmos_fv_dynamics/tools/fv_restart.F90
atmos_fv_dynamics/tools/getmax.F90
atmos_fv_dynamics/tools/gmean.F90
atmos_fv_dynamics/tools/init_dry_atm.F90
atmos_fv_dynamics/tools/init_sw_ic.F90
atmos_fv_dynamics/tools/mod_comm.F90
atmos_fv_dynamics/tools/par_vecsum.F90
atmos_fv_dynamics/tools/pmaxmin.F90
atmos_fv_dynamics/tools/pv_module.F90
atmos_fv_dynamics/tools/set_eta.F90
atmos_fv_dynamics/tools/timingModule.F90
atmos_fv_dynamics/tools/upper.F90

EOF_atmos_fv

# setup directory structure
  if ( ! -d $executable_h )             mkdir -p $executable_h
  if ( ! -d $executable_h:h/lib_atmos_phys_fv )    mkdir -p $executable_h:h/lib_atmos_fv

# compile libs
set makeFile      = Make_lib_atmos_fv
cd $executable_h:h/lib_atmos_fv

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_atmos_fv.a -c "$cppDefs" -o "-I$executable_h:h/lib_FMS -I$executable_h:h/lib_atmos_phys" $pathnames_atmos_fv $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_atmos_fv.a"
    exit
endif

endif

if( $type == CM2M | $type == ICCM ) then
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
  if ( ! -d $executable_h )             mkdir -p $executable_h
  if ( ! -d $executable_h:h/lib_land_lad )    mkdir -p $executable_h:h/lib_land_lad

# compile libs
set makeFile      = Make_lib_land_lad
cd $executable_h:h/lib_land_lad

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_land_lad.a -c "$cppDefs"  -o "-I$executable_h:h/lib_FMS" $pathnames_land_lad $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile

if( $status ) then
    echo "Make failed to create  lib_land_lad.a"
    exit
endif

endif



if( $type == CM2M ) then
#   Build the CM2M executable   
    set srcList = ( coupler )
    set makeFile      = Make_CM2M
    if ( ! -d $executable_h ) mkdir -p $executable_h
    cd $executable_h
    $mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p $executable:t -o "-I$executable_h:h/lib_FMS -I$executable_h:h/lib_ocean -I$executable_h:h/lib_atmos_fv -I$executable_h:h/lib_atmos_phys -I$executable_h:h/lib_land_lad" -l "$executable_h:h/lib_ocean/lib_ocean.a $executable_h:h/lib_atmos_fv/lib_atmos_fv.a $executable_h:h/lib_atmos_phys/lib_atmos_phys.a $executable_h:h/lib_land_lad/lib_land_lad.a $executable_h:h/lib_FMS/lib_FMS.a"  $srcList

    make -f $makeFile
    if( $status ) then
	echo "Make failed to create the CM2M executable  "
	exit
    endif    
    exit
endif



if( $type == ESM2M ) then
# Build the land_lad2 library

##LAND Component
# The list of source files that should be compiled for this component.
set pathnames_land_lad2  = $code_dir/path_names_land_lad2        # path to file containing list of source paths

cat > $pathnames_land_lad2 <<EOF_land_lad2
land_lad2/canopy_air/cana_tile.F90  
land_lad2/canopy_air/canopy_air.F90  
land_lad2/land_constants.F90  
land_lad2/land_data.F90     
land_lad2/land_model.F90    
land_lad2/land_tile.F90     
land_lad2/glacier/glac_tile.F90   
land_lad2/glacier/glacier.F90  
land_lad2/lake/lake.F90       
land_lad2/lake/lake_tile.F90   
land_lad2/river/river.F90       
land_lad2/river/river_physics.F90  
land_lad2/river/river_type.F90      
land_lad2/shared/nf_utils/getput.inc     
land_lad2/shared/nf_utils/getput_compressed.inc 
land_lad2/shared/nf_utils/nf_utils.F90
land_lad2/shared/nf_utils/nfc.F90 
land_lad2/shared/nf_utils/nfu.F90
land_lad2/shared/debug.inc  
land_lad2/shared/land_debug.F90  
land_lad2/shared/land_io.F90
land_lad2/shared/land_numerics.F90
land_lad2/shared/land_tile_diag.F90
land_lad2/shared/land_tile_diag_buff.F90
land_lad2/shared/land_tile_diag_sel.F90
land_lad2/shared/land_tile_io.F90
land_lad2/shared/land_utils.F90
land_lad2/shared/sphum.F90
land_lad2/shared/table_printer.F90
land_lad2/snow/snow.F90
land_lad2/snow/snow_tile.F90
land_lad2/soil/soil.F90
land_lad2/soil/soil_tile.F90
land_lad2/soil/uptake.F90
land_lad2/topo_rough/topo_rough.F90
land_lad2/transitions/transitions.F90
land_lad2/vegetation/read_remap_cohort_data.inc
land_lad2/vegetation/vegetation.F90
land_lad2/vegetation/vegn_cohort.F90
land_lad2/vegetation/vegn_cohort_io.F90
land_lad2/vegetation/vegn_cohort_io.inc
land_lad2/vegetation/vegn_data.F90
land_lad2/vegetation/vegn_disturbance.F90
land_lad2/vegetation/vegn_dynamics.F90
land_lad2/vegetation/vegn_harvesting.F90
land_lad2/vegetation/vegn_photosynthesis.F90
land_lad2/vegetation/vegn_radiation.F90
land_lad2/vegetation/vegn_static_override.F90
land_lad2/vegetation/vegn_tile.F90

EOF_land_lad2

set srcList = ( )

# setup directory structure
  if ( ! -d $executable_h )             mkdir -p $executable_h
  if ( ! -d $executable_h:h/lib_land_lad2 )    mkdir -p $executable_h:h/lib_land_lad2

# compile libs
set makeFile      = Make_lib_land_lad2
cd $executable_h:h/lib_land_lad2

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_land_lad2.a -c "-DUSE_LOG_DIAG_FIELD_INFO $cppDefs"  -o "-I$executable_h:h/lib_FMS" $srcList $pathnames_land_lad2 $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_land_lad2.a"
    exit
endif

endif

if( $type == ESM2M ) then
#   Build the ESM2M executable   
    set srcList = ( coupler )
    set makeFile      = Make_ESM2M
    if ( ! -d $executable_h ) mkdir -p $executable_h
    cd $executable_h
    $mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p $executable:t -o "-I$executable_h:h/lib_FMS -I$executable_h:h/lib_ocean -I$executable_h:h/lib_atmos_fv -I$executable_h:h/lib_atmos_phys -I$executable_h:h/lib_land_lad2" -l "$executable_h:h/lib_ocean/lib_ocean.a $executable_h:h/lib_atmos_fv/lib_atmos_fv.a $executable_h:h/lib_atmos_phys/lib_atmos_phys.a $executable_h:h/lib_land_lad2/lib_land_lad2.a $executable_h:h/lib_FMS/lib_FMS.a"  $srcList

    make -f $makeFile
    if( $status ) then
	echo "Make failed to create the ESM2M executable  "
	exit
    endif    
    exit
endif


if( $type == ICCM ) then
# Build the atmos_BG library

set pathnames_atmos_bg  = $code_dir/path_names_atmos_bg        # path to file containing list of source paths

cat > $pathnames_atmos_bg <<EOF_atmos_bg
atmos_coupled/atmos_model.F90
atmos_bgrid/driver/coupled/atmosphere.F90
atmos_bgrid/driver/coupled/bgrid_physics.F90
atmos_bgrid/model/bgrid_advection.F90
atmos_bgrid/model/bgrid_conserve_energy.F90
atmos_bgrid/model/bgrid_core.F90
atmos_bgrid/model/bgrid_core_driver.F90
atmos_bgrid/model/bgrid_horiz_adjust.F90
atmos_bgrid/model/bgrid_horiz_diff.F90
atmos_bgrid/model/bgrid_sponge.F90
atmos_bgrid/model/bgrid_vert_adjust.F90
atmos_bgrid/tools/bgrid_change_grid.F90
atmos_bgrid/tools/bgrid_cold_start.F90
atmos_bgrid/tools/bgrid_diagnostics.F90
atmos_bgrid/tools/bgrid_halo.F90
atmos_bgrid/tools/bgrid_horiz.F90
atmos_bgrid/tools/bgrid_integrals.F90
atmos_bgrid/tools/bgrid_masks.F90
atmos_bgrid/tools/bgrid_polar_filter.F90
atmos_bgrid/tools/bgrid_prog_var.F90
atmos_bgrid/tools/bgrid_vert.F90


EOF_atmos_bg

# setup directory structure
  if ( ! -d $executable_h )             mkdir -p $executable_h
  if ( ! -d $executable_h:h/lib_atmos_bg )    mkdir -p $executable_h:h/lib_atmos_bg

# compile libs
set makeFile      = Make_lib_atmos_bg
cd $executable_h:h/lib_atmos_bg

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_atmos_bg.a -c "$cppDefs" -o "-I$executable_h:h/lib_FMS -I$executable_h:h/lib_atmos_phys" $pathnames_atmos_bg $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_atmos_bg.a"
    exit
endif

endif




if( $type == ICCM ) then
#   Build the ICCM executable   
    set srcList = ( coupler )
    set makeFile      = Make_ICCM
    if ( ! -d $executable_h ) mkdir -p $executable_h
    cd $executable_h
    $mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p $executable:t -o "-I$executable_h:h/lib_FMS -I$executable_h:h/lib_ocean -I$executable_h:h/lib_atmos_bg -I$executable_h:h/lib_atmos_phys -I$executable_h:h/lib_land_lad" -l "$executable_h:h/lib_ocean/lib_ocean.a $executable_h:h/lib_atmos_bg/lib_atmos_bg.a $executable_h:h/lib_atmos_phys/lib_atmos_phys.a $executable_h:h/lib_land_lad/lib_land_lad.a $executable_h:h/lib_FMS/lib_FMS.a"  $srcList

    make -f $makeFile
    if( $status ) then
	echo "Make failed to create the ICCM executable  "
	exit
    endif    
    exit
endif
