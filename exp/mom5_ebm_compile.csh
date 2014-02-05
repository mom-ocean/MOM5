#!/bin/csh -f
# Minimal compile script for mom4p1 solo experiments

set echo
set platform      = ncrc2.intel     # A unique identifier for your platform
                                   # This corresponds to the mkmf templates in $root/bin dir.
#
# User does not need to change anything below!
#
set type          = EBM                               # Name of the experiment
set root          = $cwd:h                            # The directory you created when you checkout
set code_dir      = $root/src                         # source code directory
set executable    = $root/exec/$platform/$type/fms_$type.x      # executable created after compilation
set pathnames     = $code_dir/path_names_$type        # path to file containing list of source paths
set mppnccombine  = $root/bin/mppnccombine.$platform  # path to executable mppnccombine
set mkmfTemplate  = $root/bin/mkmf.template.$platform # path to template for your platform
set mkmf          = $root/bin/mkmf                    # path to executable mkmf
set cppDefs  = ( "-Duse_netCDF  -Duse_netCDF3 -Duse_libMPI -DLAND_BND_TRACERS -DOVERLOAD_C8 -DOVERLOAD_C4 -DOVERLOAD_R4" )

#
# Users must ensure the correct environment file exists for their platform.
#
source $root/bin/environs.$platform  # environment variables and loadable modules

# setup directory structure
  if ( ! -d $executable:h )    mkdir -p $executable:h

#
# compile mppnccombine.c, needed only if $npes > 1
#NOTE: On some platforms you may need to specify the location for netcdf.h and libnetcdf.a
#      by modifying the following -I and -L
  if (  ! -f $mppnccombine ) then
    cc -O -o $mppnccombine -I/usr/local/include -L/usr/local/lib $code_dir/postprocessing/mppnccombine/mppnccombine.c -lnetcdf
  endif



# The list of source files that should be compiled for this experiment.
cat > $pathnames <<EOF
atmos_coupled/atmos_model.F90
atmos_ebm/atmosphere.F90
atmos_ebm/ebm_diagnostics.F90
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
atmos_param/cosp/null/cosp_driver.F90                                
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
atmos_param/diag_integral/diag_integral.F90
atmos_param/monin_obukhov/monin_obukhov.F90
atmos_param/monin_obukhov/monin_obukhov_interfaces.h
atmos_param/monin_obukhov/monin_obukhov_kernel.F90
atmos_spectral/tools/gauss_and_legendre.F90
atmos_spectral/tools/grid_fourier.F90
atmos_spectral/tools/spec_mpp.F90
atmos_spectral/tools/spherical.F90
atmos_spectral/tools/spherical_fourier.F90
atmos_spectral/tools/transforms.F90
coupler/coupler_main.F90
coupler/flux_exchange.F90
coupler/surface_flux.F90
ice_param/ice_albedo.F90
ice_param/ocean_albedo.F90
ice_param/ocean_rough.F90
ice_sis/ice_bergs.F90
ice_sis/ice_dyn.F90
ice_sis/ice_grid.F90
ice_sis/ice_model.F90
ice_sis/ice_spec.F90
ice_sis/ice_thm.F90
ice_sis/ice_type.F90
ice_sis/mask.F90
ice_sis/rot.F90
land_lad/land_model.F90
land_lad/land_types.F90
land_lad/numerics.F90
land_lad/soil/land_properties.F90
land_lad/soil/rivers.F90
land_lad/soil/soil.F90
land_lad/vegetation/vegetation.F90
land_param/climap_albedo.F90
ocean_shared/generic_tracers/FMS_coupler_util.F90
ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90
ocean_shared/generic_tracers/generic_CFC.F90
ocean_shared/generic_tracers/generic_COBALT.F90
ocean_shared/generic_tracers/generic_TOPAZ.F90
ocean_shared/generic_tracers/generic_BLING.F90
ocean_shared/generic_tracers/generic_miniBLING.F90
ocean_shared/generic_tracers/generic_ERGOM.F90
ocean_shared/generic_tracers/generic_tracer.F90
ocean_shared/generic_tracers/generic_tracer_utils.F90
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
shared/drifters/drifters_comm.F90
shared/drifters/drifters_compute_k.h
shared/drifters/drifters_core.F90
shared/drifters/drifters.F90
shared/drifters/drifters_input.F90
shared/drifters/drifters_io.F90
shared/drifters/drifters_push.h
shared/drifters/drifters_set_field.h
shared/drifters/fms_switches.h
shared/drifters/quicksort.F90
shared/exchange/stock_constants.F90
shared/exchange/xgrid.F90
shared/fft/fft99.F90
shared/fft/fft.F90
shared/field_manager/field_manager.F90
shared/field_manager/fm_util.F90
shared/field_manager/parse.inc
shared/fms/fms.F90
shared/fms/fms_io.F90
shared/fms/read_data_2d.inc
shared/fms/read_data_3d.inc
shared/fms/read_data_4d.inc
shared/fms/test_fms_io.F90
shared/fms/write_data.inc
shared/horiz_interp/horiz_interp_bicubic.F90
shared/horiz_interp/horiz_interp_bilinear.F90
shared/horiz_interp/horiz_interp_conserve.F90
shared/horiz_interp/horiz_interp.F90
shared/horiz_interp/horiz_interp_spherical.F90
shared/horiz_interp/horiz_interp_type.F90
shared/include/fms_platform.h
shared/memutils/memuse.c
shared/memutils/memutils.F90
shared/mosaic/constant.h
shared/mosaic/create_xgrid.c
shared/mosaic/create_xgrid.h
shared/mosaic/gradient_c2l.c
shared/mosaic/gradient_c2l.h
shared/mosaic/gradient.F90
shared/mosaic/grid.F90
shared/mosaic/interp.c
shared/mosaic/interp.h
shared/mosaic/mosaic.F90
shared/mosaic/mosaic_util.c
shared/mosaic/mosaic_util.h
shared/mosaic/read_mosaic.c
shared/mosaic/read_mosaic.h
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
shared/mpp/include/mpp_domains_comm.inc
shared/mpp/include/mpp_domains_define.inc
shared/mpp/include/mpp_domains_misc.inc
shared/mpp/include/mpp_domains_reduce.inc
shared/mpp/include/mpp_domains_util.inc
shared/mpp/include/mpp_do_redistribute.h
shared/mpp/include/mpp_do_update_ad.h
shared/mpp/include/mpp_do_update.h
shared/mpp/include/mpp_do_updateV_ad.h
shared/mpp/include/mpp_do_updateV.h
shared/mpp/include/mpp_error_a_a.h
shared/mpp/include/mpp_error_a_s.h
shared/mpp/include/mpp_error_s_a.h
shared/mpp/include/mpp_error_s_s.h
shared/mpp/include/mpp_get_boundary.h
shared/mpp/include/mpp_global_field.h
shared/mpp/include/mpp_global_reduce.h
shared/mpp/include/mpp_global_sum_ad.h
shared/mpp/include/mpp_global_sum.h
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
shared/mpp/include/mpp_update_domains2D_ad.h
shared/mpp/include/mpp_update_domains2D.h
shared/mpp/include/mpp_util.inc
shared/mpp/include/mpp_util_mpi.inc
shared/mpp/include/mpp_util_nocomm.inc
shared/mpp/include/mpp_util_sma.inc
shared/mpp/include/mpp_write_2Ddecomp.h
shared/mpp/include/mpp_write.h
shared/mpp/include/system_clock.h
shared/mpp/mpp_data.F90
shared/mpp/mpp_domains.F90
shared/mpp/mpp.F90
shared/mpp/mpp_io.F90
shared/mpp/mpp_memutils.F90
shared/mpp/mpp_parameter.F90
shared/mpp/mpp_pset.F90
shared/mpp/mpp_utilities.F90
shared/mpp/nsclock.c
shared/mpp/test_mpp_domains.F90
shared/mpp/test_mpp.F90
shared/mpp/test_mpp_io.F90
shared/mpp/test_mpp_pset.F90
shared/mpp/threadloc.c
shared/platform/platform.F90
shared/random_numbers/MersenneTwister.F90
shared/random_numbers/random_numbers.F90
shared/sat_vapor_pres/sat_vapor_pres.F90
shared/sat_vapor_pres/sat_vapor_pres_k.F90
shared/time_interp/time_interp_external.F90
shared/time_interp/time_interp.F90
shared/time_manager/get_cal_time.F90
shared/time_manager/time_manager.F90
shared/topography/gaussian_topog.F90
shared/topography/topography.F90
shared/tracer_manager/tracer_manager.F90
shared/tridiagonal/tridiagonal.F90

EOF

set srcList = (  mom5/ocean_bgc mom5/ocean_core mom5/ocean_diag mom5/ocean_blobs mom5/ocean_param/neutral mom5/ocean_param/sources mom5/ocean_param/lateral mom5/ocean_param/vertical mom5/ocean_param/gotm-4.0/include mom5/ocean_param/gotm-4.0/turbulence mom5/ocean_param/gotm-4.0/util mom5/ocean_tracers mom5/ocean_wave        )


# compile the model code and create executable
  set makeFile      = Make_$type
  cd $executable:h
  $mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p $executable:t -c "$cppDefs" $srcList $pathnames $root/include $code_dir/shared/include $code_dir/shared/mpp/include 
#/usr/local/include

  make -f $makeFile

