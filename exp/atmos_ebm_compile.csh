# Build the am3 library

set pathnames_atmos_ebm  = $code_dir/path_names_atmos_ebm        # path to file containing list of source paths

cat > $pathnames_atmos_ebm <<EOF_atmos_ebm
atmos_coupled/atmos_model.F90
atmos_ebm/atmosphere.F90
atmos_ebm/ebm_diagnostics.F90
atmos_param/edt/null/edt.F90
atmos_param/diffusivity/null/diffusivity.F90
atmos_param/vert_turb_driver/vert_turb_driver.F90
atmos_param/topo_drag/null/topo_drag.F90
atmos_param/monin_obukhov/monin_obukhov.F90
atmos_param/monin_obukhov/monin_obukhov_kernel.F90
atmos_param/monin_obukhov/monin_obukhov_interfaces.h
atmos_param/diag_cloud/null/diag_cloud.F90
atmos_param/cloud_rad/cloud_rad.F90
atmos_param/vert_diff/vert_diff.F90
atmos_param/vert_diff_driver/vert_diff_driver.F90
atmos_param/cu_mo_trans/cu_mo_trans.F90
atmos_param/cloud_generator/null/cloud_generator.F90
atmos_param/cloud_generator/null/betaDistribution.F90
atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
atmos_param/physics_driver/physics_driver.F90
atmos_param/stable_bl_turb/stable_bl_turb.F90
atmos_param/rh_clouds/null/rh_clouds.F90
atmos_param/cloud_zonal/null/cloud_zonal.F90
atmos_param/shallow_cu/null/uw_conv.F90
atmos_param/cosp/null/cosp_driver.F90
atmos_param/fsrad/null/mcm_sw_driver.F90
atmos_param/fsrad/null/shortwave.F90
atmos_param/fsrad/null/hconst.F90
atmos_param/fsrad/null/longwave.F90
atmos_param/fsrad/null/fs_profile.F90
atmos_param/fsrad/null/co2_data.F90
atmos_param/fsrad/null/fsrad.F90
atmos_param/fsrad/null/co2int.F90
atmos_param/fsrad/null/mcm_lw.F90
atmos_param/fsrad/null/mcm_swtbls.F90
atmos_param/fsrad/null/mcm_swnew.F90
atmos_param/fsrad/null/rad_diag.F90
atmos_param/fsrad/null/rdparm.F90
atmos_param/moist_processes/moist_processes_utils.F90
atmos_param/moist_processes/detr_ice_num.F90
atmos_param/moist_processes/moistproc_kernels.F90
atmos_param/moist_processes/moist_processes.F90
atmos_param/cg_drag/null/cg_drag.F90
atmos_param/clouds/null/clouds.F90
atmos_param/entrain/entrain.F90
atmos_param/betts_miller/betts_miller.F90
atmos_param/betts_miller/bm_omp.F90
atmos_param/betts_miller/bm_massflux.F90
atmos_param/donner_deep/null/donner_deep.F90
atmos_param/my25_turb/null/my25_turb.F90
atmos_param/ras/ras.F90
atmos_param/damping_driver/damping_driver.F90
atmos_param/grey_radiation/grey_radiation.F90
atmos_param/lscale_cond/null/lscale_cond.F90
atmos_param/mg_drag/mg_drag.F90
atmos_param/sea_esf_rad/gas_tf.F90
atmos_param/sea_esf_rad/optical_path.F90
atmos_param/sea_esf_rad/radiative_gases.F90
atmos_param/sea_esf_rad/longwave_driver.F90
atmos_param/sea_esf_rad/cloud_spec.F90
atmos_param/sea_esf_rad/ozone.F90
atmos_param/sea_esf_rad/cloudrad_diagnostics.F90
atmos_param/sea_esf_rad/shortwave_driver.F90
atmos_param/sea_esf_rad/longwave_clouds.F90
atmos_param/sea_esf_rad/esfsw_parameters.F90
atmos_param/sea_esf_rad/sealw99.F90
atmos_param/sea_esf_rad/esfsw_driver.F90
atmos_param/sea_esf_rad/strat_clouds_W.F90
atmos_param/sea_esf_rad/null/uw_clouds_W.F90
atmos_param/sea_esf_rad/null/standalone_clouds.F90
atmos_param/sea_esf_rad/null/donner_deep_clouds_W.F90
atmos_param/sea_esf_rad/null/bulkphys_rad.F90
atmos_param/sea_esf_rad/null/lhsw_driver.F90
atmos_param/sea_esf_rad/null/diag_clouds_W.F90
atmos_param/sea_esf_rad/null/original_fms_rad.F90
atmos_param/sea_esf_rad/null/specified_clouds_W.F90
atmos_param/sea_esf_rad/null/mgrp_prscr_clds.F90
atmos_param/sea_esf_rad/null/rh_based_clouds.F90
atmos_param/sea_esf_rad/zetac_clouds_W.F90
atmos_param/sea_esf_rad/lw_gases_stdtf.F90
atmos_param/sea_esf_rad/radiation_diag.F90
atmos_param/sea_esf_rad/isccp_clouds.F90
atmos_param/sea_esf_rad/sea_esf_rad.F90
atmos_param/sea_esf_rad/longwave_params.F90
atmos_param/sea_esf_rad/microphys_rad.F90
atmos_param/sea_esf_rad/rad_utilities.F90
atmos_param/sea_esf_rad/microphys_cloud.F90
atmos_param/sea_esf_rad/longwave_tables.F90
atmos_param/sea_esf_rad/aerosol.F90
atmos_param/sea_esf_rad/longwave_fluxes.F90
atmos_param/sea_esf_rad/aerosolrad_package.F90
atmos_param/sea_esf_rad/cloudrad_package.F90
atmos_param/sea_esf_rad/rad_output_file.F90
atmos_param/strat_cloud/strat_cloud.F90
atmos_param/strat_cloud/rotstayn_klein_mp.F90
atmos_param/strat_cloud/strat_nml.h
atmos_param/strat_cloud/polysvp.F90
atmos_param/strat_cloud/simple_pdf.F90
atmos_param/strat_cloud/strat_cloud_utilities.F90
atmos_param/strat_cloud/check_nan.F90
atmos_param/strat_cloud/gamma_mg.F90
atmos_param/strat_cloud/mg_const.F90
atmos_param/strat_cloud/aerosol_cloud.F90
atmos_param/strat_cloud/strat_netcdf.F90
atmos_param/strat_cloud/nc_cond.F90
atmos_param/strat_cloud/microphysics.F90
atmos_param/strat_cloud/morrison_gettelman_microp.F90
atmos_param/strat_cloud/strat_cloud_legacy.F90
atmos_param/cloud_obs/null/cloud_obs.F90
atmos_param/diag_integral/diag_integral.F90
atmos_param/shallow_conv/null/shallow_conv.F90
atmos_param/moist_conv/null/moist_conv.F90
atmos_param/radiation_driver/radiation_driver.F90
atmos_param/dry_adj/null/dry_adj.F90
atmos_param/diag_cloud_rad/null/diag_cloud_rad.F90
atmos_shared/interpolator/interpolator.F90
atmos_shared/vert_advection/vert_advection.F90
atmos_shared/atmos_nudge/atmos_nudge.F90
atmos_shared/tracer_driver/atmos_sulfur_hex.F90
atmos_shared/tracer_driver/atmos_tracer_utilities.F90
atmos_shared/tracer_driver/atmos_soa.F90
atmos_shared/tracer_driver/atmos_sea_salt.F90
atmos_shared/tracer_driver/atmos_dust.F90
atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90
atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90
atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
atmos_shared/tracer_driver/atmos_convection_tracer.F90
atmos_shared/tracer_driver/atmos_radon.F90
atmos_shared/tracer_driver/atmos_age_tracer.F90
atmos_shared/tracer_driver/atmos_co2.F90
atmos_shared/tracer_driver/atmos_tracer_driver.F90
atmos_shared/tracer_driver/atmos_carbon_aerosol.F90
atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90
atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90
atmos_shared/tracer_driver/tropchem/moz.subs.F90
atmos_shared/tracer_driver/tropchem/m_tracname.F90
atmos_shared/tracer_driver/tropchem/mo_setinv.F90
atmos_shared/tracer_driver/tropchem/moz.mods.F90
atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
atmos_shared/tracer_driver/tropchem/mo_chemdr.F90
atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90
atmos_shared/tracer_driver/tropchem/mo_jpl.F90
atmos_shared/tracer_driver/tropchem/moz.mat.F90
atmos_shared/tracer_driver/tropchem/mo_chemini.F90
atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90
atmos_shared/tracer_driver/tropchem/mo_fastjx.F90
atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90
atmos_shared/tracer_driver/tropchem/mo_hook.F90
atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90
atmos_shared/tracer_driver/tropchem/mo_setsox.F90
atmos_shared/tracer_driver/tropchem/mo_photo.F90
atmos_shared/tracer_driver/tropchem/mo_fphoto.F90
atmos_shared/tracer_driver/tropchem/tropchem_driver.F90
atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90
atmos_shared/tracer_driver/atmos_ch3i.F90
atmos_shared/tracer_driver/atmos_sulfate.F90
atmos_spectral/tools/gauss_and_legendre.F90
atmos_spectral/tools/grid_fourier.F90
atmos_spectral/tools/spec_mpp.F90
atmos_spectral/tools/spherical.F90
atmos_spectral/tools/spherical_fourier.F90
atmos_spectral/tools/transforms.F90

EOF_atmos_ebm

set lib_name = "lib_atmos_ebm"
# setup directory structure
mkdir -p $executable:h:h/$lib_name

# compile libs
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS" $pathnames_atmos_ebm $lib_include_dirs

make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
