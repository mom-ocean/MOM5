esma_set_this ()

set (SRCS
  src/mom5/ocean_tracers/ocean_residency_range.F90 
  src/mom5/ocean_param/neutral/ocean_nphysics_util_new.F90
  src/mom5/ocean_param/lateral/ocean_sigma_transport.F90 
  src/mom5/ocean_blobs/ocean_blob_dynamic_free.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90
  src/mom5/ocean_param/neutral/ocean_nphysics.F90 
  src/mom5/ocean_param/vertical/ocean_vert_gotm.F90
  src/mom5/ocean_tracers/ocean_residency_integrand.F90 
  src/mom5/ocean_param/sources/ocean_increment_tracer.F90 
  src/mom5/ocean_core/ocean_model.F90
  src/ocean_shared/generic_tracers/generic_GEN.F90 
  src/ice_sis/ice_grid.F90 
  src/mom5/ocean_bgc/ocmip2_biotic.F90
  src/mom5/ocean_blobs/ocean_blob_diag.F90 
  src/ocean_shared/generic_tracers/FMS_coupler_util.F90 
  src/mom5/ocean_diag/ocean_velocity_diag.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 
  src/ice_param/ocean_albedo.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90
  src/mom5/ocean_param/sources/ocean_shortwave_gfdl.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 
  src/mom5/ocean_param/lateral/ocean_lap_tracer.F90 
  src/mom5/ocean_param/gotm-4.0/util/time.F90
  src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90
  src/mom5/ocean_param/gotm-4.0/util/util.F90 
  src/mom5/ocean_param/gotm-4.0/util/lagrange.F90 
  src/mom5/ocean_param/lateral/ocean_lapcst_friction.F90
  src/mom5/ocean_param/gotm-4.0/util/eqstate.F90 
  src/mom5/ocean_param/sources/ocean_xlandmix.F90 
  src/mom5/ocean_tracers/ocean_residency.F90
  src/ice_sis/ice_bergs.F90 
  src/mom5/ocean_param/vertical/ocean_form_drag.F90 
  src/mom5/ocean_param/vertical/ocean_vert_chen.F90
  src/mom5/ocean_param/neutral/ocean_nphysicsC.F90 
  src/mom5/ocean_tracers/ocean_tempsalt.F90 src/ice_sis/ice_model.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90 
  src/mom5/ocean_tracers/ocean_age_tracer.F90 
  src/ocean_shared/generic_tracers/generic_ERGOM.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90 
  src/mom5/ocean_param/sources/ocean_overflow_OFP.F90
  src/mom5/ocean_param/neutral/ocean_nphysics_flux.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90 
  src/mom5/ocean_param/sources/ocean_xlandinsert.F90
  src/mom5/ocean_param/lateral/ocean_bih_tracer.F90 
  src/mom5/ocean_core/ocean_thickness.F90 
  src/mom5/ocean_core/oda_driver.F90
  src/mom5/ocean_param/vertical/ocean_vert_kpp_test.F90 
  src/mom5/ocean_param/sources/ocean_sponges_velocity.F90 
  src/mom5/ocean_core/ocean_velocity_advect.F90
  src/mom5/ocean_diag/ocean_tracer_diag.F90 
  src/mom5/ocean_tracers/ocean_residency_meta.F90 
  src/mom5/ocean_param/vertical/ocean_vert_const.F90
  src/mom5/ocean_core/ocean_density.F90 
  src/mom5/ocean_param/neutral/ocean_nphysics_diff.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90
  src/mom5/ocean_core/ocean_workspace.F90
  src/mom5/ocean_core/ocean_operators.F90 
  src/ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 
  src/mom5/ocean_core/ocean_coriolis.F90 
  src/mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90
  src/mom5/ocean_diag/ocean_drifters.F90 
  src/mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 
  src/ice_sis/ice_dyn.F90
  src/mom5/ocean_param/sources/ocean_momentum_source.F90 
  src/mom5/ocean_diag/ocean_adv_vel_diag.F90 
  src/ocean_shared/generic_tracers/generic_tracer_utils.F90
  src/mom5/ocean_param/neutral/ocean_nphysics_tensor.F90 
  src/ocean_shared/generic_tracers/generic_CFC.F90 
  src/mom5/ocean_param/neutral/ocean_nphysicsB.F90
  src/mom5/ocean_bgc/ocean_pert_co2.F90 
  src/mom5/ocean_core/ocean_pressure.F90 
  src/mom5/ocean_param/lateral/ocean_lapgen_friction.F90
  src/ocean_shared/generic_tracers/generic_TOPAZ.F90 
  src/mom5/ocean_param/sources/ocean_overflow.F90 
  src/ice_sis/rot.F90
  src/mom5/ocean_param/vertical/ocean_vert_util.F90 
  src/mom5/ocean_param/neutral/ocean_nphysics_skew.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90 
  src/mom5/ocean_param/sources/ocean_sponges_tracer.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/production.F90 
  src/mom5/ocean_param/sources/ocean_increment_eta.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90 
  src/mom5/ocean_param/gotm-4.0/util/gridinterpol.F90 
  src/mom5/ocean_param/gotm-4.0/util/diff_center.F90
  src/mom5/ocean_param/gotm-4.0/util/adv_center.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 
  src/mom5/ocean_tracers/ocean_tpm_util.F90
  src/ice_param/ocean_rough.F90 
  src/mom5/ocean_wave/ocean_wave_type.F90 
  src/mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90
  src/mom5/ocean_tracers/ocean_passive.F90 src/ice_sis/mask.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90
  src/mom5/ocean_core/ocean_obc.F90 src/mom5/ocean_wave/ocean_wave.F90 
  src/mom5/ocean_param/neutral/ocean_nphysics_new.F90
  src/mom5/ocean_bgc/ocean_bgc_restore.F90 
  src/mom5/ocean_param/sources/ocean_increment_velocity.F90 
  src/mom5/ocean_param/vertical/ocean_vert_mix.F90
  src/mom5/ocean_core/ocean_domains.F90 
  src/mom5/ocean_core/ocean_velocity.F90 
  src/mom5/ocean_param/sources/ocean_shortwave.F90
  src/mom5/ocean_core/ocean_grids.F90 
  src/mom5/ocean_core/ocean_topog.F90 
  src/mom5/ocean_bgc/ocmip2_abiotic.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 
  src/mom5/ocean_param/gotm-4.0/util/diff_face.F90
  src/mom5/ocean_param/sources/ocean_shortwave_csiro.F90 
  src/mom5/ocean_bgc/ocmip2_he.F90 
  src/mom5/ocean_bgc/ocean_generic_tracer.F90
  src/mom5/ocean_bgc/ocmip2_cfc.F90 
  src/mom5/ocean_param/lateral/ocean_lap_friction.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90
  src/ocean_shared/generic_tracers/generic_BLING.F90
  src/ocean_shared/generic_tracers/generic_miniBLING.F90 
  src/ocean_shared/generic_tracers/generic_COBALT.F90
  src/mom5/ocean_core/ocean_util.F90 
  src/mom5/ocean_core/ocean_parameters.F90
  src/mom5/ocean_param/sources/ocean_overexchange.F90 
  src/mom5/ocean_param/sources/ocean_shortwave_jerlov.F90 
  src/mom5/ocean_param/vertical/ocean_vert_pp.F90
  src/mom5/ocean_param/sources/ocean_sponges_eta.F90 
  src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90 
  src/ice_param/ice_albedo.F90
  src/mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90 
  src/mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 
  src/mom5/ocean_blobs/ocean_blob_util.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 
  src/mom5/ocean_param/neutral/ocean_nphysicsA.F90 
  src/mom5/ocean_tracers/ocean_tracer.F90
  src/ice_sis/ice_thm.F90 
  src/mom5/ocean_tracers/ocean_residency_ml.F90 
  src/ice_sis/ice_spec.F90
  src/mom5/ocean_core/ocean_obc_barotrop.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90 
  src/mom5/ocean_param/sources/ocean_riverspread.F90
  src/mom5/ocean_param/vertical/ocean_convect.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90 
  src/mom5/ocean_core/ocean_sbc.F90
  src/mom5/ocean_param/lateral/ocean_mixdownslope.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90 
  src/mom5/ocean_blobs/ocean_blob_static_free.F90
  src/mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 
  src/mom5/ocean_param/lateral/ocean_bihcst_friction.F90 
  src/mom5/ocean_diag/ocean_tracer_util.F90
  src/mom5/ocean_core/ocean_barotropic.F90 
  src/mom5/ocean_param/lateral/ocean_bih_friction.F90 
  src/mom5/ocean_tracers/ocean_frazil.F90
  src/mom5/ocean_param/lateral/ocean_bihgen_friction.F90 
  src/mom5/ocean_bgc/ocean_po4_pre.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90
  src/mom5/ocean_tracers/ocean_tpm.F90 
  src/mom5/ocean_param/sources/ocean_rivermix.F90 
  src/ice_sis/ice_type.F90
  src/mom5/ocean_blobs/ocean_blob.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90 
  src/mom5/ocean_diag/ocean_diagnostics.F90
  src/mom5/ocean_blobs/ocean_blob_static_bottom.F90 
  src/ocean_shared/generic_tracers/generic_tracer.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90
  src/mom5/ocean_param/lateral/ocean_submesoscale.F90 
  src/mom5/ocean_tracers/ocean_tracer_advect.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90 
  src/mom5/ocean_core/ocean_bbc.F90 
  src/mom5/ocean_bgc/ocmip2_co2calc.F90
  src/mom5/ocean_param/vertical/ocean_vert_tidal.F90 
  src/mom5/ocean_core/ocean_types.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90
  src/mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90 
  src/mom5/ocean_tracers/transport_matrix.F90 
  src/mom5/ocean_param/neutral/ocean_nphysics_util.F90
  src/mom5/ocean_bgc/ocean_ibgc.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/variances.F90 
  src/mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90
  src/mom5/ocean_core/ocean_advection_velocity.F90 
  )


esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES FMS_r8
  INCLUDES src/mom5/ocean_param/gotm-4.0/include src/mom5/ocean_core
)

target_compile_definitions (${this} PRIVATE MAPL_MODE EIGHT_BYTE SPMD TIMING use_libMPI use_netCDF USE_OCEAN_BGC)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel AND CMAKE_BUILD_TYPE MATCHES Release)
   set(CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_Release_Flags} ${common_Fortran_Flags} -stack_temps -safe_cray_ptr -assume byterecl -fp-model precise -fp-model source -ftz -align all -fno-alias -align dcommons ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")
endif ()


if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   target_compile_definitions (${this} PRIVATE NO_R16 NO_CRAY_POINTERS)
endif ()

# Specs for r8 version
string(REPLACE " " ";" tmp ${FREAL8})
foreach (flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach ()

install(
   DIRECTORY data/geos5
   DESTINATION etc/MOM5
   USE_SOURCE_PERMISSIONS
   )

