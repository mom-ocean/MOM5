# Build the am3 library

set srcList = ( atmos_param/diag_integral atmos_param/monin_obukhov atmos_shared/vert_advection atmos_shared/tracer_driver/tropchem atmos_shared/tracer_driver/stratchem atmos_shared/tracer_driver atmos_shared/tracer_driver/aer_ccn_act atmos_shared/interpolator atmos_shared/atmos_nudge atmos_param/vert_turb_driver atmos_param/vert_diff_driver atmos_param/vert_diff atmos_param/topo_drag atmos_param/strat_cloud atmos_param/stable_bl_turb atmos_param/shallow_physics atmos_param/shallow_cu atmos_param/shallow_conv atmos_param/sea_esf_rad atmos_param/rh_clouds atmos_param/ras atmos_param/radiation_driver atmos_param/physics_driver atmos_param/my25_turb atmos_param/moist_processes atmos_param/moist_conv atmos_param/mg_drag atmos_param/lscale_cond atmos_param/lin_cloud_microphys atmos_param/grey_radiation atmos_param/fsrad atmos_param/entrain atmos_param/edt atmos_param/dry_adj atmos_param/donner_deep atmos_param/diffusivity atmos_param/diag_cloud_rad atmos_param/diag_cloud atmos_param/damping_driver atmos_param/cu_mo_trans atmos_param/cosp/quickbeam atmos_param/cosp/llnl atmos_param/cosp/icarus-scops-3.7 atmos_param/cosp/actsim atmos_param/cosp/MODIS_simulator atmos_param/cosp/MISR_simulator atmos_param/cosp atmos_param/clouds atmos_param/cloud_zonal atmos_param/cloud_rad atmos_param/cloud_obs atmos_param/cloud_generator atmos_param/cg_drag atmos_param/betts_miller atmos_param/clubb )

set lib_name = "lib_atmos_phys"

mkdir -p $executable:h:h/$lib_name
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS" $srcList $lib_include_dirs
make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
