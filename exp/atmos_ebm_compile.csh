# Build the EBM library

set srcList = ( atmos_spectral/tools atmos_shared/tracer_driver/stratchem atmos_shared/tracer_driver/tropchem atmos_shared/tracer_driver/aer_ccn_act atmos_shared/tracer_driver atmos_shared/atmos_nudge atmos_shared/vert_advection atmos_shared/interpolator atmos_param/diag_cloud_rad/null atmos_param/dry_adj/null atmos_param/radiation_driver atmos_param/moist_conv/null atmos_param/shallow_conv/null atmos_param/diag_integral atmos_param/cloud_obs/null atmos_param/strat_cloud atmos_param/sea_esf_rad/null atmos_param/sea_esf_rad atmos_param/mg_drag atmos_param/lscale_cond/null atmos_param/grey_radiation atmos_param/damping_driver atmos_param/ras atmos_param/my25_turb/null atmos_param/donner_deep/null atmos_param/betts_miller atmos_param/entrain atmos_param/clouds/null atmos_param/cg_drag/null atmos_param/moist_processes atmos_param/fsrad/null atmos_param/cosp/null atmos_param/shallow_cu/null atmos_param/cloud_zonal/null atmos_param/rh_clouds/null atmos_param/stable_bl_turb atmos_param/physics_driver atmos_param/lin_cloud_microphys atmos_param/cloud_generator/null atmos_param/cu_mo_trans atmos_param/vert_diff_driver atmos_param/vert_diff atmos_param/cloud_rad atmos_param/diag_cloud/null atmos_param/monin_obukhov atmos_param/topo_drag atmos_param/vert_turb_driver atmos_param/diffusivity/null atmos_param/edt/null atmos_ebm atmos_coupled )

set lib_name = "lib_atmos_ebm"
# setup directory structure
mkdir -p $executable:h:h/$lib_name

# compile libs
cd $executable:h:h/$lib_name
$mkmf_lib -p $lib_name.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS" $lib_include_dirs

make

if( $status ) then
    echo "Make failed to create $lib_name.a"
    exit 1
endif
