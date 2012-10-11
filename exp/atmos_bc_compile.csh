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
mkdir -p $executable:h:h/lib_atmos_bg

# compile libs
cd $executable:h:h/lib_atmos_bg
$mkmf_lib -p lib_atmos_bg.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_atmos_phys" $pathnames_atmos_bg $lib_include_dir

make

if( $status ) then
    echo "Make failed to create  lib_atmos_bg.a"
    exit 1
endif
