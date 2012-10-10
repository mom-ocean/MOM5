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
  if ( ! -d $executable:h )             mkdir -p $executable:h
  if ( ! -d $executable:h:h/lib_atmos_phys_fv )    mkdir -p $executable:h:h/lib_atmos_fv

# compile libs
set makeFile      = Make_lib_atmos_fv
cd $executable:h:h/lib_atmos_fv

$mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p lib_atmos_fv.a -c "$cppDefs" -o "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_atmos_phys" $pathnames_atmos_fv $root/include $code_dir/shared/include $code_dir/shared/mpp/include

make -f $makeFile 

if( $status ) then
    echo "Make failed to create  lib_atmos_fv.a"
    exit 1
endif
