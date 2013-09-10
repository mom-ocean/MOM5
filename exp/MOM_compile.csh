#!/bin/csh -f
# Minimal compile script for fully coupled model CM2M experiments

set echo
set platform      = gfortran   # A unique identifier for your platform
                               # This corresponds to the mkmf templates in $root/bin dir.
set type          = MOM_solo   # Type of the experiment
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
    echo "             EBM      : ocean-seaice-land-atmosphere coupled model with energy balance atmosphere"
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

if ( $type == EBM ) then
    set cppDefs  = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI -DLAND_BND_TRACERS -DOVERLOAD_C8 -DOVERLOAD_C4 -DOVERLOAD_R4" )
else if( $type == ACCESS-OM ) then
    set cppDefs  = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI -DACCESS" )
else if( $type == ACCESS-CM ) then
    set cppDefs  = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI -DACCESS -DACCESS_CM" )
endif

#
# Users must ensure the correct environment file exists for their platform.
#
source $root/bin/environs.$platform  # environment variables and loadable modules

#
# compile mppnccombine.c, needed only if $npes > 1
  if ( ! -f $mppnccombine ) then
    cc -O -o $mppnccombine -I/usr/local/include -L/usr/local/lib $code_dir/postprocessing/mppnccombine/mppnccombine.c -lnetcdf
  endif

set mkmf_lib = "$mkmf -f -m Makefile -a $code_dir -t $mkmfTemplate"
set lib_include_dirs = "$root/include $code_dir/shared/include $code_dir/shared/mpp/include"

# Build FMS.
source ./FMS_compile.csh
set includes = "-I$executable:h:h/lib_FMS -I$executable:h:h/lib_ocean"

# Build the core ocean.
cd $root/exp
source ./ocean_compile.csh
if ( $status ) exit $status

if( $type != MOM_solo && $type != ACCESS-OM  && $type != ACCESS-CM) then
    cd $root/exp
    source ./ice_compile.csh
    if ( $status ) exit $status
endif
if( $type == MOM_SIS) then
    cd $root/exp
    source ./land_null_compile.csh
    if ( $status ) exit $status

    cd $root/exp
    source ./atmos_null_compile.csh
    if ( $status ) exit $status
endif
if( $type == EBM) then
    cd $root/exp
    source ./atmos_ebm_compile.csh
    if ( $status ) exit $status
endif
if( $type == CM2M | $type == ESM2M | $type == ICCM ) then
    cd $root/exp
    source ./atmos_phys_compile.csh
    if ( $status ) exit $status
endif
if( $type == CM2M | $type == ESM2M ) then
    cd $root/exp
    source ./atmos_fv_compile.csh
    if ( $status ) exit $status
endif
if( $type == CM2M | $type == ICCM | $type == EBM ) then
    cd $root/exp
    source ./land_lad_compile.csh
    if ( $status ) exit $status
endif
if( $type == ESM2M ) then
    cd $root/exp
    source ./land_lad2_compile.csh
    if ( $status ) exit $status
endif
if( $type == ICCM ) then
    cd $root/exp
    source ./atmos_bg_compile.csh
    if ( $status ) exit $status
endif

# Build the executable
set mkmf_exec = "$mkmf -f -m Makefile -a $code_dir -t $mkmfTemplate -p $executable:t"
mkdir -p $executable:h
cd $executable:h
if( $type == MOM_solo ) then
    set srcList = ( mom5/drivers )
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == ACCESS-OM || $type == ACCESS-CM ) then
    set srcList = ( access_coupler )
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == MOM_SIS ) then
    set srcList = ( coupler )
    set includes = "$includes -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_null -I$executable:h:h/lib_land_null"
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_null/lib_atmos_null.a $executable:h:h/lib_land_null/lib_land_null.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == EBM ) then
    set srcList = ( coupler )
    set includes = "$includes -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_ebm  -I$executable:h:h/lib_land_lad"
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_ebm/lib_atmos_ebm.a $executable:h:h/lib_land_lad/lib_land_lad.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == CM2M ) then
    set srcList = ( coupler )
    set includes = "$includes -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_fv -I$executable:h:h/lib_atmos_phys -I$executable:h:h/lib_land_lad"
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_fv/lib_atmos_fv.a $executable:h:h/lib_atmos_phys/lib_atmos_phys.a $executable:h:h/lib_land_lad/lib_land_lad.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == ESM2M ) then
    set srcList = ( coupler )
    set includes = "$includes -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_fv -I$executable:h:h/lib_atmos_phys -I$executable:h:h/lib_land_lad2"
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_fv/lib_atmos_fv.a $executable:h:h/lib_atmos_phys/lib_atmos_phys.a $executable:h:h/lib_land_lad2/lib_land_lad2.a $executable:h:h/lib_FMS/lib_FMS.a"
else if( $type == ICCM ) then
    set srcList = ( coupler )
    set includes = "$includes -I$executable:h:h/lib_ice -I$executable:h:h/lib_atmos_bg -I$executable:h:h/lib_atmos_phys -I$executable:h:h/lib_land_lad" 
    set libs = "$executable:h:h/lib_ocean/lib_ocean.a $executable:h:h/lib_ice/lib_ice.a $executable:h:h/lib_atmos_bg/lib_atmos_bg.a $executable:h:h/lib_atmos_phys/lib_atmos_phys.a $executable:h:h/lib_land_lad/lib_land_lad.a $executable:h:h/lib_FMS/lib_FMS.a"
endif
$mkmf_exec -o "$includes" -c "$cppDefs" -l "$libs"  $srcList
make
if( $status ) then
    echo "Make failed to create the $type executable"
    exit 1
endif    

exit
