#!/bin/tcsh -f
#=======================================================================
#      preprocessing : ocean_grid_generator
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscritp can be used to generate horizontal grid and/or vertical grid
#  and/or topography. It creates Makefile, compile if necessary and  saves output.
#  The generated data sets are grid_spec.nc (default value of namelist varible
#  output_file of grid_generator_nml) and are saved in the directory $work_dir.
#  When you want to generate topography (not idealized or special topography), a
#  topography source data file is needed, which is specified by the namelist
#  variable topog_file of topog_nml.
#  *******Important notice: You need to configure hgrid_nml, vgrid_nml, topog_nml
#  and grid_generator_nml, which are listed in this runscripts, to obtain desired grid.
#=======================================================================
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform     = "ncrc.intel"                                    # A unique identifier for your platform
  set npes         = 1                                      # number of processors
#
# input data file and destination grid
#
  set type         = "tripolar"                             # type of grid "tripolar" or "simple"
#
  set root         = $cwd:h:h:h:h                       # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 

  set topog_file   = $FMS_ARCHIVE/mom4/input_data/OCCAM_p5degree.nc

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################

  set name         = "ocean_grid"                           # name of the grid file will be generated
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include                          # fms include directory  
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/ocean_grid_generator.exe # executable created after compilation
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/ocean_grid_generator.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )        # list of cpp #defines to be passed to the source files

# list the source code
  
  set CORE         = "$tooldir/{grids_type.f90,grids_util.f90,ocean_grid_generator.f90}"
  set CORE         = "$CORE $tooldir/{hgrid.f90,vgrid.f90,topog.f90,check_mask.F90}"
  set UTILITIES    = "$sharedir/{axis_utils,constants,fms,mpp,horiz_interp,platform,memutils,mosaic}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir " 
  set srclist      =  ( $CORE $UTILITIES )                       # list of the source code

# compile the code and create executable
  if( ! -d $executable:h ) mkdir $executable:h
  cd $executable:h
  $mkmf -t $mkmfTemplate -p $executable:t -c "$cppDefs" $srclist 

  make $executable:t

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir

  cd $workdir

# get executable  
  cp $executable $executable:t

# --- set up namelist
  if($type == "simple") then
    cat >input.nml <<EOF
     &ocean_grid_generator_nml
       grid_type   = 'hgrid_vgrid_topog'
       output_file = '$name.nc'  /    
    &hgrid_nml
       nxlons=2,x_lon=0.,360.,dx_lon=2.5,2.5,
       nylats=2,y_lat=-90.,90.,dy_lat=2.0,2.0,
       tripolar_grid=.false.,lat_join=65, 
       debug = .true. /
    &vgrid_nml
       nzdepths=3,z_depth=0.0,100.0,5600.0,dz_depth=25.0,25.0,975.0
       debug = .true. /
    &topog_nml
       topography = 'rectangular_basin' 
       topog_depend_on_vgrid = .TRUE.
       fill_first_row = .true.
       kmt_min=2,
       filter_topog=.false.,
       num_filter_pass=5,
       scale_factor=-1,
       interp_method = "conservative"
       debug = .true. /
    &fms_nml
       domains_stack_size = 105560 /

EOF

  else if($type == "tripolar") then
    cat >input.nml <<EOF
     &ocean_grid_generator_nml
       grid_type   = 'hgrid_vgrid_topog'
       output_file = '$name.nc'  /    
    &hgrid_nml
       nxlons=2,x_lon=-280.,80.,dx_lon=1.0,1.0,
       nylats=7,y_lat=-82.0,-30,-10.,0.,10.,30.0,90.
       dy_lat= 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0,
       tripolar_grid=.true.,lat_join=65, 
       debug = .true. /
    &vgrid_nml
       nzdepths=3,z_depth=0.,220.,5500.,dz_depth=10.,10.,367.14286
       debug = .true. /
    &topog_nml
       topography = 'from_file' 
       topog_depend_on_vgrid = .TRUE.
       topog_file='$topog_file',
       topog_field='topo',
       lon_field = 'geolon_t'
       lat_field = 'geolat_t'
       fill_first_row = .true.
       num_nbrs = 4
       kmt_min=4,
       max_dist = 0.05,
       filter_topog=.false.,
       adjust_topo = .true.,
       num_filter_pass=5,
       scale_factor= -1,
       src_is_spherical = .true.
       interp_method = "bilinear"
       debug = .true. /
    &fms_nml
       domains_stack_size = 576000 /

EOF

  else
	echo Unknown grid type $type
	exit 1
  endif


#  run the executable
  $mpirunCommand $npes ./$executable:t >fms.out
  
  cat fms.out     

#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo

