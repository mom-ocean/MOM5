#!/bin/tcsh 
#=======================================================================
#      idealized_bc : generate idealized boundary condition
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscritp can be used to generate idealized surface boundary condition
#  for mom4. It creates Makefile, compile if necessary and  saves output.The generated
#  data sets are saved in the directory $output_dir/data. The generated data sets will depend
#  on the namelist boundary_nml. A grid_spec file is needed for this program.
#  ******Important notice: to obtain desired initial condition, you need to
#  specify boundary_nml which is listed in this runscript.
#=======================================================================
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform     = "ncrc.intel"                           # A unique identifier for your platform
  set npes         = 1                                      # number of processors
#
# input data file and destination grid
#
  set grid_spec_file  = $cwd:h:h/generate_grids/ocean/workdir/ocean_grid.nc  # destination grid
#
  set root         = $cwd:h:h:h:h                         # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set name         = "idealized_bc"                         # name of the tool
  set tooldir      = $cwd                                   # directory of the tool
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include                      # fms include directory
  set workdir      = $tooldir/workdir                       # where the tool is run and output is produced
  set executable   = $tooldir/exec/idealized_bc.exe            # executable created after compilation
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/idealized_bc.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                         # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )        # list of cpp #defines to be passed to the source files

# list the source code
  set CORE      = $tooldir/{idealized_bc.f90,idealized_bc_driver.f90}
  set UTILITIES = "$sharedir/{fms,mpp,platform,constants,memutils}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir "
  set srclist   = ( $CORE $UTILITIES )


# compile the code and create executable
  if( ! -d $executable:h ) mkdir $executable:h
  cd $executable:h
  $mkmf -t $mkmfTemplate -p $executable:t -c "$cppDefs" $srclist 

  make $executable:t

#--------------------------------------------------------------------------------------------------------
# setup directory structure
  if ( ! -d $workdir )         mkdir $workdir
  if ( ! -d $workdir/INPUT )   mkdir $workdir/INPUT

  cd $workdir

# get executable  
  cp $executable $executable:t

# get grid file
  cp $grid_spec_file INPUT/grid_spec.nc

# --- set up namelist

    cat >input.nml <<!
    &idealized_bc_nml
       wind_type = 'frank_bryan_winds_compress' 
/
    &fms_nml
       domains_stack_size = 5000000
/
!


#  run the executable
  $mpirunCommand $npes $executable:t >fms.out
  
  cat fms.out
       
#---------------------------------------------------------------------
# rename the ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo

