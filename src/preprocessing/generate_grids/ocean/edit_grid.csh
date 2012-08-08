#!/bin/tcsh 
#=======================================================================
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscritp can can edit the topography of input grid_spec file "orig_grid" 
#  according to the ascii input file "grid_edits". Then it will output the 
#  new grid_spec file "mod_grid". 
#=======================================================================
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################
  set echo
  set platform     = "ncrc.intel"                      # A unique identifier for your platform
#
  set root         = $cwd:h:h:h:h                       # The directory that contains src/ and bin/
  set npes         = 1
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set name         = "edit_grid"                       # name of the grid file will be generated
  set tooldir      = $cwd                              # directory of the tool
  set sharedir     = $root/src/shared                  # directory of the shared code.
  set includedir   = $sharedir/include                 # fms include directory   
  set workdir      = $tooldir/workdir                  # where the tool is run and output is produced
  set executable   = $tooldir/exec/edit_grid.exe       # executable created after compilation
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/edit_grid.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                    # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )                # list of cpp #defines to be passed to the source files

# grid to be edited and grid_edits test file
  set grid_edits = grid_edits.txt                      # text file to specify the edit region and new depth.
  set orig_grid  = $workdir/ocean_grid.nc         

# list the source code
  
  set CORE      = " $tooldir/{edit_grid.F90,topog.f90,grids_type.f90,grids_util.f90,check_mask.F90} "
  set UTILITIES    = "$sharedir/{axis_utils,constants,fms,mpp,horiz_interp,platform,memutils,mosaic}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir " 
  set srclist   = ( $CORE $UTILITIES )


# compile the model code and create executable
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

# if the grid_edits file does not exist, create one here.
  if( ! -f $grid_edits) then
     cat >$grid_edits <<EOF
     100:105, 50:65, 1000   #lon_start:lon_end, lat_start:lat_end, new height.
     65, 75, 200
     5,16, 0
EOF
  endif

# --- set up namelist  

  cat >input.nml <<!
    &edit_grid_nml
       orig_grid   = '$orig_grid' 
       mod_grid    = '$name.nc'
       grid_edits  = '$grid_edits'
       /
!

  $mpirunCommand $npes ./$executable:t >fms.out

  cat fms.out

# rename ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out

  unset echo
