#!/bin/tcsh 
#=======================================================================
#      compare_grid : compare topography of two grid files.
#         Contact : Zhi Liang   email : Zhi.Liang
#
#  This runscritp compares two grid descriptor files (generated via ocean_grid_generator) 
#  and creates a text file output listing line-by-line differences between the
#  two files. Output file format is the same as the grid_edits file used by
#  edit_grid.F90. These two compared files should have same grid size and same grid
#  resolution.
#=======================================================================
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################

  set echo
  set platform     = "ncrc.intel"                              # A unique identifier for your platform
#
  set root         = $cwd:h:h:h:h                       # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 
  # grid data will be compared
  set grid_file  = $FMS_ARCHIVE/mom4/mom4p0/mom4p0c/mom4_test5/preprocessing/grid_spec_v7.nc
#
#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set name         = "check_mask"                      # name of tool and name of the output
  set tooldir      = $cwd                              # directory of the tool
  set sharedir     = $root/src/shared                  # directory of the shared code.   
  set includedir   = $sharedir/include
  set workdir      = $tooldir/workdir                  # where the tool is run and output is produced
  set executable   = $tooldir/exec/check_mask.exe      # executable created after compilation
  if ( $?DEBUG ) then
    set executable    = $tooldir/debug/check_mask.exe 
    set mkmfTemplate  = $root/bin/mkmf.debugtemplate.$platform         
  endif
  set mkmf         = $root/bin/mkmf                    # path to executable mkmf
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI -Dcheck_mask" )                # list of cpp #defines to be passed to the source files

# list the source code
  
  set CORE      = "$tooldir/check_mask.F90"
  set UTILITIES    = "$sharedir/{constants,fms,mpp,platform,memutils}"
  set UTILITIES    = " $UTILITIES $sharedir/mpp/include $includedir" 
  set srclist   = "$CORE $UTILITIES"

#---------------------------------------------------------------------
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

# --- set up namelist

    cat >input.nml <<!
    &check_mask_driver_nml
       grid_file   = '$grid_file' 
       /
!

  $executable:t >fms.out
  cat fms.out
#   --- rename ascii output file
  mv fms.out $name.fms.out
  mv logfile*.out $name.logfile.out
  
  unset echo  
