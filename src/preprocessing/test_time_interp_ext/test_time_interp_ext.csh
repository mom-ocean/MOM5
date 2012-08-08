#!/bin/tcsh -f
#
# This is a TCSH program to test files used by time_interp_external_mod
# and/or data_override_mod.  This is intended to help debug runtime errors
# in FMS based codes where the above modules fail reading an INPUT file. 
# This is a quick way to make sure the file works without having to run the
# entire model.
#
# USER INPUT:
# filename = name of file to be tested
# fieldname = name of field in file
# year0 = initial year to use (typically start time of model)
# month0 = initial month to use
# day0 = initial day to use
# days_inc = increment in days
# ntime = number of time steps
# cal_type = calendar ('julian' , 'no_leap' or '360_day')
#
#
#  This preprocessing program was tested with Intel Fortran compiler on ia64 at GFDL.
#  In order to run on other system, some changes may be needed.
#######################################################################
# Users need to modify the following variables according to their need. 
# These particular values are here for testing purposes in GFDL only.  
#######################################################################
  set echo
  set platform     = "ncrc.intel"                                    
  set npes = 1
#
  set root         = $cwd:h:h:h                         # The directory that contains src/ and bin/
#
# Users must ensure the following  files exist for their platform.
#
  source $root/bin/environs.$platform                   # environment variables and compiler options

  set mkmfTemplate = $root/bin/mkmf.template.$platform  # path to mkmf template 

#############################################################################
# Users need not change anything below this line except the namelists values.
#############################################################################
  set mkmf         = $root/bin/mkmf                         
  set cppDefs      = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )        # list of cpp #defines to be passed to the source files
  set sharedir     = $root/src/shared                       # directory of the shared code.
  set includedir   = $sharedir/include
  set UTILITIES    = "$sharedir/{time_manager,fms,mpp,clocks,time_interp,axis_utils,platform,horiz_interp,constants,memutils,mosaic}"
  set srclist    = " $UTILITIES $sharedir/mpp/include $includedir "   
#
# Set the compiler options and environment variables. User must ensure the file exists for the platform.
#
  source $root/bin/environs.$platform

# compile the code and create executable
$mkmf -p test_time_interp_ext.exe -t $mkmfTemplate -c "$cppDefs -Dtest_time_interp_external" $srclist 
make  test_time_interp_ext.exe


cat > input.nml <<EOF
 &test_time_interp_ext_nml
 filename='foo.nc'
 fieldname='foo'
 year0=1900
 month0=1
 day0=1
 days_inc=1
 ntime=1
 cal_type='julian'
 /
EOF

$mpirunCommand $npes `./test_time_interp_ext.exe
unset echo
