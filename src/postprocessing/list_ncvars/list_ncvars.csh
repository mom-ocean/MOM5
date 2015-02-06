#!/bin/tcsh -f
#
# $Id: list_ncvars.csh,v 20.0 2013/12/14 00:29:09 fms Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Script to Call the "list_ncvars" Program
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Initial version (copied from 1.1.2.1)          June 10
# afy    Ver   2.00  Don't source the "init.csh" script             June 10
# afy    Ver   2.01  Use "which" to locate the "list_ncvars"        June 10
# afy    Ver   3.00  Return to the original executable name         June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2000-2010
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#

set static    = .false.
set nonstatic = .false.
set var0d  = .false.
set var1d  = .false.
set var2d  = .false.
set var3d  = .false.
set var4d  = .false.

##################################################################
#  ----- parse input argument list ------

set argv = (`getopt st01234 $*`)

while ("$argv[1]" != "--")
    switch ($argv[1])
        case -s:
            set static = .true.; breaksw
        case -t:
            set nonstatic = .true.; breaksw
        case -0:
            set var0d = .true.; breaksw
        case -1:
            set var1d = .true.; breaksw
        case -2:
            set var2d = .true.; breaksw
        case -3:
            set var3d = .true.; breaksw
        case -4:
            set var4d = .true.; breaksw
    endsw
    shift argv
end
shift argv
set list = ( $argv )

# usage message?
if ($#list == 0) set help
if ($static == ".false." && $nonstatic == ".false.") set help
if ($var0d == ".false." && $var1d == ".false." && \
    $var2d == ".false." && $var3d == ".false." && $var4d == ".false.") set help

##################################################################
#----- help message -----

if ($?help) then
set name = `basename $0`
cat << EOF

Usage:  $name -s | -t [-0123] files...

        -s   list names of static fields (no time dimension)
        -t   list names of non-static fields (with time dimension)
        -0   list names of scalar fields
        -1   list names of 1-d fields
        -2   list names of 2-d fields
        -3   list names of 3-d fields
        -4   list names of 4-d fields
        files... list of netcdf files

   Options -0,1,2,3,4 exclude the time dimension.
   Must specify either "-s" or "-t".
   Also must specify one (or more) of the following: -0, -1, -2, -3, -4.
   A file name must be supplied.
   Variables that are also dimensions are skipped.
   Variables that end in "_T1", "_T2", and "_DT" are skipped.

Example:  $name -st -0123 myfile.nc
          
          This will list all variables (static and non-static) that
          have a rank of 0, 1, 2, or 3.

EOF
exit 1
endif

##################################################################

alias list_ncvars `which list_ncvars.exe`

# unique namelist name
set nml_name = nml`date '+%j%H%M%S'`

#-------------------------
# loop thru netcdf files
foreach file ($list)

cat > $nml_name << EOF
 &input
   filename = '$file',
   list_static = $static,
   list_nonstatic = $nonstatic,
   var0d = $var0d,
   var1d = $var1d,
   var2d = $var2d,
   var3d = $var3d,
   var4d = $var4d,
 &end
EOF

list_ncvars < $nml_name

end
#-------------------------

# clean up
rm -f $nml_name

