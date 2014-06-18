#!/bin/tcsh -f
#
# $Id: timavg.csh,v 20.0 2013/12/14 00:30:18 fms Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Script to Call Time-Averaging Executables
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Copied from ~fms/local/ia64/netcdf4.fix        June 10
# afy    Ver   2.00  Don't source the 'init.csh' script             June 10
# afy    Ver   2.01  Use 'which' to locate executables              June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2000-2010
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#

#-----------------------------------------------------------------------
#
#  This script takes multiple netcdf files that all contain the same
#  dimensions and variables but may have a different number of time
#  records and writes a new netcdf file that has a time record 
#  averaged from each input file.
#
#-----------------------------------------------------------------------

set ofile
set ifiles
set etime = .true.
set do_verbose = .false.
set do_bounds  = .false.
set do_errors  = .false.
set no_warning = .false.
set vers = v08
set precision = 8

#  ----- parse input argument list ------

set argv = (`getopt abdmWo:v:w:r: $*`)

while ("$argv[1]" != "--")
    switch ($argv[1])
        case -a:
            set do_errors = .true.; breaksw
        case -b:
            set do_bounds = .true.; breaksw
        case -d:
            set debug; set do_verbose = .true.; breaksw
        case -m:
            set etime = .false.; breaksw
        case -r:
            set precision = $argv[2]; shift argv; breaksw
        case -W:
            set no_warning = .true.; breaksw
        case -w:
            set weight = $argv[2]; shift argv; breaksw
        case -o:
            set ofile = $argv[2]; shift argv; breaksw
        case -v:
            set vers = $argv[2]; shift argv; breaksw
    endsw
    shift argv
end
shift argv
if ($?debug) set echo
set ifiles = ( $argv )

##################################################################
#  ----- help message -----

if ($ofile == "" || $#ifiles == 0) then
set name = `basename $0`
cat << EOF

Time averaging script

Usage:  $name [-a] [-b] [-d] [-m] [-r prec] [-v vers] -o ofile  files.....

        -a       = skips "average information does not agree" errors
        -b       = adds time axis bounds and cell methods (CF convention)
        -d       = turns on command echo (for debugging)
        -m       = average time (instead of end time) for t-axis values
        -r prec  = precision used for time averaging, either -r4 or -r8 (default)
        -W       = suppress warning messages (use with caution)
        -w wght  = minimum fraction of missing data needed for valid data
        -o ofile = name of the output file

        files... = list netcdf files, each file will be a time record
                   in the output file (the files should be in 
                   chronological order)
                  
EOF
exit 1
endif

#        -v vers  = executable version (TAVG.$vers.exe)
##################################################################

# executable name depends on precision 
if ($precision == 4) then
    set executable = `which TAVG.r4.exe`
else if ($precision == 8) then
    set executable = `which TAVG.exe`
else
    echo "ERROR: use -r4 or -r8"; exit 1
endif

set error_flag = 0

#-- check existence of executable --
if (! -e $executable) then
   echo Executable does not exist
   set error_flag = 1
endif

#-- check existence of files --
foreach file ($ifiles)
   if (! -e $file) then
       echo File $file does not exist
       set error_flag = 1
   endif
end
if ($error_flag != 0) then
   echo ERROR: input files do not exist
   exit 1
endif

#-- namelist (create unique name) --
set nml_name = nml`date '+%j%H%M%S'`
    echo " &input" > $nml_name
set i = 1
foreach file ($ifiles)
    echo "    file_names("$i") = " \'$file\' , >> $nml_name
    @ i++
end
    echo "    file_name_out = " \'$ofile\' ,     >> $nml_name
    echo "    use_end_time  = " $etime ,         >> $nml_name
    echo "    verbose  = " $do_verbose ,         >> $nml_name
    echo "    add_cell_methods = " $do_bounds ,  >> $nml_name
    echo "    skip_tavg_errors = " $do_errors ,  >> $nml_name
    echo "    suppress_warnings = " $no_warning, >> $nml_name

    if ($?weight) then
       echo "    frac_valid_data = " $weight ,   >> $nml_name
    endif

    echo " &end"                                 >> $nml_name

#-- run averaging program --
    $executable < $nml_name
    set exit_status = $status

#-- clean up --
rm -f $nml_name

exit $exit_status

