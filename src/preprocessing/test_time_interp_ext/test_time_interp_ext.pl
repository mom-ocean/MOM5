#!/usr/bin/perl

use Getopt::Long;
use Cwd;


#-------------------------------------------------------------------------------
#get command-line arguments

my $filename;
my $fieldname;
my $ntime=12;
my $year0=1990;
my $month0=1;
my $day0=1;
my $inc_days=31;
my $NPES = 2;
my $help = 0;
my $nbuf = -1;
my $cal_type = 'julian';

GetOptions("f=s" => \$filename , "v=s" => \$fieldname, "nt=i" => \$ntime, "npes=i" => \$NPES, "year0=i" => \$year0, "month0=i" => \$month0, "day0=i" => \$day0, "inc_days=i" => \$inc_days, "nbuf=i" => \$nbuf, "cal_type=s" => \$cal_type, "h" => \$help);

exit print "Usage: test_time_interp_ext -f filename -v variable [--nt ntimes] [--npes #] [--year0 ] [--month0 ] [--day0 ] [--inc_days ]
	   \n"
           
if($help || ! defined $filename || ! defined $fieldname);


open(FH, ">input.nml") || die print "cannot open file\n";
print FH " &test_time_interp_ext_nml\n";
print FH "filename='$filename',\n";
print FH "fieldname='$fieldname',\n";
print FH "year0=$year0,\n";
print FH "month0=$month0,\n";
print FH "day0=$day0,\n";
print FH "days_inc=$inc_days,\n";
print FH "ntime=$ntime\n";
print FH "cal_type='$cal_type',\n";
print FH "/\n";
print FH;
print FH " &time_interp_external_nml\n";
print FH "num_io_buffers=$nbuf";
print FH "/";
close FH;

my $current_directory = cwd;

my $status = system ("mpirun -np $NPES $ENV{FMS_BIN_DIR}/test_time_interp_ext.exe  2> test_error");

if ($status != 0) {
print "\nTest failed.  Examine output and test_error files to diagnose\n";exit}


system("rm test_error");


