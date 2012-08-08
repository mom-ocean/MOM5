#!/bin/tcsh

switch (`uname -m`)
case 'ia64':
   source /opt/modules/default/init/tcsh
   module purge
   module load ifort.9.1.041
   module load netcdf-3.6.2
   breaksw
case 'i686':
   breaksw
default:
   echo Unsupported architecture `uname -m` 
endsw

module list

make
