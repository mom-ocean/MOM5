#
# $Id: fre-nctools.mk,v 20.0 2013/12/14 00:30:17 fms Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build Regridding Executables
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Initial version (Makefile, ver 17.0.4.2)       June 10
# afy    Ver   1.01  Add rules to build MPI-based executable        June 10
# afy    Ver   1.02  Simplified according to fre-nctools standards  June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
#
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#
include ./env.$(SITE)

FC       := ifort
FFLAGS   := -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -real_size 64 -traceback
FFLAGS_r4:= -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -real_size 32 -traceback
INCLUDES := -I${NETCDF_HOME}/include
LIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib $(LIBNETCDFF) -lnetcdf -lhdf5_hl -lhdf5 -lz -limf $(LIBS2) $(STATIC)

TARGETS  := TAVG.exe TAVG.r4.exe

OBJECTS  := time_average.o
OBJECTS_r4= time_average.r4.o

HEADERS = fre-nctools.mk

all: $(TARGETS)

TAVG.exe: time_average.o
	$(FC) -o $@ $(OBJECTS) $(LIBS)

time_average.o: time_average.f90 $(HEADERS)
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $< 

TAVG.r4.exe: time_average.r4.o
	$(FC) -o $@ $(OBJECTS_r4) $(LIBS)

time_average.r4.o: time_average.f90 $(HEADERS)
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $< 

%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<

%.r4..o: %.f90
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $<

clean:
	-rm -f *.o *.mod $(TARGETS)
