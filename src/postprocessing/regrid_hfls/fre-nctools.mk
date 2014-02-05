#
# $Id: fre-nctools.mk,v 20.0 2013/12/14 00:29:57 fms Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build Regridding Executables
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Initial version (Makefile, ver 17.0.4.2)       June 10
# afy    Ver   1.01  Add rules to build MPI-based executable        June 10
# afy    Ver   1.02  Simplified according to fre-nctools standards  June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2010
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#

FC       := ifort
FFLAGS   := -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -r8 -traceback
INCLUDES := -I${NETCDF_HOME}/include
LIBS     := -L${NETCDF_HOME}/lib/shared -L${HDF5_HOME}/lib/shared -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -limf

TARGETS  := regrid_hfls.exe

SOURCES  := regrid_hfls.F90

OBJECTS  := $(SOURCES:F90=o)

HEADERS = fre-nctools.mk

all: $(TARGETS)

regrid_hfls.exe: regrid_hfls.o
	$(FC) -o $@ $^ $(OBJECTS) $(LIBS)

regrid_hfls.o: regrid_hfls.F90 $(HEADERS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)
