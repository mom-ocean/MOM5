#
# $Id: fre-nctools.mk,v 20.0 2013/12/14 00:29:08 fms Exp $
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
FFLAGS   := -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -r8 -traceback
INCLUDES := -I${NETCDF_HOME}/include
LIBS     :=  -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib $(LIBNETCDFF) -lnetcdf -lhdf5_hl -lhdf5 -lz -limf $(LIBS2) $(STATIC)

TARGETS  := list_ncvars.exe

SOURCES  := list_ncvars.f90

OBJECTS  := $(SOURCES:f90=o)

HEADERS = fre-nctools.mk

all: $(TARGETS)

list_ncvars.exe: list_ncvars.o
	$(FC) -o $@ $(OBJECTS) $(LIBS)

list_ncvars.o: list_ncvars.f90 $(HEADERS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)
