#
# $Id: fre-nctools.mk,v 20.0 2013/12/14 00:28:52 fms Exp $
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

CC       := icc
CFLAGS   := -O3 -g -traceback
CLIBS    := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz -limf $(CLIBS2) $(STATIC)

FC       := ifort
FFLAGS   := -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -real_size 64 -traceback
INCLUDES := -I${NETCDF_HOME}/include
LIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib $(LIBNETCDFF) -lnetcdf -lhdf5_hl -lhdf5 -lz -limf $(LIBS2) $(STATIC)

TARGETS  := combine-ncc decompress-ncc is-compressed

SOURCES  := nfu.F90 nfu_compress.F90

OBJECTS  := $(SOURCES:F90=o)

HEADERS = fre-nctools.mk

all: $(TARGETS)

combine-ncc: combine-ncc.o
	$(FC) -o $@ $^ $(OBJECTS) $(LIBS)

combine-ncc.o: combine-ncc.F90 $(OBJECTS) $(HEADERS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

decompress-ncc: decompress-ncc.o
	$(FC) -o $@ $^ $(OBJECTS) $(LIBS)

decompress-ncc.o: decompress-ncc.F90 $(OBJECTS) $(HEADERS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

nfu_compress.o: nfu_compress.F90 nfu.o $(HEADERS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

nfu.o: nfu.F90 $(HEADERS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

is-compressed: is-compressed.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(CLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o *.mod $(TARGETS)
