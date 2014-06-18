#
# $Id: fre-nctools.mk,v 20.0 2013/12/14 00:29:24 fms Exp $
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
CFLAGS_O2:= -O2 -g -traceback
INCLUDES := -I${NETCDF_HOME}/include
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz -limf $(CLIBS2) $(STATIC)

TARGETS  := mppnccombine

SOURCES  := mppnccombine.c

OBJECTS  := $(SOURCES:c=o)

HEADERS = fre-nctools.mk

all: $(TARGETS)

mppnccombine: $(OBJECTS)
	$(CC) -o $@ $^ $(CLIBS)

mppnccombine.o: mppnccombine.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)
