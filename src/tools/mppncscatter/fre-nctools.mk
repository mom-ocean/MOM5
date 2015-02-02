#
# $Id: fre-nctools.mk,v 19.4 2013/12/04 17:14:04 fms Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build Regridding Executables
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
CLIBS    := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz -limf $(CLIBS2) $(STATIC)

TARGETS  := mppncscatter

SOURCES  := mppncscatter.c opt.c scatterdim.c strlist.c xmalloc.c

OBJECTS  := $(SOURCES:c=o)

HEADERS = fre-nctools.mk

all: $(TARGETS)

mppncscatter: $(OBJECTS) main.c
	$(CC) -o $@ $^ $(CLIBS) $(INCLUDES)

mppncscatter.o: mppncscatter.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)

