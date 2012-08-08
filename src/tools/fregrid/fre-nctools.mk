#
# $Id: fre-nctools.mk,v 18.0 2010/07/28 17:44:11 fms Exp $
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

CC       := icc
CFLAGS   := -O3
CFLAGS_O2:= -O2
INCLUDES := -I${NETCDF_HOME}/include
LIBS     := -L${NETCDF_HOME}/lib/shared -L${HDF5_HOME}/lib/shared -lnetcdf -lhdf5_hl -lhdf5 -lmpi -lz

TARGETS  := fregrid fregrid_parallel

SOURCES  := fregrid.c bilinear_interp.c conserve_interp.c fregrid_util.c
SOURCES  += create_xgrid.c gradient_c2l.c interp.c mosaic_util.c
SOURCES  += mpp_domain.c mpp_io.c tool_util.c

OBJECTS  := $(SOURCES:c=o)

all: $(TARGETS)

fregrid: $(OBJECTS) read_mosaic.o mpp.o
	$(CC) -o $@ $^ $(LIBS)

fregrid_parallel: $(OBJECTS) read_mosaic_parallel.o mpp_parallel.o
	$(CC) -o $@ $^ $(LIBS)

read_mosaic.o: read_mosaic.c
	$(CC) -Duse_netCDF $(CFLAGS) $(INCLUDES) -c $< 

read_mosaic_parallel.o: read_mosaic.c
	$(CC) -Duse_netCDF -Duse_libMPI $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp_parallel.o: mpp.c
	$(CC) -Duse_libMPI $(CFLAGS) $(INCLUDES) -o $@ -c $< 

create_xgrid.o: create_xgrid.c
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)
