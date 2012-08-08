# The following three directory may need to set.
NETCDFPATH = /usr/local/netcdf-3.6.2
SHAREDIR = $(PWD)/../../shared/mosaic
COREDIR  = $(PWD)
TOOLSHAREDIR = $(PWD)/../shared
CFLAGS = -O2 -fast -I$(TOOLSHAREDIR) -I$(COREDIR) -I$(SHAREDIR) -I${NETCDFPATH}/include -I/usr/include
LDFLAGS =  -L${NETCDFPATH}/lib -lm -lnetcdf 
LNFLAGS = -v
CC      = icc

OBJS = mosaic_util.o interp.o create_xgrid.o mpp.o mpp_domain.o mpp_io.o tool_util.o create_lonlat_grid.o \
       create_conformal_cubic_grid.o create_gnomonic_cubic_grid.o create_grid_from_file.o make_hgrid.o 

HEADERS = Makefile $(TOOLSHAREDIR)/mpp.h  $(TOOLSHAREDIR)/mpp_domain.h   \
          $(TOOLSHAREDIR)/mpp_io.h $(SHAREDIR)/mosaic_util.h  $(SHAREDIR)/create_xgrid.h  \
          $(SHAREDIR)/interp.h $(TOOLSHAREDIR)/tool_util.h               \
          $(COREDIR)/create_hgrid.h

make_hgrid: $(OBJS)
	$(CC) $(LNFLAGS) -o $@ $(OBJS) $(LDFLAGS)

make_hgrid.o: $(COREDIR)/make_hgrid.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(COREDIR)/make_hgrid.c 

mosaic_util.o: $(SHAREDIR)/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/mosaic_util.c

interp.o: $(SHAREDIR)/interp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/interp.c 

create_xgrid.o: $(SHAREDIR)/create_xgrid.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/create_xgrid.c

mpp.o: $(TOOLSHAREDIR)/mpp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp.c 

mpp_domain.o: $(TOOLSHAREDIR)/mpp_domain.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_domain.c 

mpp_io.o: $(TOOLSHAREDIR)/mpp_io.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_io.c

tool_util.o: $(TOOLSHAREDIR)/tool_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/tool_util.c

create_gnomonic_cubic_grid.o: $(COREDIR)/create_gnomonic_cubic_grid.c 
	$(CC) $(CFLAGS) -c $(COREDIR)/create_gnomonic_cubic_grid.c 

create_conformal_cubic_grid.o: $(COREDIR)/create_conformal_cubic_grid.c 
	$(CC) $(CFLAGS) -c $(COREDIR)/create_conformal_cubic_grid.c

create_lonlat_grid.o: $(COREDIR)/create_lonlat_grid.c 
	$(CC) $(CFLAGS) -c $(COREDIR)/create_lonlat_grid.c 

create_grid_from_file.o: $(COREDIR)/create_grid_from_file.c 
	$(CC) $(CFLAGS) -c $(COREDIR)/create_grid_from_file.c 
