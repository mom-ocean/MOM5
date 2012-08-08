# The following three directory may need to set.
NETCDFPATH = /usr/local/netcdf-4.0.1
SHAREDIR = $(PWD)/../../shared/mosaic
COREDIR  = $(PWD)
TOOLSHAREDIR = $(PWD)/../shared
CFLAGS = -O3 -I$(TOOLSHAREDIR) -I$(COREDIR) -I$(SHAREDIR) -I${NETCDFPATH}/include
CFLAGS_O2 = -O2 -I$(TOOLSHAREDIR) -I$(COREDIR) -I$(SHAREDIR) -I${NETCDFPATH}/include
LDFLAGS =  -L${NETCDFPATH}/lib/shared -lm -lnetcdf -lhdf5_hl -lhdf5 -lmpi -lz
DEFFLAG = -Duse_netCDF
LNFLAGS = -v
CC      = icc

OBJS = mosaic_util.o gradient_c2l.o create_xgrid.o interp.o read_mosaic.o mpp.o mpp_domain.o mpp_io.o tool_util.o \
       conserve_interp.o bilinear_interp.o fregrid_util.o fregrid.o

HEADERS = Makefile $(TOOLSHAREDIR)/mpp.h  $(TOOLSHAREDIR)/mpp_domain.h   \
          $(TOOLSHAREDIR)/mpp_io.h $(SHAREDIR)/mosaic_util.h                    \
          $(SHAREDIR)/interp.h $(SHAREDIR)/create_xgrid.h
fregrid: $(OBJS)
	$(CC) $(LNFLAGS) -o $@ $(OBJS) $(LDFLAGS)

fregrid.o: $(COREDIR)/fregrid.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(COREDIR)/fregrid.c 

fregrid_util.o: $(COREDIR)/fregrid_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(COREDIR)/fregrid_util.c 

conserve_interp.o: $(COREDIR)/conserve_interp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(COREDIR)/conserve_interp.c

bilinear_interp.o: $(COREDIR)/bilinear_interp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(COREDIR)/bilinear_interp.c

mosaic_util.o: $(SHAREDIR)/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/mosaic_util.c

gradient_c2l.o: $(SHAREDIR)/gradient_c2l.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/gradient_c2l.c

create_xgrid.o: $(SHAREDIR)/create_xgrid.c $(HEADERS)
	$(CC) $(CFLAGS_O2) -c $(SHAREDIR)/create_xgrid.c 

interp.o: $(SHAREDIR)/interp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/interp.c 

read_mosaic.o: $(SHAREDIR)/read_mosaic.c $(HEADERS)
	$(CC) $(DEFFLAG) $(CFLAGS) -c $(SHAREDIR)/read_mosaic.c 

mpp.o: $(TOOLSHAREDIR)/mpp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp.c 

mpp_domain.o: $(TOOLSHAREDIR)/mpp_domain.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_domain.c 

mpp_io.o: $(TOOLSHAREDIR)/mpp_io.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_io.c

tool_util.o: $(TOOLSHAREDIR)/tool_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/tool_util.c


