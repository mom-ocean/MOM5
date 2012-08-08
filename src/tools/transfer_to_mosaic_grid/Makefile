# The following three directory may need to set.
NETCDFPATH = /usr/local/netcdf-3.6.2
SHAREDIR = $(PWD)/../../shared/mosaic
COREDIR  = $(PWD)
TOOLSHAREDIR = $(PWD)/../shared
#CFLAGS = -O0 -g -I$(TOOLSHAREDIR) -I$(COREDIR) -I$(SHAREDIR) -I${NETCDFPATH}/include -I/usr/include
CFLAGS = -O2 -fast -I$(TOOLSHAREDIR) -I$(COREDIR) -I$(SHAREDIR) -I${NETCDFPATH}/include -I/usr/include
LDFLAGS =  -L${NETCDFPATH}/lib -lm -lnetcdf 
LNFLAGS = -v
CC      = icc

OBJSM = mpp.o mpp_domain.o mpp_io.o transfer_to_mosaic.o create_xgrid.o mosaic_util.o

HEADERS = Makefile $(TOOLSHAREDIR)/mpp.h  $(TOOLSHAREDIR)/mpp_domain.h   \
          $(TOOLSHAREDIR)/mpp_io.h $(SHAREDIR)/create_xgrid.h  \
	  $(SHAREDIR)/mosaic_util.h

transfer_to_mosaic_grid: $(OBJSM)
	$(CC) $(LNFLAGS) -o $@ $(OBJSM) $(LDFLAGS)

run: 
	transfer_to_mosaic
clean: 
	rm -f *.o transfer_to_mosaic

transfer_to_mosaic.o: $(COREDIR)/transfer_to_mosaic.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(COREDIR)/transfer_to_mosaic.c 

mpp.o: $(TOOLSHAREDIR)/mpp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp.c 

mpp_domain.o: $(TOOLSHAREDIR)/mpp_domain.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_domain.c 

mpp_io.o: $(TOOLSHAREDIR)/mpp_io.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_io.c

create_xgrid.o: $(SHAREDIR)/create_xgrid.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/create_xgrid.c

tool_util.o: $(TOOLSHAREDIR)/tool_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/tool_util.c

mosaic_util.o: $(SHAREDIR)/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/mosaic_util.c
