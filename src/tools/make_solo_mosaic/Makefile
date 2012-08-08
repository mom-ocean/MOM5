# The following three directory may need to set.
NETCDFPATH = /usr/local/netcdf-3.6.2
SHAREDIR = $(PWD)/../../shared/mosaic
COREDIR  = $(PWD)
TOOLSHAREDIR = $(PWD)/../shared
CFLAGS = -O2 -fast -I$(TOOLSHAREDIR) -I$(COREDIR) -I$(SHAREDIR) -I${NETCDFPATH}/include -I/usr/include
LDFLAGS =  -L${NETCDFPATH}/lib -lm -lnetcdf 
LNFLAGS = -v
CC      = icc

OBJS = mosaic_util.o mpp.o mpp_domain.o mpp_io.o tool_util.o get_contact.o make_solo_mosaic.o \
       interp.o create_xgrid.o 

HEADERS = Makefile $(TOOLSHAREDIR)/mpp.h  $(TOOLSHAREDIR)/mpp_domain.h   \
          $(TOOLSHAREDIR)/mpp_io.h $(TOOLSHAREDIR)/tool_util.h               \
          $(COREDIR)/get_contact.h $(SHAREDIR)/interp.h  $(SHAREDIR)/create_xgrid.h

make_solo_mosaic: $(OBJS)
	$(CC) $(LNFLAGS) -o $@ $(OBJS) $(LDFLAGS)

make_solo_mosaic.o: $(COREDIR)/make_solo_mosaic.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(COREDIR)/make_solo_mosaic.c 

mpp.o: $(TOOLSHAREDIR)/mpp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp.c 

mpp_domain.o: $(TOOLSHAREDIR)/mpp_domain.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_domain.c 

mpp_io.o: $(TOOLSHAREDIR)/mpp_io.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_io.c

mosaic_util.o: $(SHAREDIR)/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/mosaic_util.c

tool_util.o: $(TOOLSHAREDIR)/tool_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/tool_util.c

get_contact.o: $(COREDIR)/get_contact.c 
	$(CC) $(CFLAGS) -c $(COREDIR)/get_contact.c 

interp.o: $(SHAREDIR)/interp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/interp.c 

create_xgrid.o: $(SHAREDIR)/create_xgrid.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/create_xgrid.c
