NC=/usr/local/netcdf-3.6.1
CFLAGS=-g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE 
INC=-I$(NC)/include
LIB=-L$(NC)/lib -lnetcdf -lm

all:
	$(CC) xmalloc.c -c $(CFLAGS) $(INC)
	$(CC) strlist.c xmalloc.c -c $(CFLAGS) $(INC)
	$(CC) mppncscatter.c strlist.o xmalloc.o -o mppncscatter $(CFLAGS) $(INC) $(LIB) -static 

test:
	cd test
	../mppncscatter -x 4 -y 4 test3.nc
	../mppncscatter -x 4 -y 4 fv_core.res.tile1.nc
	../mppncscatter -x 2 -y 2 -X nx,nxp -Y ny,nyp C48_grid.tile1.nc

clean:
	rm -f *.o
	cd test; rm -f *.nc.????

