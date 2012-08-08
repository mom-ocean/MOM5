#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Makefile for LM3 data processing utilities
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# files: 
# combine-ncc -- combines lm3 netcdf restart files (in general, any compressed-
#    by-gathering set of files) into one.
#    v 2.0.1 -- a minor bug fix: number of input files was used incorrectly in
#       argument sanity check at the beginninhg of the program
#    v 2.0.2 -- fixed an error in inq_compressed_dim (nfu_compress.F90), where 
#       value of 'compress' attribute was not initialized to empty string every 
#       time before reading, resulting in possible confusion when there is more
#       than one compressed dimension per file
#    v.2.0.3 -- fixed an error where an unlimited dimension in input file was
#       converted to regular dimension of corresponding size in the output
#    v.2.0.4 -- introduced an ability to skip negative indices in compressed
#       dimension -- important for the cases when compressed index is zero
#       length (e.g., no cohorts in compute domain), but the file structure
#       must be preserved
#    v.2.0.5 -- limited minimum size of the compressed dimension to 1 -- important
#       in cases when the files don't have any valid compressed points in them, 
#       for example lake model restarts if there is no lakes at all. 
# decompress-ncc -- decompresses a compressed-by-gathering netcdf file into
#    a regular netcdf file
# gentest -- generates test data for restart combination utility
# nfu_c_test -- test of nfu compressed utilities

GARBAGE  = *.o *.mod *~ ferret.jnl* 
EXECUTABLES = combine-ncc decompress-ncc 
NETCDFPATH = /usr/local/netcdf-4.0.1

# ---- find out the system architecture ----------------------------------------
ARCH = $(shell uname -m)

ifeq ($(ARCH),ia64)
FC       = ifort
FFLAGS      = -r8 -O2 -g
INCLUDES    = -I. -I${NETCDFPATH}/include
LIBS        = -L${NETCDFPATH}/lib/shared -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lmpi -lz
F90         = ifort
CXX         = c++ 
CXXFLAGS    = -LANG:std 
endif

ifeq ($(ARCH),i686)
FC          = ifort
FFLAGS      = -r8
INCLUDES    = -I. -I/usr/local/include
LIBS        = -L/usr/local/lib -lnetcdf
F90         = ifort
CXX         = c++ 
CXXFLAGS    = -LANG:std 
endif

all: $(EXECUTABLES)

combine-ncc : combine-ncc.o nfu.o nfu_compress.o
	$(F90) $(FFLAGS) -o $@ combine-ncc.o nfu.o nfu_compress.o $(LIBS) 
decompress-ncc : decompress-ncc.o nfu.o nfu_compress.o
	$(F90) $(FFLAGS) -o $@ decompress-ncc.o nfu.o nfu_compress.o $(LIBS) 

# ------------------------ additional dependencies ----------------------------
average-res.o combine-ncc.o compress-ncc.o decompress-ncc.o \
   combine-res.o gentest.o nfu_c_test.o nfu_compress.o : nfu.o
combine-ncc.o compress-ncc.o decompress-ncc.o combine-res.o nfu_c_test.o : nfu_compress.o

%.o : %.F90
	$(F90) $(FFLAGS) $(INCLUDES) -c $<
%.o : %.f90
	$(F90) $(FFLAGS) $(INCLUDES) -c $<
%.o : %.f
	$(FC) $(FFLAGS) $(INCLUDES) -c $<
%.o : %.F
	$(FC) $(FFLAGS) $(INCLUDES) -c $<
%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<
%.o : %.c
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<

# ----------------------------- auxilary targets -------------------------------
.PHONY: clean realclean

clean:
	-rm -f $(GARBAGE)

realclean: clean
	-rm -f $(EXECUTABLES) ferret.jnl* *~
