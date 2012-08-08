# The following three directory may need to set.
NETCDFPATH = /usr/local/netcdf-4.0.1
POSTPPATH  = .
FFLAGS = -I${NETCDFPATH}/include -I/usr/include -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -r8 -nowarn -Wp,-w
LDFLAGS =  -L${NETCDFPATH}/lib/shared -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lmpi -lz
DEFFLAG = -Duse_netCDF
LNFLAGS = -v
FC      = ifort

OBJS = list_ncvars.o

list_ncvars.exe: $(OBJS)
	$(FC) $(LNFLAGS) -o $@ $(OBJS) $(LDFLAGS)

list_ncvars.o: $(POSTPPATH)/list_ncvars.f90
	$(FC) $(FFLAGS) -c $(POSTPPATH)/list_ncvars.f90
