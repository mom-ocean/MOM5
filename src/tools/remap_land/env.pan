# GFDL did not build shared libs for 4.1.1. So we have to add explicit link to libcurl
LIBS2  := -lcurl
CLIBS2 := -lcurl

# GFDL uses the mpicc wrapper and Intel icc
MPICC    := mpicc
CC       := icc
