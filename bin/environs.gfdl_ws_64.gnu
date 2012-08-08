         source $MODULESHOME/init/csh
         module purge
	 module rm netcdf hdf5
         module load mpich2/1.2.1p1
         module use -a /home/fms/local/modulefiles
         module load hdf5/1.8.5-patch1-gnu-4
         module load netcdf/4.1.1-gnu-4
#
         setenv PATH ${PATH}:.
         setenv mpirunCommand   "/net2/nnz/opt/mpich2-1.3_ifort11_x64/bin/mpirun -np"
         setenv FMS_ARCHIVE /archive/fms
         setenv PATH ${PATH}:.


