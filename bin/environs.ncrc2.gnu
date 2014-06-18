          source $MODULESHOME/init/csh
          module use -a /ncrc/home2/fms/local/modulefiles
          module unload PrgEnv-pgi PrgEnv-pathscale PrgEnv-intel PrgEnv-gnu PrgEnv-cray
          module unload netcdf fre fre-commands
          module load PrgEnv-gnu
          module load hdf5/1.8.8
          module load netcdf/4.2.0
          module list
          setenv MPICH_MAX_SHORT_MSG_SIZE 8000
          setenv KMP_STACKSIZE 512m
          setenv NC_BLKSZ 1M

        setenv mpirunCommand   "aprun -n"
        setenv PATH ${PATH}:.

