#PBS -P v45 
#PBS -l wd,ncpus=960,mem=1920Gb,walltime=5:00:00 
#PBS -W block=true

module purge

module use /g/data3/hh5/public/modules
module load conda/analysis3-unstable
module load pbs

nosetests --with-xunit test_run.py -s
