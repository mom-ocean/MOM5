#PBS -P v45
#PBS -l ncpus=960
#PBS -l mem=1920Gb
#PBS -l walltime=5:00:00 
#PBS -l wd
#PBS -l storage=gdata/hh5+scratch/v45
#PBS -W block=true

module purge

module use /g/data3/hh5/public/modules
module load conda/analysis3-unstable
module load pbs

nosetests --with-xunit test_run.py -s
