module purge

module use /g/data3/hh5/public/modules
module load conda/analysis3-unstable
module load pbs

nosetests --with-xunit test_run.py -s
