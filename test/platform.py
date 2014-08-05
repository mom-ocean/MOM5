run_scripts = {}
run_scripts['nci'] = """
#!/bin/csh -f

#PBS -P x77
#PBS -q normal
#PBS -l walltime={walltime}
#PBS -l ncpus={ncpus}
#PBS -l mem={mem}
#PBS -l wd
#PBS -o {stdout_file}
#PBS -e {stderr_file}
#PBS -N {run_name}
#PBS -W block=true

./MOM_run.csh --platform nci --type {type} --experiment {exp}
"""

build_cmd = " ./MOM_compile.csh --platform {platform} --type {type}"

local_data_repos = {}
local_data_repos['nci'] = "/short/v45/mom/"
