run_script = """
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

build_cmd = " ./MOM_compile.csh --platform nci --type %s"
