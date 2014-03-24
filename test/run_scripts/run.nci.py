run_script = """
#!/bin/csh -f

#PBS -P x77
#PBS -q %s
#PBS -l walltime=%s
#PBS -l ncpus=%s
#PBS -l mem=%s
#PBS -l wd
#PBS -N %s

./MOM_run.csh --platform nci --type %s --experiment %s --download_input_data %s
"""

build_script = """

./MOM_build.csh --platform nci --type %s --experiment %s --download_input_data %s
"""
