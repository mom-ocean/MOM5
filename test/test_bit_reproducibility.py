
import os
import sys
import subprocess
import time

# qsub doesn't like the -l wd option on the command line, so need to make a script.

run_script = """
#!/bin/csh -f

#PBS -P v45
#PBS -q express
#PBS -l walltime=00:10:00
#PBS -l ncpus=32
#PBS -l mem=63Gb
#PBS -l wd

./MOM_run.csh --platform nci --type %s --experiment %s --download_input_data
"""

class ModelTestSetup:

    def __init__(self): 

        self.exp_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../', 'exp')

    def build(self, model_type):

        os.chdir(self.exp_path)
        ret = subprocess.check_call(['./MOM_compile.csh', '--platform', 'nci', '--type', model_type])

        if ret == 0:
            return True
        else:
            return False

    def run(self, model_type, exp):

        os.chdir(self.exp_path)
        
        # Write script out as a file.
        with open('run_script.sh', 'w+') as f:
            f.write(run_script % (model_type, exp))

        # Submit the experiment
        run_id = subprocess.check_output(['qsub', 'run_script.sh'])
        run_id = run_id.rstrip()

        # Wait for termination.
        self.wait(run_id)

        # Read the output file and check that run suceeded.
        output = 'run_script.sh.o%s' % run_id.split('.')[0]
        s = ''
        with open(output, 'r') as f:
            s = f.read()
        assert 'NOTE: Natural end-of-script.' in s

        return s

    def wait(self, run_id):

        while True:
            time.sleep(1)
            try: 
                qsub_out = subprocess.check_output(['qstat', run_id])
            except subprocess.CalledProcessError:
                break

            if 'Job has finished' in qsub_out:
                break
            

class TestBitReproducibility(ModelTestSetup):

    def get_norms(self, output):
        """
        Extract the norms from model run output.
        """

        return None

    def test_MOM_SIS(self):
        """
        Test whether MOM_SIS model runs give reproducible and expected results.

        Run the model for a single day (or month, for now), read a norm from the output, compare to existing norm. 
        """

        expected_norms = []

        assert self.build('MOM_SIS')

        output = self.run('MOM_SIS', 'om3_core3')

        assert get_norms(output) == expected_norms, "Norms do not match."

