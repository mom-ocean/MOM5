
import os
import sys
import subprocess
import time

# qsub doesn't like the -l wd option on the command line, so need to make a script.
# NOTE: this is only an 8 core job, however needs a lot of memory, hence the oversubscription.

run_script = """
#!/bin/csh -f

#PBS -P v45
#PBS -q express
#PBS -l walltime=00:10:00
#PBS -l ncpus=16
#PBS -l mem=32Gb
#PBS -l wd

./MOM_run.csh --platform nci --type %s --experiment %s --download_input_data
"""

class ModelTestSetup(object):

    def __init__(self): 

        self.my_path = os.path.dirname(os.path.realpath(__file__))
        self.exp_path = os.path.join(self.my_path, '../', 'exp')

    def build(self, model_type):

        os.chdir(self.exp_path)
        ret = subprocess.check_call(['./MOM_compile.csh', '--platform', 'nci', '--type', model_type])
        os.chdir(self.my_path)

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

        os.chdir(self.my_path)

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

    def __init__(self):
        super(TestBitReproducibility, self).__init__()

    def expected_checksums(self, type, experiment):

        s = ''
        with open(os.path.join(self.my_path, '%s.%s.checksums.txt' % (type, experiment))) as f:
            s = f.read()

        return s

    def get_checksums(self, output):
        """
        Extract checksums from model run output.
        """

        s = ''
        for line in output.splitlines(True):
            if '[chksum]' in line:
                s += line

        return s

    def test_om3_core3(self):
        """
        Test whether MOM_SIS model runs give reproducible and expected results.

        Run the model for a single day (or month, for now), read a norm from the output, compare to existing norm. 
        """

        type = 'MOM_SIS'
        experiment = 'om3_core3'

        assert self.build(type)

        output = self.run(type, experiment)

        assert self.get_checksums(output) == self.expected_checksums(type, experiment), "Checksums do not match."
