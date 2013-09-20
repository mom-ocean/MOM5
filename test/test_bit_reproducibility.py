
import os
import sys
import subprocess
import time

# qsub doesn't like the -l wd option on the command line, so need to make a script.

run_script = """
#!/bin/csh -f

#PBS -P v45
#PBS -q %s
#PBS -l walltime=%s
#PBS -l ncpus=%s
#PBS -l mem=%s
#PBS -l wd
#PBS -N %s

./MOM_run.csh --platform nci --type %s --experiment %s --download_input_data --npes %s
"""

class ModelTestSetup(object):

    def __init__(self): 

        self.my_path = os.path.dirname(os.path.realpath(__file__))
        self.exp_path = os.path.join(self.my_path, '../', 'exp')

    def run(self, model_type, exp, queue='normal', walltime='00:10:00', ncpus='16', mem='32Gb',):

        os.chdir(self.exp_path)

        run_name = "TC_%s" % exp
        
        # Write script out as a file.
        with open('run_script.sh', 'w+') as f:
            f.write(run_script % (queue, walltime, ncpus, mem, run_name, model_type, exp, ncpus))

        # Submit the experiment
        run_id = subprocess.check_output(['qsub', 'run_script.sh'])
        run_id = run_id.rstrip()

        # Wait for termination.
        self.wait(run_id)

        # Read the output file and check that run suceeded.
        output = '%s.o%s' % (run_name, run_id.split('.')[0])
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

        output = self.run(type, experiment)

        assert self.get_checksums(output) == self.expected_checksums(type, experiment), "Checksums do not match."

    def test_om3_core1(self):

        self.run('MOM_SIS', 'om3_core1')

    def test_atlantic1(self):

        self.run('MOM_SIS', 'atlantic1', ncpus='32', mem='64Gb')

    def test_MOM_SIS_TOPAZ(self):

        self.run('MOM_SIS', 'MOM_SIS_TOPAZ')

    def test_MOM_SIS_BLING(self):

        self.run('MOM_SIS', 'MOM_SIS_BLING')

    def test_CM2_1_p1(self):

        self.run('CM2M', 'CM2.1.P1', ncpus='64', mem='128Gb')

    def test_CM2M_coarse_BLING(self):

        self.run('CM2M', 'CM2M_coarse_BLING', ncpus='64', mem='128Gb')

    def test_ESM2M_pi_control_C2(self):

        self.run('ESM2M', 'ESM2M_pi_control_C2', ncpus='128', mem='256Gb')

    def test_ICCMp1(self):

        self.run('ICCM', 'ICCMp1', ncpus='64', mem='128Gb')


