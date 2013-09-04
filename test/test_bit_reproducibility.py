
import os
import subprocess
import time

class ModelTestSetup:

    def __init__(self): 
        pass

    def build(self, model_type):

        os.chdir(os.path.join(os.path.realpath(__file__), '../', 'exp'))
        return subprocess.check_call(['./MOM_compile.csh', '--platform', 'nci', '--type', model_type])

    def run(self, model_type, exp):

        os.chdir(os.path.join(os.path.realpath(__file__), '../', 'exp'))
        
        # Submit the experiment
        qsub_args = ['qsub', '-P', 'v45', '-q', 'normal', '-k', 'oe', '-l', 'walltime=00:10:00,ncpus=16,mem=31Gb', '-wd']
        command_args = ['./MOM_run.csh', '--platform', 'nci', '--type', model_type, '--experiment', exp]
        run_id = subprocess.check_output(qsub_args + ['--'] + command_args)

        # Wait for termination.
        self.wait(run_id)

        # Read the output file. 

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

    def test_MOM_SIS(self):
        """
        Test whether MOM_SIS model runs give reproducible and expected results.

        Run the model for a single day (or month, for now), read a norm from the output, compare to existing norm. 
        """

        expected_norms = [
                ]
        assert self.build('MOM_SIS')

        output = self.run('MOM_SIS', 'om3_core3')

        assert get_norms(output) == expected_norms, "Norms do not match."

