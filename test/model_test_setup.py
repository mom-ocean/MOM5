
import os
import sys
import subprocess
import time

class ModelTestSetup(object):

    def __init__(self): 

        self.my_path = os.path.dirname(os.path.realpath(__file__))
        self.exp_path = os.path.join(self.my_path, '../', 'exp')

    def run(self, model_type, exp, walltime='00:10:00', ncpus='32', npes=None, mem='64Gb'):
        """
        ncpus is for requested cpus, npes is for how many mom uses.
        """

        os.chdir(self.exp_path)

        run_name = "TC_%s" % exp
        # -N value is a maximum of 15 chars.
        run_name = run_name[0:15]

        if npes != None:
            npes = '--npes %s' % npes
        else:
            npes = ''
        
        # Write script out as a file.
        with open('run_script.sh', 'w+') as f:
            if npes != None:
                f.write(run_script % (walltime, ncpus, mem, run_name, model_type, exp, npes))

        # Submit the experiment
        run_id = subprocess.check_output(['qsub', 'run_script.sh'])
        run_id = run_id.rstrip()

        # Wait for termination.
        self.wait(run_id)

        # Read the output file and check that run suceeded.
        output = '%s.o%s' % (run_name, run_id.split('.')[0])
        error = '%s.e%s' % (run_name, run_id.split('.')[0])
        so = ''
        with open(output, 'r') as f:
            so = f.read()
        with open(error, 'r') as f:
            se = f.read()

        os.chdir(self.my_path)

        print so
        print se
        assert 'NOTE: Natural end-of-script.' in so

        return (so, se)

    def wait(self, run_id):

        while True:
            time.sleep(10)
            qsub_out = ''
            try: 
                qsub_out = subprocess.check_output(['qstat', run_id], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as err:
                qsub_out = err.output

            if 'Job has finished' in qsub_out:
                break

