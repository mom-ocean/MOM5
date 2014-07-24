
from __future__ import print_function

import os
import sys
import subprocess as sp
import shlex
import tempfile
import time
import run_scripts.nci as plat

class ModelTestSetup(object):

    def __init__(self): 

        self.my_path = os.path.dirname(os.path.realpath(__file__))
        self.exp_path = os.path.join(self.my_path, '../', 'exp')
        self.archive_path = os.path.join(self.my_path, '../data/archives')

    def download_input_data(self, exp):

        os.chdir(self.archive_path)

        cmd = '/usr/bin/git annex get {}.input.tar.gz'.format(exp)
        ret = sp.call(shlex.split(cmd))

        os.chdir(self.my_path)

        return ret

    def get_output(self, fo, fe):

        # The command has finished. Read output and write stdout.
        # We don't know when output has stopped so just keep trying
        # until it is all gone. 
        empty_reads = 0
        stderr = ''
        stdout = ''
        while True:
            so = os.read(fo, 1024*1024)
            se = os.read(fe, 1024*1024)

            if so == '' and se == '':
                empty_reads += 1
            else:
                stdout += so
                stderr += se
                empty_reads = 0

            if empty_reads > 10:
                break

            time.sleep(2)

        return (stdout, stderr)


    def run(self, model_type, exp, walltime='00:10:00', ncpus='32',
            npes=None, mem='64Gb'):
        """
        ncpus is for requested cpus, npes is for how many mom uses.
        """

        ret = self.download_input_data(exp)
        if ret != 0:
            print('Error: could not download input data.', file=sys.stderr)
            return (ret, None, None)

        os.chdir(self.exp_path)

        run_name = "CI_%s" % exp
        # -N value is a maximum of 15 chars.
        run_name = run_name[0:15]

        if npes != None:
            npes = '--npes %s' % npes
        else:
            npes = ''

        # Get temporary file names for the stdout, stderr.
        fo, stdout_file = tempfile.mkstemp(dir=self.exp_path)
        fe, stderr_file = tempfile.mkstemp(dir=self.exp_path)
        
        # Write script out as a file.
        run_script = plat.run_script.format(walltime=walltime, ncpus=ncpus,
                                            mem=mem, stdout_file=stdout_file,
                                            stderr_file=stderr_file,
                                            run_name=run_name,
                                            type=model_type, exp=exp)
        # Write out run script
        frun, run_file = tempfile.mkstemp(dir=self.exp_path)
        os.write(frun, run_script)
        os.close(frun)
        os.chmod(run_file, 0755)

        # Submit the experiment. This will block until it has finished.
        ret = sp.call(['qsub', run_file])
        stdout, stderr = self.get_output(fo, fe)

        # Clean up temporary files. 
        os.remove(stdout_file)
        os.remove(stderr_file)
        os.remove(run_script)

        # Change back to test dir. 
        os.chdir(self.my_path)

        return (ret, stdout, stderr)


    def build(self, model_type):

        os.chdir(self.exp_path)

        build_cmd = plat.build_cmd.format(model_type)
        # Build the model.
        ret = sp.call(shlex.split(build_cmd))

        os.chdir(self.my_path)

        return ret
