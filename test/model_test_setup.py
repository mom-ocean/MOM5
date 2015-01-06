
from __future__ import print_function

import os
import sys
import subprocess as sp
import shlex
import shutil
import tempfile
import time
import platform as plat

class ModelTestSetup(object):

    def __init__(self): 

        self.my_dir = os.path.dirname(os.path.realpath(__file__))
        self.exp_dir = os.path.join(self.my_dir, '../', 'exp')
        self.archive_dir = os.path.join(self.my_dir, '../data/archives')
        self.work_dir = os.path.join(self.my_dir, '../', 'work')

    def download_input_data(self, exp):
        """
        Download the experiment input data.

        This needs to be done before submitting the MOM_run.sh script because
        the compute nodes may not have Internet access.
        """

        os.chdir(self.archive_dir)

        input = '{}.input.tar.gz'.format(exp)

        ret = 0
        if not os.path.exists(input):
            cmd = '/usr/bin/git annex get {}'.format(input)
            ret = sp.call(shlex.split(cmd))

        if not os.path.exists(input):
            cmd = 'wget ftp.gfdl.noaa.gov:/perm/MOM4/mom5_pubrel_dec2013/exp/{}'.format(input)
            ret = sp.call(shlex.split(cmd))

        if ret != 0:
            return ret

        # Unzip into work directory.
        if not os.path.exists(self.work_dir):
            os.mkdir(self.work_dir)

        if not os.path.exists(os.path.join(self.work_dir, input)):
            shutil.copy(input, self.work_dir)

        if not os.path.exists(os.path.join(self.work_dir, exp)):
            os.chdir(self.work_dir)
            cmd = '/bin/tar -xvf {}'.format(input)
            ret += sp.call(shlex.split(cmd))

        os.chdir(self.my_dir)

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


    def get_platform(self):

        # We need to get the node/platform - see if Jenkins has this set.
        platform = 'nci'

        try:
            platform = os.environ['label']
        except KeyError:
            pass

        return platform


    def run(self, model_type, exp, walltime='00:30:00', ncpus='32',
            npes=None, mem='64Gb', qsub=True):
        """
        ncpus is for requested cpus, npes is for how many mom uses.
        """

        ret = self.download_input_data(exp)
        if ret != 0:
            print('Error: could not download input data.', file=sys.stderr)
            return (ret, None, None)

        os.chdir(self.exp_dir)

        run_name = "CI_%s" % exp
        # -N value is a maximum of 15 chars.
        run_name = run_name[0:15]

        if npes != None:
            npes = '--npes %s' % npes
        else:
            npes = ''

        # Get temporary file names for the stdout, stderr.
        fo, stdout_file = tempfile.mkstemp(dir=self.exp_dir)
        fe, stderr_file = tempfile.mkstemp(dir=self.exp_dir)

        # Write script out as a file.
        run_script = plat.run_scripts[self.get_platform()]
        run_script = run_script.format(walltime=walltime, ncpus=ncpus,
                                       mem=mem, stdout_file=stdout_file,
                                       stderr_file=stderr_file,
                                       run_name=run_name,
                                       type=model_type, exp=exp, npes=npes)

        # Write out run script
        frun, run_file = tempfile.mkstemp(dir=self.exp_dir)
        os.write(frun, run_script)
        os.close(frun)
        os.chmod(run_file, 0755)

        # Submit the experiment. This will block until it has finished.
        if qsub:
            ret = sp.call(['qsub', run_file])
        else:
            ret = sp.call([run_file])

        stdout, stderr = self.get_output(fo, fe)

        # Move temporary files to experiment directory.
        shutil.move(stdout_file, os.path.join(self.work_dir, exp, 'fms.out'))
        shutil.move(stderr_file, os.path.join(self.work_dir, exp, 'fms.err'))
        shutil.move(run_file, os.path.join(self.work_dir, exp, 'run.sh'))

        # Change back to test dir.
        os.chdir(self.my_dir)

        return (ret, stdout, stderr)


    def build(self, model_type):

        os.chdir(self.exp_dir)

        platform = self.get_platform()
        build_cmd = plat.build_cmd.format(type=model_type, platform=platform)
        # Build the model.
        ret = sp.call(shlex.split(build_cmd))

        os.chdir(self.my_dir)

        return ret
