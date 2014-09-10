
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

        # Set the local remote if there is one.
        if plat.local_data_repos.has_key(self.get_platform()):
            remote = plat.local_data_repos[self.get_platform()]
            cmd = '/usr/bin/git remote add local_data {}'.format(remote)
            try: 
                sp.check_output(shlex.split(cmd), stderr=sp.STDOUT)
            except sp.CalledProcessError as err:
                # This is allowed to fail in this case
                assert('remote local_data already exists' in err.output)

            cmd = '/usr/bin/git annex get {} --from local_data'.format(input)
            ret = sp.call(shlex.split(cmd))
        else:
            # Otherwise data will be download from Amazon S3. 
            cmd = '/usr/bin/git annex get {}'.format(input)
            ret = sp.call(shlex.split(cmd))

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
 

    def run(self, model_type, exp, walltime='00:10:00', ncpus='32',
            npes=None, mem='64Gb'):
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
        ret = sp.call(['qsub', run_file])
        stdout, stderr = self.get_output(fo, fe)

        # Clean up temporary files. 
        os.remove(stdout_file)
        os.remove(stderr_file)
        os.remove(run_file)

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
