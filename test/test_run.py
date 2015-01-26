
from __future__ import print_function

from model_test_setup import ModelTestSetup

import os
import shutil

# This defines the different tests. To run an individual test on the command
# line type (for example):
# $ python -c "import test_run ; tc = test_run.TestRun() ; test_run.TestRun.check_run(tc, 'MOM_SIS.om3_core3')"
#
# If you are in an interactive qsub job (i.e. do not want the test harness to
# submit a job, then type:
# $ python -c "import test_run ; tc = test_run.TestRun() ; test_run.TestRun.check_run(tc, 'om3_core3', qsub=False)"

tests = {'om3_core3' : (('MOM_SIS', 'om3_core3'), {'ncpus' : '32', 'npes' : '24'}),
         'om3_core1' : (('MOM_SIS', 'om3_core1'), {'ncpus' : '32', 'npes' : '24'}),
         'atlantic1' : (('MOM_SIS', 'atlantic1'), {'ncpus' : '32', 'npes' : '24', 'mem' : '64Gb'}),
         'mom4p1_ebm1' : (('EBM', 'mom4p1_ebm1'), {'ncpus' : '32', 'npes' : '17', 'mem' : '64Gb'}),
         'MOM_SIS_TOPAZ' : (('MOM_SIS', 'MOM_SIS_TOPAZ'), {'ncpus' : '32', 'npes' : '24'}),
         'MOM_SIS_BLING' : (('MOM_SIS', 'MOM_SIS_BLING'), {'ncpus' : '32', 'npes' : '24'}),
         'CM2.1p1' : (('CM2M', 'CM2.1p1'), {'ncpus' : '64', 'npes' : '45', 'mem' : '128Gb'}),
         'CM2M_coarse_BLING' : (('CM2M', 'CM2M_coarse_BLING'), {'ncpus' : '64', 'npes' : '45', 'mem' : '128Gb'}),
         'ICCMp1' :  (('ICCM', 'ICCMp1'), {'ncpus' : '64', 'npes' : '54', 'mem' : '128Gb'}),
         'ESM2M_pi-control_C2' : (('ESM2M', 'ESM2M_pi-control_C2'), {'ncpus' : '128', 'npes' : '90', 'mem' : '256Gb'}),
         'global_0.25_degree_NYF' : (('MOM_SIS', 'global_0.25_degree_NYF'), {'ncpus' : '960', 'npes' : '960', 'mem' : '1900Gb'})}

class TestRun(ModelTestSetup):
    """
    Run all test cases and check for successful output.
    """

    # Run tests in parallel.
    # Run with nosetests test_run.py --processes=<n>
    _multiprocess_can_split_ = True

    def __init__(self):
        super(TestRun, self).__init__()

    def check_run(self, key, qsub=True):

        args = tests[key][0]
        kwargs = tests[key][1]
        kwargs['qsub'] = qsub

        print('############ Running {}.{} ############'.format(args[0], args[1]))
        # Clean out the work directory.
        if os.path.exists(os.path.join(self.work_dir, key)):
            shutil.rmtree(os.path.join(self.work_dir, key))
        r, so, se = self.run(*args, **kwargs)
        if r != 0:
            print(so)
            print(se)
            print('Run exited with error {}'.format(r))
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_experiments(self):
        for k in tests.keys():
            yield self.check_run, k
