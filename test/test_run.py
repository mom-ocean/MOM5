
from __future__ import print_function

from model_test_setup import ModelTestSetup

test_args = [(('MOM_SIS', 'om3_core3'), {'ncpus' : '32'}),
             (('MOM_SIS', 'om3_core1'), {'ncpus' : '32'}),
             (('MOM_SIS', 'atlantic1'),
                         {'ncpus' : '32', 'npes' : '24', 'mem' : '64Gb'}),
             (('EBM', 'mom4p1_ebm1'),
                     {'ncpus' : '32', 'npes' : '17', 'mem' : '64Gb'}),
             (('MOM_SIS', 'MOM_SIS_TOPAZ'), {'ncpus' : '24'}),
             (('MOM_SIS', 'MOM_SIS_BLING'), {'ncpus' : '24'}),
             (('CM2M', 'CM2.1p1'),
                      {'ncpus' : '64', 'npes' : '45', 'mem' : '128Gb'}),
             (('CM2M', 'CM2M_coarse_BLING'),
                     {'ncpus' : '64', 'npes' : '45', 'mem' : '128Gb'}),
             (('ICCM', 'ICCMp1'), {'ncpus' : '64', 'npes' : '54', 'mem' : '128Gb'}),
             (('ESM2M', 'ESM2M_pi-control_C2'),
                      {'ncpus' : '128', 'npes' : '120', 'mem' : '256Gb'})]


class TestRun(ModelTestSetup):
    """
    Run all test cases and check for successful output.

    FIXME: use a test generator.
    """

    # Run tests in parallel.
    # Run with nosetests test_run.py --processes=<n>
    _multiprocess_can_split_ = True

    def __init__(self):
        super(TestRun, self).__init__()

    def check_run(self, args, kwargs):

        print('############ Running {}.{} ############'.format(args[0], args[1]))
        r, so, se = self.run(*args, **kwargs)
        if r != 0:
            print(so)
            print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_experiments(self):
        for t in test_args: 
            yield self.check_run, t[0], t[1]

