
from __future__ import print_function

from model_test_setup import ModelTestSetup

class TestRun(ModelTestSetup):

    # Run tests in parallel.
    # Run with nosetests test_run.py --processes=<n>
    _multiprocess_can_split_ = True

    def __init__(self):
        super(TestRun, self).__init__()

    def test_om3_core3(self):

        r, so, se = self.run('MOM_SIS', 'om3_core3', npes=24)
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_om3_core1(self):

        r, so, se = self.run('MOM_SIS', 'om3_core1', npes=24)
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_atlantic1(self):

        r, so, se = self.run('MOM_SIS', 'atlantic1', ncpus='32', npes='24', mem='64Gb')
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_mom4p1_ebm1(self):

        r, so, se = self.run('EBM', 'mom4p1_ebm1', ncpus='32', npes='17', mem='64Gb')
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_MOM_SIS_TOPAZ(self):

        r, so, se = self.run('MOM_SIS', 'MOM_SIS_TOPAZ', npes='24')
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_MOM_SIS_BLING(self):

        r, so, se = self.run('MOM_SIS', 'MOM_SIS_BLING', npes=24)
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_CM2_1p1(self):

        r, so, se = self.run('CM2M', 'CM2.1p1', ncpus='64', npes='45', mem='128Gb')
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_CM2M_coarse_BLING(self):

        r, so, se = self.run('CM2M', 'CM2M_coarse_BLING', ncpus='64',
                             npes='45', mem='128Gb')
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_ICCMp1(self):

        r, so, se = self.run('ICCM', 'ICCMp1', ncpus='64', npes='45', mem='128Gb')
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

    def test_ESM2M_pi_control_C2(self):

        r, so, se = self.run('ESM2M', 'ESM2M_pi-control_C2', ncpus='128',
                             npes='120', mem='256Gb')
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

