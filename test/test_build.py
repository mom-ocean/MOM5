
from model_test_setup import ModelTestSetup

class TestBuild(ModelTestSetup):

    def __init__(self):
        super(TestBuild, self).__init__()

    def test_build_mom_sis(self):

        ret = self.build('MOM_SIS')
        assert(ret == 0)

    def test_build_ebm(self):

        ret = self.build('EBM')
        assert(ret == 0)

    def test_build_CM2(self):

        ret = self.build('CM2M')
        assert(ret == 0)

    def test_build_ESM2M(self):

        ret = self.build('ESM2M')
        assert(ret == 0)

    def test_build_ICCM(self):

        ret = self.build('ICCM')
        assert(ret == 0)
