
from model_test_setup import ModelTestSetup

class TestBuild(ModelTestSetup):

    def __init__(self):
        super(TestBuild, self).__init__()

    def test_MOM_SIS(self):

        ret = self.build('MOM_SIS')
        assert(ret == 0)

    def test_EBM(self):

        ret = self.build('EBM')
        assert(ret == 0)

    def test_CM2M(self):

        ret = self.build('CM2M')
        assert(ret == 0)

    def test_ESM2M(self):

        ret = self.build('ESM2M')
        assert(ret == 0)

    def test_ICCM(self):

        ret = self.build('ICCM')
        assert(ret == 0)
