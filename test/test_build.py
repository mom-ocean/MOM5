
from model_test_setup import ModelTestSetup

class TestBuild(ModelTestSetup):

    def __init__(self):
        super(TestBuild, self).__init__()

    def test_build_mom_sis(self):

        self.build('MOM_SIS')

    def test_build_ebm(self):

        self.build('EBM')

    def test_build_CM2(self):

        self.build('CM2M')

    def test_build_ESM2M(self):

        self.build('ESM2M')

    def test_build_ICCM(self):

        self.build('ICCM')
