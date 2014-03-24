
from model_test_setup import ModelTestSetup

class TestRun(ModelTestSetup):

    def __init__(self):
        super(TestRun, self).__init__()

    def test_om3_core3(self):

        self.run('MOM_SIS', 'om3_core3')

    def test_om3_core1(self):

        self.run('MOM_SIS', 'om3_core1')

    def test_atlantic1(self):

        self.run('MOM_SIS', 'atlantic1', ncpus='32', npes='24', mem='64Gb')

    def test_mom4p1_ebm1(self):

        self.run('EBM', 'mom4p1_ebm1', ncpus='32', npes='17', mem='64Gb')

    def test_MOM_SIS_TOPAZ(self):

        self.run('MOM_SIS', 'MOM_SIS_TOPAZ')

    def test_MOM_SIS_BLING(self):

        self.run('MOM_SIS', 'MOM_SIS_BLING')

    def test_CM2_1p1(self):

        self.run('CM2M', 'CM2.1p1', ncpus='64', npes='45', mem='128Gb')

    def test_CM2M_coarse_BLING(self):

        self.run('CM2M', 'CM2M_coarse_BLING', ncpus='64', npes='45', mem='128Gb')

    def test_ESM2M_pi_control_C2(self):

        self.run('ESM2M', 'ESM2M_pi-control_C2', ncpus='128', npes='120', mem='256Gb')

    def test_ICCMp1(self):

        self.run('ICCM', 'ICCMp1', ncpus='64', npes='45', mem='128Gb')
