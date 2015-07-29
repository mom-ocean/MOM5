
from model_test_setup import ModelTestSetup

"""
To run a single test, e.g. MOM_SIS:

python -c "import test_build ; tc = test_build.TestBuild() ; test_build.TestBuild.check_build(tc, 'MOM_SIS')"
"""

test_args = ['MOM_solo', 'MOM_SIS', 'EBM', 'CM2M', 'ESM2M', 'ICCM',
             'ACCESS-CM', 'ACCESS-OM']

class TestBuild(ModelTestSetup):
    """
    Build all model types.
    """

    def __init__(self):
        super(TestBuild, self).__init__()

    def check_build(self, model):
        ret = self.build(model)
        assert(ret == 0)

    def test_builds(self):
        for t in test_args:
            yield self.check_build, t
