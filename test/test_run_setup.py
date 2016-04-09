
from __future__ import print_function

import os
import shutil

from model_test_setup import ModelTestSetup
from test_run import tests as experiments

class TestRunSetup(ModelTestSetup):
    """
    Do setup to run all experiments.

    This is kep separate from test_run.py because it can be
    run on a single core (i.e. no need for qsub)
    """

    def __init__(self):
        super(TestRunSetup, self).__init__()

    def check_build(self, model):
        ret = self.build(model, unit_testing=False)
        assert(ret == 0)

    def check_download_data(self, experiment):
        ret = self.download_input_data(experiment)
        assert(ret == 0)

    def test_setup(self):

        # Clean out the work directory.
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)

        # Download data
        for e in experiments.keys():
            yield self.check_download_data, e

        # Do build
        models = list(set([v[0][0] for v in experiments.values()]))
        for m in models:
            print(m)
            yield self.check_build, m

