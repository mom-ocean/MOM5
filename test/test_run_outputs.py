
from __future__ import print_function

from model_test_setup import ModelTestSetup

import os
import sys
import shutil

from test_run import tests as experiments

class TestRunOutputs(ModelTestSetup):
    """
    Check the output of all runs only. Don't actually do the runs.
    """

    # Run tests in parallel.
    # Run with nosetests test_run.py --processes=<n>
    _multiprocess_can_split_ = True

    def __init__(self):
        super(TestRunOutputs, self).__init__()

    def check_run_output(self, key):

        output_filename = os.path.join(self.work_dir, key, 'fms.out')
        assert(os.path.exists(output_filename))

        with open(output_filename) as f:
            s = f.read()

        assert('NOTE: Natural end-of-script for experiment {} with model {}'.format(key, experiments[key][0][0]) in s)

    def test_experiments(self):
        for k in experiments.keys():
            yield self.check_run_output, k
