
from __future__ import print_function

import os
import sys
import shutil

from model_test_setup import ModelTestSetup
from test_run import tests as test_specs

class TestValgrind(ModelTestSetup):
    """
    Run all test cases under valgrind.
    """

    # Run tests in parallel.
    # Run with nosetests test_valgrind.py --processes=<n>
    _multiprocess_can_split_ = True

    def __init__(self):
        super(TestValgrind, self).__init__()

    def check_run(self, key, qsub=True):

        args = test_specs[key][0]
        kwargs = test_specs[key][1]
        kwargs['qsub'] = qsub
        kwargs['walltime'] = '05:00:00'
        kwargs['valgrind'] = True

        print('####### Running {}.{} in Valgrind ########'.format(args[0], args[1]))
        # Clean out the work directory.
        if os.path.exists(os.path.join(self.work_dir, key)):
            shutil.rmtree(os.path.join(self.work_dir, key))
        r, so, se = self.run(*args, **kwargs)
        print(so)
        print(se)
        assert(r == 0)
        assert('NOTE: Natural end-of-script.' in so)

        # FIXME: check that valgrind doesn't produce any errors.

    def test_experiments(self):
        for k in test_specs.keys():
            yield self.check_run, k
