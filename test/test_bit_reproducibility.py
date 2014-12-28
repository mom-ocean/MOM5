
from __future__ import print_function

import os
import sys
import re

from model_test_setup import ModelTestSetup
from test_run import tests as test_specs

class TestBitReproducibility(ModelTestSetup):

    def __init__(self):
        super(TestBitReproducibility, self).__init__()

    def checksums_to_dict(self, filename):
        """
        Look at each line an make a dictionary entry.
        """

        regex = re.compile(r'\[chksum\]\s+(.*)\s+(-?[0-9]+)$')

        dict = {}
        with open(filename) as f:
            for line in f:
                m = regex.match(line)
                if m is not None:
                    dict[m.group(1).rstrip()] = int(m.group(2))

        return dict


    def expected_checksums(self, test_name):

        filename = os.path.join(self.my_dir, 'checksums',
                                '{}.txt'.format(test_name))
        return self.checksums_to_dict(filename)


    def produced_checksums(self, test_name):
        """
        Extract checksums from model run output.
        """

        filename = os.path.join(self.work_dir, test_name, 'fms.out')
        return self.checksums_to_dict(filename)


    def check_run(self, key):

        # Compare expected to produced.
        expected = self.expected_checksums(key)
        produced = self.produced_checksums(key)

        for k in expected:
            assert(produced.has_key(k))
            if expected[k] != produced[k]:
                print('{}: expected {}, produced {}'.format(key, expected[k],
                                                            produced[k]))
            assert(expected[k] == produced[k])

    def test_checksums(self):
        for k in test_specs.keys():
            yield self.check_run, k
