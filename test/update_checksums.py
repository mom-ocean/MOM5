#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import re

from test_run import tests as test_specs

"""
This script updates the files in checksums/ based on the latest run.
It should only be used when the checksum changes are expected.
"""

def main():

    my_dir = os.path.dirname(os.path.realpath(__file__))
    work_dir = os.path.join(my_dir, '../', 'work')
    checksum_dir = os.path.join(my_dir, 'checksums')

    regex = re.compile(r'\[chksum\]\s+(.*)\s+(-?[0-9]+)$')

    for test_name in test_specs.keys():
        model_out = os.path.join(work_dir, test_name, 'fms.out')
        checksum_file = os.path.join(checksum_dir, test_name + '.txt')

        with open(model_out, 'r') as m_f:
            with open(checksum_file, 'w') as c_f:
                for line in m_f:
                    m = regex.match(line)
                    if m is not None:
                        print(m.group(0), file=c_f)

if __name__ == '__main__':
    sys.exit(main())
