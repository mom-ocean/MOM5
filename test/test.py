#!/usr/bin/env python

import sys
import os
import nose
import nose.loader as loader
import subprocess
import test_bit_reproducibility as tb
import optparse

def main():

    parser = optparse.OptionParser()
    parser.add_option("--platform")
    parser.add_option("--type")
    parser.add_option("--experiment")

    (opts, _) = parser.parse_args()

    experiment = opts.experiment.replace('.', '_').replace('-', '_')

    if 'test_%s' % experiment in dir(tb.TestBitReproducibility):
        l = loader.TestLoader()
        suite = l.loadTestsFromName("test_bit_reproducibility.py:TestBitReproducibility.test_%s" % experiment)
        # nose returns False on failure - want 1. 
        ret = not nose.run(suite=suite, argv=[sys.argv[0], '-s']) 
    else:

        my_path = os.path.dirname(os.path.realpath(__file__))
        exp_path = os.path.join(my_path, '../', 'exp')
        os.chdir(exp_path)

        # Specific test was not found. Try to just run the experiment directly.
        ret = subprocess.check_call(['./MOM_run.csh', '--platform', opts.platform, '--type', opts.type, '--experiment', opts.experiment, '--download_input_data'])

    return ret

if __name__ == "__main__":
    sys.exit(main())
