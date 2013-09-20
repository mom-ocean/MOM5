#!/usr/bin/env python

import sys
import nose
import nose.loader as loader
import optparse

def main():

    parser = optparse.OptionParser()
    parser.add_option("--platform", dest="platform")
    parser.add_option("--type", dest="type")
    parser.add_option("--experiment", dest="experiment")

    (opts, _) = parser.parse_args()

    l = loader.TestLoader()
    suite = l.loadTestsFromName("test_bit_reproducibility.py:TestBitReproducibility.test_%s" % opts.experiment)

    return nose.run(suite=suite, argv=[sys.argv[0]]) 

if __name__ == "__main__":
    sys.exit(main())
