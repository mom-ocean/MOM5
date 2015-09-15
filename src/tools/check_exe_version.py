#!/usr/bin/env python

"""
This script checks whether a MOM executable was build with the HEAD of a
particular repo. This is useful 

It returns 0 if they are the same, 
"""

from __future__ import print_function

import sys
import argparse
import os
import subprocess as sp
import shlex
import re


def main():

    # Path to this script
    my_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser()
    parser.add_argument("exe", help="""Path the MOM executable.""")
    parser.add_argument("--git_repo", default=my_path, help="""
                        The path to the git repo.""")
    parser.add_argument("--verbose", default=False, action='store_true',
                        help="""Verbose output, will print the hashes.""")
    args = parser.parse_args()

    readelf_out = sp.check_output(shlex.split('readelf -p .rodata {}'.format(args.exe)))
    m = re.search('MOM_VERSION=(\w{40})', readelf_out)
    exe_hash = m.group(1)

    if args.verbose:
        print('Exe hash {}'.format(exe_hash), file=sys.stderr)

    curr_dir = os.getcwd()
    os.chdir(args.git_repo)
    git_out = sp.check_output(shlex.split('git rev-parse HEAD'))
    os.chdir(curr_dir)

    if args.verbose:
        print('Repo hash {}'.format(git_out.strip()), file=sys.stderr)

    if exe_hash == git_out.strip():
        return 0
    else:
        return 1

if __name__ == '__main__':
    sys.exit(main())
