#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import re
import argparse
import tempfile
import subprocess as sp

base_url = 'https://climate-cms.nci.org.au/repository/entry/{}/Data+Repository/Other+Data+at+NCI/MOM+Test+Data/'

def get_file_list(verbose=False):

    tmp_f = tempfile.NamedTemporaryFile()

    out = sp.check_output(['wget', '-O', tmp_f.name, base_url.format('show')],
                          stderr=sp.STDOUT)
    if verbose:
        print(out, file=sys.stderr)

    filenames = re.findall('MOM Test Data/(.+?\.tar\.gz)', tmp_f.read())
    tmp_f.close()
    assert(filenames != [])

    return filenames

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs='?', default=None, help="""
                        Name of the file to get.""")
    parser.add_argument('--list', action='store_true', default=False,
                        help="List all available files.")
    parser.add_argument('--verbose', action='store_true', default=False,
                        help='Verbose output')

    args = parser.parse_args()

    if args.list:
        available_files = get_file_list(args.verbose)
        print('\n'.join(available_files))
        return 0

    if args.filename is None:
        parser.print_help()
        return 1

    # Set up file source and destination
    my_dir = os.path.dirname(os.path.realpath(__file__))
    dest_dir = os.path.join(my_dir, 'archives')
    dest = os.path.join(dest_dir, args.filename)
    src = base_url.format('get') + args.filename

    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    if os.path.exists(dest):
        print('Error: destination {} already exists.'.format(dest))
        return 1

    ret = sp.call(['wget', '-P', dest_dir, src])
    if ret != 0:
        print('Error: wget of {} failed. Does it exist?'.format(args.filename),
              file=sys.stderr)
        parser.print_help()
        return ret


if __name__ == '__main__':
    sys.exit(main())
