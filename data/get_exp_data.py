#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import re
import argparse
import tempfile
import subprocess as sp

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs='?', default=None, help="""
                        Name of the file to get.""")
    parser.add_argument('--list', action='store_true', default=False,
                        help="List all available files.")
    parser.add_argument('--sources', default=None,
                        help='File specifying the source of each data file.')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help='Verbose output')

    args = parser.parse_args()

    src_dict = {}

    if args.sources is None:
        dir_path = os.path.dirname(os.path.realpath(__file__))
        args.sources = os.path.join(dir_path, 'data_sources.txt')

    # Read in the sources and convert to a dictionary.
    with open(args.sources) as sf:
        for line in sf:
            filename = (line.split(',')[0]).strip()
            urls = [l.strip() for l in line.split(',')[1:]]
            src_dict[filename] = urls

    if args.list:
        print('\n'.join(src_dict.keys()))
        return 0

    if args.filename is None:
        parser.print_help()
        return 1

    # Set up file source and destination
    my_dir = os.path.dirname(os.path.realpath(__file__))
    dest_dir = os.path.join(my_dir, 'archives')
    dest = os.path.join(dest_dir, args.filename)

    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    if os.path.exists(dest):
        print('Error: destination {} already exists.'.format(dest))
        return 1

    # Possibly try multiple sources to get file. 
    downloaded = False
    for url in src_dict[args.filename]:
        ret = sp.call(['wget', '--quiet', '-O',
                        os.path.join(dest_dir, args.filename), url])
        if ret == 0:
            downloaded = True
            break

    if not downloaded:
        print('Error: wget of {} failed. Does it exist?'.format(args.filename),
              file=sys.stderr)
        parser.print_help()
        return ret


if __name__ == '__main__':
    sys.exit(main())
