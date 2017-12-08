#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import re
import argparse
import tempfile
import subprocess as sp
import shutil

def get_local_path(url):
    """
    Convert the url to a local file path if it exists.

    Return whether or not this is a local path and it's value.
    """

    prefix = url[:7]
    path = url[7:]

    if prefix == 'file://':
        if os.path.exists(path):
            return True, path
        else:
            return True, None

    return False, None

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
        args.sources = os.path.join(dir_path, 'data_sources.csv')

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
    ret = 0
    for url in src_dict[args.filename]:
        is_local_path, local_path = get_local_path(url)

        if local_path is not None:
            shutil.copy(local_path, dest_dir)
            break
        elif not is_local_path:
            ret = sp.call(['wget', '--quiet', '-O',
                            os.path.join(dest_dir, args.filename), url])
            if ret == 0:
                break

    if not os.path.exists(os.path.join(dest_dir, args.filename)) or ret != 0:
        print('Error: wget of {} failed. Does it exist?'.format(args.filename),
              file=sys.stderr)
        parser.print_help()
        return ret


if __name__ == '__main__':
    sys.exit(main())
