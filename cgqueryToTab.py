#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-10-06 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
cgqueryToTab.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: cgquery's output format is horrendous. this reformats it")
    # parser.add_argument('-a', '--argA', metavar='argA', type=str, required=True, help='description of argument')
    # parser.add_argument('-b', '--argB', metavar='argB', required=False, help='description of argument B')
    # parser.add_argument('-c', '--flagC', required=False, action='store_true', help='sets flagC to true')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='cgquery data to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    # send back the user input
    return args

# primary function
def parse_cgquery(file):
    col_keys = ['analysis_id', 'state', 'last_modified', 'upload_date', 'published_date', 'center_name', 'study', 'aliquot_id', 'sample_accession', 'legacy_sample_id', 'disease_abbr', 'tss_id', 'participant_id', 'sample_id', 'analyte_code', 'sample_type', 'library_strategy', 'platform', 'refassem_short_name', 'analysis_submission_uri', 'analysis_full_uri', 'analysis_data_uri']
    print '#' + '\t'.join(['%s_%s' % (i+1, col_keys[i]) for i in xrange(len(col_keys))])
    analysis = {} # dict
    for line in file:
        v = line.strip().split(' : ')
        if v[0][:9]=='Analysis ':
            write_out(analysis, col_keys)
            analysis = {} # dict
            continue
        for i in xrange(len(v)):
            v[i] = v[i].strip()
            if len(v) == 2 and v[0] in col_keys:
                analysis[v[0]] = v[1]
    return

def write_out(analysis, col_keys):
    if analysis:
        for key in col_keys:
            if key not in analysis:
                analysis[key] = ''
        print '\t'.join([analysis[x] for x in col_keys])
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    parse_cgquery(args.input)

    # close the input file
    args.input.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
