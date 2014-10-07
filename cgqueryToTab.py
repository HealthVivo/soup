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
    analysis = None
    first_entry = True
    for line in file:
        v = line.strip().split(' : ')
        if v[0][:9]=='Analysis ':
            if first_entry and analysis:
                header = []
                name_index = 1
                for name in analysis:
                    header.append("%s_%s" % (name_index, name))
                    name_index += 1
                print '#' + '\t'.join(header)
                first_entry = False
            write_out(analysis)
            analysis = {} # dict
            continue
        if type(analysis) is dict:
            for i in xrange(len(v)):
                v[i] = v[i].strip()
                if len(v) == 2:
                    analysis[v[0]] = v[1]
    return

def write_out(analysis):
    if analysis:
        print '\t'.join([analysis[x] for x in analysis])
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
