#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam
import sys
import getopt
import string
from string import *
from optparse import OptionParser


def bamtofastq(bamfile):
# get file and header
#   bam = pysam.Samfile(bamfile, "rb")
    if bamfile == "stdin": 
        bam = pysam.Samfile("-", "r")
    else:
        bam = pysam.Samfile(bamfile, "rb")
    num = 0
    d = {}
    for al in bam.fetch():
        if al.is_secondary: continue
        key = al.qname
        if key not in d:
            d.setdefault(key,al)
        else:
            for tag, value in d[key].tags:
                if tag=="RG": RG1=value
            for tag, value in al.tags:
                if tag=="RG": RG2=value
            # RG:Z:ID
            if al.is_read1:
                printfastq_rg(al,1,RG2)
                printfastq_rg(d[key],2,RG1)
            else:
                printfastq_rg(d[key],1,RG1)
                printfastq_rg(al,2,RG2)
            del d[key]

#===================================================================================================================================================
# functions
#===================================================================================================================================================

def printfastq(al,read):
    if(al.is_reverse):
        print "@" + str(al.qname) + "/" + str(read) + "\n" + str(revcomp(al.seq)) + "\n" + "+" + "\n" + str(al.qual[::-1])
    else: 
        print "@" + str(al.qname) + "/" + str(read) + "\n" + str(al.seq) + "\n" + "+" + "\n" + str(al.qual)

def printfastq_rg(al,read, rg):
    if(al.is_reverse):
        print "@" + str(al.qname) + "/" + str(read) + " " + "RG:Z:" + str(rg) + "\n" + str(revcomp(al.seq)) + "\n" + "+" + "\n" + str(al.qual[::-1])
    else: 
        print "@" + str(al.qname) + "/" + str(read) + " " + "RG:Z:" + str(rg) + "\n" + str(al.seq) + "\n" + "+" + "\n" + str(al.qual)

def revcomp(seq):
    seq1 = seq.translate(maketrans("AGCTagct", "TCGAtcga"))
    seq2 = seq1[::-1]
    return seq2

#===================================================================================================================================================
# parsing
#===================================================================================================================================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    
    usage = """%prog -i <file>

bamtofastq


    """
    parser = OptionParser(usage)
    
    parser.add_option("-i", "--bamfile", dest="bamfile", 
        help="A BAM file",
        metavar="FILE")

    (opts, args) = parser.parse_args()

    if opts.bamfile is None:
        parser.print_help()
        print
    else:
        try:
            bamtofastq(opts.bamfile)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return
if __name__ == "__main__":
    sys.exit(main()) 
    
# new features needed
# don't forget to select for primary alignments before using on bwa-mem files - DONE
# query memory usage to print out
# write output files based on library name
# include sample name as part of output file
# include sequencing center as part of output file
# stop opening and closing file the whole time - DONE
# close files at end
# add strict error checking ; need to catch all mistakes!!!!
# option for interleaved versus two fastq formats
# option for compressed versus uncompressed
# write out sequencing center - DONE
# need to make read-group based
# need to actually align output:
# need to have output directory prefix - default blank
