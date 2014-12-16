#!/usr/bin/env python

#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python

import pysam
import sys
import getopt
import string
# import gzip
from string import *
from optparse import OptionParser
import subprocess
import numpy

def bamProfiler(bamfile, concsize, quitnum):
    print "#---------------------------------------------------------------------------------------------------------------------------------------------------------------"
    fname = str(bamfile)
    fsize = strip(str(subprocess.check_output(["du",'-h',str(bamfile)]))).split("\t")[0]

# get file and header
    bam = pysam.Samfile(bamfile, "rb")
    bamheader=bam.header['RG']
        
# initiate variables
    num = 0
    rg_hash = {}
    frag_hash = {}
    read_hash = {}
    hist_hash = {}
# extract info from header
    for i in bamheader:
        rg_hash.setdefault(i['ID'],i)
        frag_hash.setdefault(i['ID'],[])
        read_hash.setdefault(i['ID'],[])

    for key in rg_hash:
        rg_hash[key]["num_reads"] = 0
        rg_hash[key]["is_mapped"] = 0
        rg_hash[key]["is_discordant"] = 0
        rg_hash[key]["is_secondary"] = 0
        rg_hash[key]["is_paired"] = 0
        rg_hash[key]["is_duplicate"] = 0
#       print "# rg_hash = " + str(rg_hash[key])

# fill potentially missing values that are later printed:
    for key in rg_hash:
        if 'SM' not in rg_hash[key]:
            rg_hash[key]['SM'] = "NA"
        if 'LB' not in rg_hash[key]:
            rg_hash[key]['LB'] = "NA"
        if 'ID' not in rg_hash[key]:
            rg_hash[key]['ID'] = "NA"
        if 'CN' not in rg_hash[key]:
            rg_hash[key]['CN'] = "NA"
        if 'DT' not in rg_hash[key]:
            rg_hash[key]['DT'] = "NA"

# go get alignments and collect stats
# things to count: number aligmnents each read group; number mapped; number discordants; number secondaries; readlen
    for al in bam.fetch():
        if num > quitnum: break
        for tag, val in al.tags:
            if tag=="RG": RG=val
        rg_hash[RG]['num_reads'] = rg_hash[RG]['num_reads'] + 1
        if al.is_unmapped == False:
            rg_hash[RG]['is_mapped'] = rg_hash[RG]['is_mapped'] + 1
        if al.is_proper_pair == False and al.is_unmapped == False:
            rg_hash[RG]['is_discordant'] = rg_hash[RG]['is_discordant'] + 1
        if al.is_secondary == True:
            rg_hash[RG]['is_secondary'] = rg_hash[RG]['is_secondary'] + 1
        if al.is_duplicate == True:
            rg_hash[RG]['is_duplicate'] = rg_hash[RG]['is_duplicate'] + 1
        if al.is_paired == True:
            rg_hash[RG]['is_paired'] = rg_hash[RG]['is_paired'] + 1
        if al.is_reverse == False and al.mate_is_reverse == True and al.tlen > 0 and al.tlen < concsize:
            frag_hash[RG].append(al.tlen)
        if al.is_secondary == False:
            read_hash[RG].append(len(al.seq))
        num = num + 1
        
    for key in frag_hash:
        rg_hash[key]["fraglen_mean"] = numpy.mean(frag_hash[key])
        rg_hash[key]["fraglen_SD"] = numpy.std(frag_hash[key])
        rg_hash[key]["fraglen_median"] = numpy.median(frag_hash[key])
        rg_hash[key]["fraglen_MAD"] = mad(frag_hash[key])
        hist_hash[key] = numpy.histogram(frag_hash[key],bins=1000,range=(0,1000))

    for key in read_hash:
        rg_hash[key]["readlen_mean"] = numpy.mean(read_hash[key])

    for key in frag_hash:
        for i in range(0,1000): 
            print "@INS" + "\t" + str(key) + "\t" + str(hist_hash[key][0][i]) + "\t" + str(hist_hash[key][1][i])

    for key in rg_hash:
        temp=[]
        for tag in rg_hash[key]:
            temp.append(str(tag) + ":" + str(rg_hash[key][tag]))
        print "@RG" + "\t" + '\t'.join([str(x) for x in temp]) + "\t" + "INPUT:" + str(bamfile)+ "\t" + "INPUT_SIZE:" + fsize

    print "@INFO" + "\t" + "#" + '1=SM' + "\t" + '2=LB' + "\t" + '3=ID' + "\t" + '4=CN' + "\t" + '5=DT' + "\t" + '6=PL' + "\t" + '7=INPUT' + "\t" + '8=INPUT_SIZE' \
    + "\t" + '9=num_reads' + "\t" + '10_is_paired' + "\t" + '11=is_mapped' + "\t" + '12=is_secondary' + "\t" + '13=is_discordant' + "\t" + '14=is_duplicate' \
    + "\t" + '15=fraglen_mean' + "\t" + '16=fraglen_SD' + "\t" + '17=fraglen_median' + "\t" + '18=fraglen_MAD' + "\t" + '19=readlen_mean'

    for key in rg_hash:
        x = rg_hash[key]
        p = [x['SM'],x['LB'],x['ID'],x['CN'],x['DT'],x['PL'],str(bamfile),fsize,x['num_reads'],x['is_paired'],x['is_mapped'],x['is_secondary'],x['is_discordant'],x['is_duplicate'],x['fraglen_mean'],x['fraglen_SD'],x['fraglen_median'],x['fraglen_MAD'],x['readlen_mean']]
        print "@INFO" + "\t" + "\t".join([str(y) for y in p])
            
#===================================================================================================================================================
# functions
#===================================================================================================================================================

def mad(a):
    residuals = []
    m = numpy.median(a)
    for val in a:
        residuals.append(abs(val - m))
    output = numpy.median(residuals)
    return output
    
#===================================================================================================================================================
# parsing
#===================================================================================================================================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    
    usage = """%prog -i <file>

extractBamInfo
Author: Ira Hall    
Description: per library bam to fastq from position-sorted bam

    """
    parser = OptionParser(usage)
    
    parser.add_option("-i", "--bamfile", dest="bamfile", 
        help="A BAM file",
        metavar="FILE")

    parser.add_option("-c", "--concsize", dest="concsize", default=1000, type = "int",
        help="max size of concordant readpairs for fragment size determination; default = 1000",
        metavar="INT")
        
    parser.add_option("-q", "--quitnum", dest="quitnum", default=10000000, type = "int",
        help="number of reads to examine to collect read group stats; default = 10000000",
        metavar="STR")
        
    (opts, args) = parser.parse_args()

    if opts.bamfile is None:
        parser.print_help()
        print
    else:
        try:
            bamProfiler(opts.bamfile, opts.concsize, opts.quitnum)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n");
            return
if __name__ == "__main__":
    sys.exit(main()) 

