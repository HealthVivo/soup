#!/bin/sh
export ROOTSYS=../../root
export SAMDIR=samtools
export DYLD_LIBRARY_PATH=../../root/lib
./cnvnator -tree /Users/hiaips/Desktop/research/NA12878_S1.19.bam -genome GRCh37 -root tree.root -unique
