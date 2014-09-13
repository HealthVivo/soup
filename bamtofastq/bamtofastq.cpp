#include <iostream>
#include <fstream>

#include "kstring.h"
#include "sam.h"

using namespace std;

char compBase (char base) {
  switch ( base ) {
  case 'A':
    return 'T';
    break;
  case 'C':
    return 'G';
    break;
  case 'G':
    return 'C';
    break;
  case 'T':
    return 'A';
    break;
  case 'N':
    return 'N';
    break;
  case 'a':
    return 't';
    break;
  case 'c':
    return 'g';
    break;
  case 'g':
    return 'c';
    break;
  case 't':
    return 'a';
    break;
  default:
    return 'N';
    break;
  }
}

int min (int a, int b) {
  if (a < b) { return a; }
  else { return b; }
}

int main( int argc, char *argv[] ) {
  // argc should be 2 for correct execution
  if ( argc != 2 ) {
    cout << "\n\
Program: bamtofastq\n\
Version: 0.0.1\n\
Authors: Ira Hall <ihall@genome.wustl.edu> and Colby Chiang <cc2qe@virginia.edu>\n\
Description: Convert a coordinate sorted BAM file to FASTQ, while preserving\n\
    read group information as a comment" << endl;
    cout << "\n\
Usage: " << argv[0] << " [bamFile]\n" << endl;
  }
  
  else {
    bam1_t *bamrec = bam_init1();
    samfile_t *samfile;
    if ( (samfile = samopen(argv[1], "rb", NULL)) == 0) {
      cout << "error: could not open bam file" << endl;
      return 1;
    }
    
    // open output streams
    ofstream fq1;
    fq1.open("1_ex.fq");
    
    ofstream fq2;
    fq2.open("2_ex.fq");
    
    int flag;

    while (samread(samfile, bamrec) > 0) {
      flag = bamrec->core.flag;
      // fprintf(stderr, "flag: %d", flag);

      int notPrimary = flag & 256;
      int firstInPair = flag & 64;
      int revStrand = flag & 16;
      
      if (! notPrimary) {
	if (firstInPair) {
	  fq1 << "@" << bam1_qname(bamrec) << endl;

	  if (!revStrand) {
	    int i;
	    for (i = 0; i < bamrec->core.l_qseq; ++i) {
	      fq1 << bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)] ;
	    }
	    fq1 << "\n+" << endl;

	    for (i = 0; i < bamrec->core.l_qseq; ++i) {
	      fq1 << (char) (bam1_qual(bamrec)[i] + 33) ;
	    }
	    fq1 << endl;
	  }

	  else {
	    int i;
	    for (i = bamrec->core.l_qseq - 1; i >= 0 ; --i) {
	      fq1 << compBase(bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)]) ;
	    }
	    fq1 << "\n+" << endl;

	    for (i = bamrec->core.l_qseq - 1; i >= 0; --i) {
	      fq1 << (char) (bam1_qual(bamrec)[i] + 33) ;
	    }
	    fq1 << endl;
	  }
	}
	else { // not first in pair
	  fq2 << "@" << bam1_qname(bamrec) << endl;
	  if (!revStrand) {
	    int i;
	    for (i = 0; i < bamrec->core.l_qseq; ++i) {
	      fq2 << bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)] ;
	    }
	    fq2 << "\n+" << endl;

	    for (i = 0; i < bamrec->core.l_qseq; ++i) {
	      fq2 << (char) (bam1_qual(bamrec)[i] + 33) ;
	    }
	    fq2 << endl;
	  }

	  else {
	    int i;
	    for (i = bamrec->core.l_qseq - 1; i >= 0 ; --i) {
	      fq2 << compBase(bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)]) ;
	    }
	    fq2 << "\n+" << endl;

	    for (i = bamrec->core.l_qseq - 1; i >= 0; --i) {
	      fq2 << (char) (bam1_qual(bamrec)[i] + 33) ;
	    }
	    fq2 << endl;
	  }
	}
      }
    }
    fq1.close();
    fq2.close();
  }
  return 0;
}

