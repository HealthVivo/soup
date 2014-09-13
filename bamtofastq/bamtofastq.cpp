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


// Dictionary from http://stackoverflow.com/a/4384446
struct nlist { /* table entry: */
  struct nlist *next; /* next entry in the chain */
  char *name; /* defined name */
  bam1_t *defn; /* replacement test */
};

#define HASHSIZE 1000001
  static struct nlist *hashtab[HASHSIZE]; /* pointer table */

/* hash: form hash value for string s */
unsigned hash(char *s) {
  unsigned hashval;
  for (hashval = 0; *s != '\0'; s++)
    hashval = *s + 31 * hashval;
  return hashval % HASHSIZE;
}

/* lookup: look for s in hashtab */
struct nlist *lookup(char *s) {
  struct nlist *np;
  for (np = hashtab[hash(s)]; np != NULL; np = np->next)
    if (strcmp(s, np->name) == 0)
      return np; /* found */
  return NULL; /* not found */
}

bam1_t *bamrec_dup(bam1_t *);
/* install: put (name, defn) in hashtab */
struct nlist *install(char *name, bam1_t *defn) {
  struct nlist *np;
  unsigned hashval;
  if ((np = lookup(name)) == NULL) { /* not found */
    np = (struct nlist *) malloc(sizeof(*np));
    if (np == NULL || (np->name = strdup(name)) == NULL)
      return NULL;
    hashval = hash(name);
    np->next = hashtab[hashval];
    hashtab[hashval] = np;
  } else /* already there */
    free((void *) np->defn); /* free previous defn */
  if ((np->defn = bamrec_dup(defn)) == NULL)
    return NULL;
  return np;
}

bam1_t *bamrec_dup(bam1_t *b) { /* make a duplicate of b */
  bam1_t *p;
  p = (bam1_t *) malloc(sizeof(bam1_t)); /* +1 for '\0' */
  if (p != NULL)
    p = b;
  return p;
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


    int flag, notPrimary, firstInPair, revStrand, i;
    int j = 0;
    while (samread(samfile, bamrec) > 0) {
      flag = bamrec->core.flag;
      // fprintf(stderr, "flag: %d", flag);

      notPrimary = flag & 256;
      firstInPair = flag & 64;
      revStrand = flag & 16;

      if (! notPrimary) {
	// install("bam1_qname(bamrec)", 'string');

	install(bam1_qname(bamrec), bamrec);
	// cout << lookup("hello")->defn << endl;
	j += 1;

	cout << lookup("D2D6HACXX130815:4:2309:17632:73110")->defn->core.l_qseq << endl;

	if (firstInPair) {
	  cout << "@" << bam1_qname(bamrec) << " RG:Z:" << bam_aux_get(bamrec, "RG") << endl;
	  if (!revStrand) {
	    for (i = 0; i < bamrec->core.l_qseq; ++i) {
	      cout << bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)] ;
	    }
	    cout << "\n+" << endl;

	    for (i = 0; i < bamrec->core.l_qseq; ++i) {
	      cout << (char) (bam1_qual(bamrec)[i] + 33) ;
	    }
	    cout << endl;
	  }

	  else {
	    for (i = bamrec->core.l_qseq - 1; i >= 0 ; --i) {
	      cout << compBase(bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)]) ;
	    }
	    cout << "\n+" << endl;

	    for (i = bamrec->core.l_qseq - 1; i >= 0; --i) {
	      cout << (char) (bam1_qual(bamrec)[i] + 33) ;
	    }
	    cout << endl;
	  }
	}
	else { // not first in pair
	  cout << "@" << bam1_qname(bamrec) << " RG:Z:" << bam_aux_get(bamrec, "RG") << endl;
	  if (!revStrand) {
	    for (i = 0; i < bamrec->core.l_qseq; ++i) {
	      cout << bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)] ;
	    }
	    cout << "\n+" << endl;

	    for (i = 0; i < bamrec->core.l_qseq; ++i) {
	      cout << (char) (bam1_qual(bamrec)[i] + 33) ;
	    }
	    cout << endl;
	  }

	  else {
	    for (i = bamrec->core.l_qseq - 1; i >= 0 ; --i) {
	      cout << compBase(bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)]) ;
	    }
	    cout << "\n+" << endl;

	    for (i = bamrec->core.l_qseq - 1; i >= 0; --i) {
	      cout << (char) (bam1_qual(bamrec)[i] + 33) ;
	    }
	    cout << endl;
	  }
	}
      }
    }
  }
  return 0;
}

