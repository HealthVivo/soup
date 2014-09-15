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

// Hash table from http://stackoverflow.com/a/4384446
struct nlist { /* table entry: */
  struct nlist *next; /* next entry in the chain */
  char *name; /* defined name */
  bam1_t *defn; /* replacement test */
};

#define HASHSIZE 15485863 /* best if a prime number */
static struct nlist *hashtab[HASHSIZE]; /* pointer table */

/* hash: form hash value for string s */
unsigned hash(char *s) {
  unsigned hashval;
  for (hashval = 0; *s != '\0'; s++)
    hashval = *s + 31 * hashval;
  return hashval % HASHSIZE;
}

unsigned next_hash(unsigned h) {
  unsigned next_h;
  while (hashtab[h] != NULL) {
    h = (h+1) % HASHSIZE;
  }
  return h;
}

/* lookup: look for s in hashtab */
struct nlist *lookup(char *s) {
  struct nlist *np;
  for (np = hashtab[hash(s)]; np != NULL; np = np->next)
    if (strcmp(s, np->name) == 0)
      return np; /* found */
    // else {
    //   cout << "collision \n" << np->name << "\n";
    //   cout << np->defn->core.flag << "\n";
    // }
  return NULL; /* not found */
}

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
    // free((void *) np->name); /* free previous name */
    bam_destroy1(np->defn); /* free previous defn */
    // free((void *) np->defn); /* free previous defn */
  if ((np->defn = bam_dup1(defn)) == NULL)
    return NULL;
  return np;
}

int *remove(char *name) {
  struct nlist *np;
  if ((np = lookup(name)) == NULL) { /* not found */
    return NULL;
  } else
    free((void *) np->defn);
}

int print_fastq(bam1_t *bamrec, char read) {
  int revStrand = bamrec->core.flag & 16;
  int i;
  uint8_t *readgroup = bam_aux_get(bamrec, "RG");
  cout << "@" << bam1_qname(bamrec) << '/' << read << " RG:Z:" << readgroup+1 << '\n'; // +1 omits the first character (Z) indicating string
  if (!revStrand) { // forward strand
    for (i = 0; i < bamrec->core.l_qseq; ++i) {
      cout << bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)] ;
    }
    cout << "\n+\n";
    
    for (i = 0; i < bamrec->core.l_qseq; ++i) {
      cout << (char) (bam1_qual(bamrec)[i] + 33) ;
    }
    cout << '\n';
  }
  else { // reverse strand
    for (i = bamrec->core.l_qseq - 1; i >= 0 ; --i) {
      cout << compBase(bam_nt16_rev_table[bam1_seqi(bam1_seq(bamrec),i)]) ;
    }
    cout << "\n+\n";
	  
    for (i = bamrec->core.l_qseq - 1; i >= 0; --i) {
      cout << (char) (bam1_qual(bamrec)[i] + 33) ;
    }
    cout << '\n';
  }
  return 0;
}

int main( int argc, char *argv[] ) {
  // argc should be 2 for correct execution
  if ( argc != 2 ) {
    cout << "\n\
Program: bamtofastq\n\
Version: 0.0.1\n\
Authors: Ira Hall <ihall@genome.wustl.edu> and Colby Chiang <cc2qe@virginia.edu>\n\
Description: Convert a coordinate sorted BAM file to FASTQ, while preserving\n\
    read group information as a comment" << '\n';
    cout << "\n\
Usage: " << argv[0] << " <bamFile>\n" << '\n';
  }
  
  else {
    bam1_t *bamrec = bam_init1();
    samfile_t *samfile;
    if ( (samfile = samopen(argv[1], "rb", NULL)) == 0) {
      cout << "error: could not open bam file\n";
      return 1;
    }

    int notPrimary, firstInPair;
    while (samread(samfile, bamrec) > 0) {
      notPrimary = bamrec->core.flag & 256;

      // bam1_t *bamrec = malloc(sizeof(bam1_t));
      if (notPrimary) continue; // skip secondary alignments

      // add to hash table if mate already there
      nlist *mate = lookup(bam1_qname(bamrec));
      // cout << "lookedup" << bam1_qname(bamrec) << "\n";
      if (! mate) {
      	install(bam1_qname(bamrec), bamrec);

      	continue;
      }
      // print the FASTQ for matched pairs
      firstInPair = bamrec->core.flag & 64;
      if (firstInPair) {
      	print_fastq(bamrec, '1');
      	print_fastq(mate->defn, '2');
      }
      else {
      	print_fastq(mate->defn, '1');
      	print_fastq(bamrec, '2');
      }
      
      bam_destroy1(mate->defn);
      

      // maybe this will free memory?
      // install(bam1_qname(bamrec), bamrec);

    }
  }
  return 0;
}

