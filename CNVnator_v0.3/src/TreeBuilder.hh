#ifndef __TreeBuilder__
#define __TreeBuilder__

#include <iostream>
#include <string>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THashTable.h>
#include <stdint.h>
#include "Genome.hh"
using std::string;

class TreeBuilder
{
    public:
    TreeBuilder(Genome *genome, string root_fn, int max_threads = 32, bool forUnique = true);
    ~TreeBuilder();
    int build(string *user_chroms, int n_chroms, string user_files);
    bool syncData();
    
    private:
    void writeTreeForChromosome(string chrom, short *arr_p, short *arr_u, int len);
    int findIndex(string *arr, int n, string name);
    int findIndex(const vector<string> &arr, string name);
    void buildSerial(const string &fn, const vector<pair<string, size_t> > &chroms);
    void buildParallel(const string &fn, const vector<pair<string, size_t> > &chroms);
    
    private:
    Genome *refGenome_;
    string root_;
    int max_threads_;
    bool forUnique_;
    
    // histograms
    vector<void *> counts_;
};

#endif
