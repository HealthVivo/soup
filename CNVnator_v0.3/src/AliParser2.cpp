#include "AliParser2.hh"
#include <iostream>
#include <bam.h>
#include <vector>

using namespace std;

AliParser2::AliParser2(const string &filename)
{
    this->bamidx_ = NULL;
    this->samfp_ = NULL;
    
    // load index
    this->bamidx_ = bam_index_load(filename.c_str());
    if (this->bamidx_ == NULL) {
        cerr << "Error loading BAM index file" << endl;
        return;
    }
    
    // load SAM file
    this->samfp_ =  samopen(filename.c_str(), "rb", "");
    if (this->samfp_ == NULL) {
        cerr << "Error loading BAM file" << endl;
        return;
    }
}

AliParser2::~AliParser2()
{
    if (this->bamidx_ != NULL) {
        bam_index_destroy(this->bamidx_);
        this->bamidx_ = NULL;
    }
    if (this->samfp_ != NULL) {
        samclose(this->samfp_);
        this->samfp_ = NULL;
    }
}

int AliParser2::parseRegion(bam_fetch_f callback, const string &chunk, void *data)
{
    int tid = 0, beg = 0, end = 0;
    if (bam_parse_region(this->samfp_->header, chunk.c_str(), &tid, &beg, &end) < 0) {
        return -1;
    }
    bam_fetch(this->samfp_->x.bam, this->bamidx_, tid, beg, end, data, callback);
    return 0;
}
