#ifndef __AliParser2__
#define __AliParser2__

#include <string>
#include <sam.h>
#include <sys/types.h>
using std::string;

class AlignmentData
{
    public:
    AlignmentData(const bam1_t *data) { this->data = data; }
    
    inline string getQueryName() { return (data) ? bam1_qname(data) : ""; }
    inline bool isUnmapped()     { return data->core.flag & 0x4; }
    inline bool isNextUnmapped() { return data->core.flag & 0x8; }
    inline bool isReversed()     { return data->core.flag & 0x10; }
    inline bool isNextReversed() { return data->core.flag & 0x20; }
    inline bool isDuplicate()    { return data->core.flag & 0x400; }
    inline int32_t getChromosomeIndex() { return data->core.tid; }
    inline int32_t getStart() { return data->core.pos + 1; }
    inline int32_t getEnd() { return bam_calend(&data->core, bam1_cigar(data)); }
    inline int32_t getReadLength() { return data->core.l_qseq; }
    inline int32_t getFragmentLength() { return data->core.isize; }
    inline int32_t getQuality() { return data->core.qual; }
    inline bool isQ0() { return (data->core.qual <= 1); }
    
    private:
    const bam1_t *data;
};


class AliParser2
{
    public:
    AliParser2(const string &filename);
    virtual ~AliParser2();
    int parseRegion(bam_fetch_f callback, const string &chunk, void *data);
    int produceTrees(string *user_chroms,int n_chroms, string *user_files,int n_files, bool forUnique);
    
    private:
    bam_index_t *bamidx_;
    samfile_t *samfp_;
};

#endif
