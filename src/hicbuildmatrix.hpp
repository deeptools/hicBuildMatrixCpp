#ifndef HICBUILDMATRIX_HPP
#define HICBUILDMATRIX_HPP

#include <iostream>
#include <vector> 
#include <string>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <assert.h> 

class HiCBuildMatrix
{
  public:
    HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead);
    size_t readBamFile();
    bool is_duplicacted(chrom1, start1, chrom2, start2); 
  private:
    std::string mForwardRead;
    std::string mReverseRead;
    std::string chrom1;
    std::string chrom2;
    std::int start1;;
    std::int start2;

};

#endif