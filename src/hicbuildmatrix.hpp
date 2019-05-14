#ifndef HICBUILDMATRIX_HPP
#define HICBUILDMATRIX_HPP

#include <iostream>
#include <vector> 
#include <string>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <assert.h> 
#include <tuple> 
#include <unordered_map> 
#include "IntervalTree.h"

class HiCBuildMatrix
{
  public:
    HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead);
    size_t readBamFile(int &pNumberOfItemsPerBuffer, bool &pSkipDuplicationCheck, std::vector<int > &pRefId2name);

    bool is_duplicacted(const std::string &pchrom1, const int &pstart1,const std::string &pchrom2,  const int &pstart2); 
    std::unordered_map<std::string, IntervalTree<size_t, size_t> > createInitialStructures(seqan::BamStream pBamStream, int pBinSize);
  private:
    std::string mForwardRead;
    std::string mReverseRead;
    std::string mchrom1;
    std::string mchrom2;
    int mstart1;
    int mstart2;
    int mNumberOfItemsPerBuffer;
    bool mSkipDublicationCheck;
};

#endif