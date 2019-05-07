#ifndef HICBUILDMATRIX_HPP
#define HICBUILDMATRIX_HPP

#include <iostream>
#include <vector> 
#include <string>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <assert.h> 
#include <tuple> 
#include <unordered_map> 

class HiCBuildMatrix
{
  public:
    HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead);
    size_t readBamFile(std::int &pNumberOfItemsPerBuffer,std::bool &pSkipDuplicationCheck, std::tuple<T...> &pRefId2name);
);
    bool is_duplicacted(const std::string &pchrom1, const std::int &pstart1,const std::string &pchrom2,  const std::int &pstart2); 
  private:
    std::string mForwardRead;
    std::string mReverseRead;
    std::string mchrom1;
    std::string mchrom2;
    std::int mstart1;;
    std::int mstart2;
    std::int mNumberOfItemsPerBuffer;
    Std::bool mSkipDublicationCheck;
};

#endif