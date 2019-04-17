#ifndef HICBUILDMATRIX_HPP
#define HICBUILDMATRIX_HPP

// #include <omp.h>
// #include <seqan/file.h>
// using namespace seqan;
// #include <iostream>

#include <iostream>
#include <vector> 
#include <string>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

class HiCBuildMatrix
{
  public:
    HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead);
    size_t readBamFile();

  private:
    std::string mForwardRead;
    std::string mReverseRead;

};

#endif