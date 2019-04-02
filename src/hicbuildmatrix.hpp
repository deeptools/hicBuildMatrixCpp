#ifndef HICBUILDMATRIX_HPP
#define HICBUILDMATRIX_HPP

#include <omp.h>

#include <iostream>
#include <string>
#include <seqan/bam_io.h>

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