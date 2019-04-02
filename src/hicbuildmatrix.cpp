#include "hicbuildmatrix.hpp"

HiCBuildMatrix::HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead) {
    mForwardRead = pForwardRead;
    mReverseRead = pReverseRead;
}

size_t HiCBuildMatrix::readBamFile() {
    std::cout << "forward read: " << mForwardRead << std::endl;
    std::cout << "reverse read: " << mReverseRead << std::endl;
    return 0;
}