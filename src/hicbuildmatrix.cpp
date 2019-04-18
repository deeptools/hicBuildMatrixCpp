#include "hicbuildmatrix.hpp"

HiCBuildMatrix::HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead) {
    mForwardRead = pForwardRead;
    mReverseRead = pReverseRead;
}

size_t HiCBuildMatrix::readBamFile() {
    
   // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn1("R1.sam");
    seqan::BamStream bamStreamIn2("R2.sam");
    int duplicated_pairs = 0;
    int one_mate_unmapped = 0;
    int one_mate_not_unique = 0;
    int one_mate_low_quality = 0;
    bool  all_data_read=1;
    if (!isGood(bamStreamIn1))
    {
        std::cerr << "ERROR: Could not open the first file!\n";
        return 1;
    }
     else if (!isGood(bamStreamIn2))
    {
        std::cerr << "ERROR: Could not open the second file!\n";
        return 1;
    }
    
    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamStreamOut1("-", seqan::BamStream::WRITE);
    seqan::BamStream bamStreamOut2("-", seqan::BamStream::WRITE);
  
    // Copy header.  The header is automatically written out before
    // the first record.
    bamStreamOut1.header = bamStreamIn1.header;
    bamStreamOut2.header = bamStreamIn2.header;

    seqan::BamAlignmentRecord record1;
 seqan::BamAlignmentRecord record2;
    while (!atEnd(bamStreamIn1) && !atEnd(bamStreamIn2))
    {
        readRecord(record1, bamStreamIn1);
        while ((hasFlagAllProper(record1)==1) & 256 == 256) {
            readRecord(record1, bamStreamIn1);
            
        }
          all_data_read =1;
            break; 
        readRecord(record2, bamStreamIn2);
        while ((hasFlagAllProper(record2)==1)  & 256 == 256) {
            readRecord(record2, bamStreamIn1);    
        }
        all_data_read = 1;
            break; 

    }
assert(record1.qName == record2.qName && "Be sure that the sam files have the same read order ");
// if any of the reads is not mapped
   if(((hasFlagAllProper(record1)==1) && 0x4 == 4)|| ((hasFlagAllProper(record2)==1) && 0x4 == 4)){
       one_mate_unmapped += 1;
   }
   // if the read quality is low, higher probabilities of error
    if(record1.mapQ < pMinMappingQuality || record2.mapQ < pMinMappingQuality){
        if(record1.mapQ ==0 && record2.mapQ==0) {
            one_mate_not_unique += 1;
        }
    }
    
    return 0; 
}

