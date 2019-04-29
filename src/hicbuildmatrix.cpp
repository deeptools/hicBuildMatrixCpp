#include "hicbuildmatrix.hpp"

HiCBuildMatrix::HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead) {
    mForwardRead = pForwardRead;
    mReverseRead = pReverseRead;
}
std::set<string> pos_matrix;
bool HiCBuildMatrix::is_duplicacted(const std::string &pchrom1, const std::int &pstart1, const std::string &pchrom2,const std::int &pstart2) {
   mchrom1 = pchrom1;
   mstart1 = pstart1;
   mchrom2 = pchrom2;
   mstart2 = pstrt2;
    if(mchrom1< mchrom2) {
        string id_string= Format("___ - ___", "mchrom1","mchrom2");
    }
    else {
        string id_string=Format("___ - ___", "mchrom2","mchrom1");
    }
    if (mstart1 < mstart2){
        string id_string= Format("___ - ___", "mstart1","mstart2");
    }
    else {
         string id_string= Format("___ - ___", "mstart2","mstart1");
    }
    if (pos_matrix.count(id_string)) {
        return 1;
   // id_string is in the set, count is 1
} else {
    pos_matrix.insert(id_string);
   // count zero, i.e. id_string not in the set
}
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
    bool pSkipDuplicationCheck;
    read_pos_matrix = HiCBuildMatrix();
    pReadPosMatrix = read_pos_matrix;

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
 seqan::
  record2;
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
  if(pSkipDuplicationCheck == 0){
            if pReadPosMatrix.is_duplicated(pRefId2name[record1.rname], record1.pos,  
            pRefId2name[record2.rname],record2.pos){
               duplicated_pairs += 1; 

            }

  }  
    return 0; 
}

