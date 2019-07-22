#include "hicbuildmatrix.hpp"
 //  Hi-C quantifies interactions between all possible pairs of fragments simultaneously.
HiCBuildMatrix::HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead)
{
    mForwardRead = pForwardRead;
    mReverseRead = pReverseRead;
}

bool HiCBuildMatrix::is_duplicated(seqan::CharString pChrom1, size_t pStart1, seqan::CharString pChrom2, size_t pStart2)
{

    std::string sequence;
    if (pChrom1 < pChrom2)
    {
        sequence.append(seqan::toCString(pChrom1));
        sequence.append("-");

        sequence.append(seqan::toCString(pChrom2));
    }
    else
    {
        sequence.append(seqan::toCString(pChrom2));
        sequence.append("-");
        sequence.append(seqan::toCString(pChrom1));
    }
    sequence.append("-");

    if (pStart1 < pStart2)
    {
        sequence.append(std::to_string(pStart1));
        sequence.append("-");
        sequence.append(std::to_string(pStart2));
    }
    else
    {
        sequence.append(std::to_string(pStart2));
        sequence.append("-");
        sequence.append(std::to_string(pStart1));
    }
    
    std::cout << 'id_string' << sequence << std::endl;
    auto it = mPosMatrix.find(sequence);
    if (it != mPosMatrix.end()) return true;
    mPosMatrix.insert(sequence);

    return false;
}  

std::vector<std::string> HiCBuildMatrix::createRefId2name(seqan::BamStream pBamStream)
{
    std::vector<std::string> chromosome_refName; // vector with chromosomes
    std::vector<size_t> chromosome_size;
    for (size_t i = 0; i < length(pBamStream.header.sequenceInfos); ++i)
    {

        chromosome_refName.push_back(seqan::toCString(pBamStream.header.sequenceInfos[i].i1));

        chromosome_size.push_back(pBamStream.header.sequenceInfos[i].i2);
    }
    for (size_t j = 0; j < chromosome_refName.size(); j++) // read names of chromosomes
    {
        std::cout << chromosome_refName[j] << ' ' << chromosome_size[j] << std::endl;
        
    }
}
      
size_t HiCBuildMatrix::readBamFile(int pNumberOfItemsPerBuffer, bool pSkipDuplicationCheck, std::vector<std::string> pRefId2name)
{

    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn1("R1.sam");
    seqan::BamStream bamStreamIn2("R2.sam");
    std::vector<seqan::BamAlignmentRecord> buffer_mate1;
    std::vector<seqan::BamAlignmentRecord> buffer_mate2;
    std::vector<seqan::BamAlignmentRecord> mate_bins;
    //std::unordered_map<std::string,IntervalTree<size_t, size_t>> pSharedBinIntvalTree;
    int duplicated_pairs = 0; 
    int one_mate_unmapped = 0;
    int one_mate_not_unique = 0;
    int one_mate_low_quality = 0;
    bool all_data_read = 0;
    int j = 0;
    int iter_num = 0;
    mSkipDublicationCheck = pSkipDuplicationCheck;
    bool mate_is_unasigned = 0;
    int pMinMappingQuality;
    bool pKeepSelfCircles;
    std::string pRestrictionSequence;
    int pMatrixSize;
    int pBinSize = 10000;
    int pResultIndex;
    int start = 0;
    int end;

    int pair_added = 0;
    int  dangling_end = 0;
    int self_circle = 0;
    int self_ligation = 0;
    int same_fragment = 0;
    int mate_not_close_to_rf = 0;
    int count_inward = 0;
    int count_outward = 0;
    int count_left = 0;
    int count_right = 0;
    int inter_chromosomal = 0;
    int short_range = 0;
    int long_range = 0;
    int row[] = {};
    int col[] = {};
    int data[] = {};

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
    std::cout << __LINE__ << std::endl;
    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamStreamOut1("-", seqan::BamStream::WRITE);
    seqan::BamStream bamStreamOut2("-", seqan::BamStream::WRITE);
    std::cout << __LINE__ << std::endl;

    // Copy header.  The header is automatically written out before
    // the first record.
    bamStreamOut1.header = bamStreamIn1.header;
    bamStreamOut2.header = bamStreamIn2.header;

    seqan::BamAlignmentRecord record1;
    seqan::BamAlignmentRecord record2;
     while (j < pNumberOfItemsPerBuffer)
     {
         readRecord(record1, bamStreamIn1);
         readRecord(record2, bamStreamIn2);
        j++;
     }
    all_data_read = 1;
    iter_num += 1;
    std::cout << __LINE__ << std::endl; 
    std::unordered_map<std::string,IntervalTree<size_t,size_t>> pSharedBinIntvalTree = HiCBuildMatrix::createInitialStructures(bamStreamIn1,pBinSize);

    while (!atEnd(bamStreamIn1) && !atEnd(bamStreamIn2))
    {
        readRecord(record1, bamStreamIn1);
        while ((hasFlagAllProper(record1) == 1) & 256 == 256)
        {
            readRecord(record1, bamStreamIn1);
        }
        // all_data_read = 1;
        // break;
        readRecord(record2, bamStreamIn2);
        while ((hasFlagAllProper(record2) == 1) & 256 == 256)
        {
            readRecord(record2, bamStreamIn2);
        }
    
        std::cout << __LINE__ << std::endl;

        assert(record1.qName == record2.qName && "Be sure that the sam files have the same read order ");
        // if any of the reads is not mapped
        if (hasFlagAllProper(record1) || hasFlagAllProper(record2))
        {
            one_mate_unmapped += 1;
            continue;
        }
        // if the read quality is low, higher probabilities of error
        if (record1.mapQ < pMinMappingQuality || record2.mapQ < pMinMappingQuality)
        {
            if (record1.mapQ == 0 && record2.mapQ == 0)
            {
                one_mate_not_unique += 1;
            }
        }
        if (pSkipDuplicationCheck == 0)
        {
            
            if (is_duplicated(record1.qName, record1.beginPos, // chromosome_refName
                            record2.qName, record2.beginPos))
            {
                duplicated_pairs += 1;
            }
        }
        std::cout << __LINE__ << std::endl;
        std::cout << duplicated_pairs << std::endl;

         
         buffer_mate1.push_back(record1);
         buffer_mate2.push_back(record2);
         j += 1;
         if (((all_data_read != 0) && (getAlignmentLengthInRef(record1) != 0)) && getAlignmentLengthInRef(record2) != 0)
         {
             return buffer_mate1, buffer_mate2, 1, duplicated_pairs,
                    one_mate_unmapped, one_mate_not_unique, one_mate_low_quality,
                   iter_num - buffer_mate1.size();
         }
         if ((all_data_read == 0 && buffer_mate1.size() == 0) || (buffer_mate2.size() == 0))
         {
             return 0, 0, 1, duplicated_pairs,
                   one_mate_unmapped, one_mate_not_unique, one_mate_low_quality,
                    iter_num - buffer_mate1.size();
         }
         else
             return buffer_mate1, buffer_mate2, 0, duplicated_pairs, one_mate_unmapped,
                    one_mate_not_unique, one_mate_low_quality, iter_num - buffer_mate1.size();
         while (!atEnd(bamStreamIn1) && !atEnd(bamStreamIn2))
        {
            readRecord(record1, bamStreamIn1);
            seqan::CharString mate_ref1 = record1.qName; //get<record1.rname>(pRefId2name);
            readRecord(record2, bamStreamIn2);
            seqan::CharString mate_ref2 = record2.qName; //get<record2.rname>(pRefId2name);
            int read_middle = record1.beginPos + int(getAlignmentLengthInRef(record1) / 2);
            int middle_pos = int((start + end) / 2);

             /*std::vector<int> middle_position_element;
             int k = pSharedBinIntvalTree[middle_pos];
             for (int i = 0; i < middle_position_element.size(); i++)
             {
                 for (int j = 0; j < pSharedBinIntvalTree.size(); j++)
                 {
                    if (pSharedBinIntvalTree[j] == k)
                    {
                         middle_position_element.insert(k);
                    }
                 }
             }
             second attemp
              interval_T = pSharedBinIntvalTree.createInitialStructures(bamStreamIn1,pBinSize);
              unordered_map<std::string, IntervalTree<size_t, size_t>>:: iterator p; 
              p=pSharedBinIntvalTree.begin()
              while(p!=pSharedBinIntvalTree.size()/2) 
              {
                  p++;
              }
              p->second;
             */

             while (!start > end)
             {   int i;
                if ( (middle_position_element[1]) <= read_middle && (read_middle <= middle_position_element[i]))
                {
                     std::string mate_bin = std::to_string(pSharedBinIntvalTree[middle_pos]);
                     mate_is_unasigned = 0;
                }
                 else if (middle_position_element[1] <= read_middle)
                 {
                     end = middle_pos - 1;
                 middle_pos = int((start + end) / 2);
                     mate_is_unasigned = 1;
                 }
                 else
                {
                    start = middle_pos + 1;
                     middle_pos = int((start + end) / 2);
                     mate_is_unasigned = 1;
                 }
             }
             try
             {
                 throw mate_is_unasigned = 1;
             }
            catch (...)
             {
                 std::cout << "An exception occurred.";
             }
             if (mate_bins.empty())
             {
                 mate_is_unasigned = 1;
                 break;
             }  
        }
     }

    return 0; 
}  

std::unordered_map<std::string, IntervalTree<size_t, size_t>> HiCBuildMatrix::createInitialStructures(seqan::BamStream pBamStream, int pBinSize)
{


    std::vector<std::string> chromosome_refName; // vector with chromosomes
    std::vector<size_t> chromosome_size;         // vector with sizes

    for (size_t i = 0; i < length(pBamStream.header.sequenceInfos); ++i)
    {

        chromosome_refName.push_back(seqan::toCString(pBamStream.header.sequenceInfos[i].i1));

        chromosome_size.push_back(pBamStream.header.sequenceInfos[i].i2);
        
    }
    for (size_t j = 0; j < chromosome_refName.size(); j++) // read names of chromosomes
    {
        std::cout << chromosome_refName[j] << ' ' << chromosome_size[j] << std::endl;
       
    }

    std::unordered_map<std::string, IntervalTree<size_t, size_t>> intervalTree;
    size_t numberOfChromosomes = chromosome_size.size(); //  get number of chromosomes
    size_t lengthOfChromosome = 0;
    size_t binSize = 10000; // TODO: get real bin size as a parameter from python
    size_t counter = 0;
    std::string chromosomeName;

    std::cout << "numberOfChromosomes: " << numberOfChromosomes << std::endl;
    for (size_t i = 0; i < numberOfChromosomes; ++i)
    {
        std::vector<Interval<size_t, size_t>> intervals;

        lengthOfChromosome = chromosome_size[i]; //YES// get real value for the specific chromosome

        counter = 0;

        chromosomeName = chromosome_refName[i]; //get real chromosome name for specific chromosome

        for (size_t j = 0; j < lengthOfChromosome; j = j + binSize)
        {

            intervals.push_back(Interval<size_t, size_t>(j, j + binSize, counter));

            counter++;
        }
        // std::cout << "intervals siye: " << intervals.size() << std::endl;   NO
        IntervalTree<size_t, int> tree;
        intervalTree[chromosomeName] = IntervalTree<size_t, size_t>(std::move(intervals));
    }

    // intervals.push_back(Interval<size_t, int>(2, 10, 1));  NO
    // intervals.push_back(Interval<size_t, int>(3, 4, 2));
    // intervals.push_back(Interval<size_t,int>(20, 100, 3));
    // IntervalTree<size_t, int> tree;
    // tree = IntervalTree<size_t, int>(std::move(intervals));

    // std::vector<Interval<size_t, int> > results;
    std::string chromosome("chr2L");
    auto results = intervalTree[chromosome].findOverlapping(999, 99999);
    std::cout << "found " << results.size() << " overlapping intervals" << std::endl;
    for (size_t i = 0; i < results.size(); i++)
    {
        std::cout << results[i] << std::endl;
    }
    return intervalTree;
} 
