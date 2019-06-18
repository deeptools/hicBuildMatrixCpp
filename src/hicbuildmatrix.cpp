#include "hicbuildmatrix.hpp"

HiCBuildMatrix::HiCBuildMatrix(const std::string &pForwardRead, const std::string &pReverseRead)
{
    mForwardRead = pForwardRead;
    mReverseRead = pReverseRead;
}

bool HiCBuildMatrix::is_duplicated(seqan::CharString pChrom1, size_t pStart1, seqan::CharString pChrom2, size_t pStart2)
{

    std::string sequence;
    // seqan::CharString mchrom1 = pChrom1;
    // int mstart1 = pStart1;
    //std::string  = std::to_string(mstart1);
    // seqan::CharString mchrom2 = pChrom2;
    // int mstart2 = pstart2;
    //std::string  = std::to_string(mstart2);
    // seqan::CharString id_string;

    if (pChrom1 < pChrom2)
    {
        sequence.append(seqan::toCString(pChrom1));
        // sequence.insert(pChrom2);
        sequence.append("-");

        sequence.append(seqan::toCString(pChrom2));
    }
    else
    {
        sequence.append(seqan::toCString(pChrom2));
        // sequence.insert(pChrom2);
        sequence.append("-");

        sequence.append(seqan::toCString(pChrom1));
        // seqan::CharString id_string = concat(pchrom2, '-', pchrom1); //Format("___ - ___", "mchrom2","mchrom1");
    }
    sequence.append("-");

    if (pStart1 < pStart2)
    {
        sequence.append(std::to_string(pStart1));
        sequence.append("-");
        sequence.append(std::to_string(pStart2));

        // sequence.insert(pStart2);

        // seqan::CharString id_string =  - mstart2; //std::to_string(mstart1 - mstart2); //Format("___ - ___", "mstart1","mstart2");
    }
    else
    {
        sequence.append(std::to_string(pStart2));
        sequence.append("-");
        sequence.append(std::to_string(pStart1));

        // sequence.insert(std::to_string(pStart1));
        // sequence.insert('-');
        // sequence.insert(std::to_string(pStart2));

        // seqan::CharString id_string = mstart2 - mstart1; //std::to_string(mstart2 - mstart1); //Format("___ - ___", "mstart2","mstart1");
    }
    // seqan::CharString id_string = concat(sequence); //Format("___ - ___", "mchrom1","mchrom2");
    std::cout << 'id_string' << sequence << std::endl;
    auto it = mPosMatrix.find(sequence);
    if (it != mPosMatrix.end()) return true;
    mPosMatrix.insert(sequence);

    return false;
    // for (it = mPosMatrix.begin(); it != mPosMatrix.end(); ++it)
    // {
    //     return true;
    // }
    // //if (pos_matrix.count(id_string)) {
    // // id_string is in the set, count is 1

    // // count zero, i.e. id_string not in the set

    // return true;
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
        //   } // read sizes of chromosomes
    }
}
size_t HiCBuildMatrix::readBamFile(int pNumberOfItemsPerBuffer, bool pSkipDuplicationCheck, std::vector<std::string> pRefId2name)
{

    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn1("R1.sam");
    seqan::BamStream bamStreamIn2("R2.sam");
    // std::vector<seqan::BamAlignmentRecord> buffer_mate1;
    // std::vector<seqan::BamAlignmentRecord> buffer_mate2;
    std::vector<seqan::BamAlignmentRecord> mate_bins;
    std::vector<int> pSharedBinIntvalTree;
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
    // unordered_map<string, int> pDanglingSequences ;
    int pBinSize = 10000;
    int pResultIndex;
    int start = 0;
    int end;

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
    // while (j < pNumberOfItemsPerBuffer)
    // {
    //     readRecord(record1, bamStreamIn1);
    //     readRecord(record2, bamStreamIn2);
    //     j++;
    // }
    all_data_read = 1;
    iter_num += 1;
    std::cout << __LINE__ << std::endl;

    while (!atEnd(bamStreamIn1) && !atEnd(bamStreamIn2))
    {
        readRecord(record1, bamStreamIn1);
        while ((hasFlagAllProper(record1) == 1) & 256 == 256)
        {
            readRecord(record1, bamStreamIn1);
        }
        all_data_read = 1;
        break;
        readRecord(record2, bamStreamIn2);
        while ((hasFlagAllProper(record2) == 1) & 256 == 256)
        {
            readRecord(record2, bamStreamIn2);
        }
    }
    std::cout << __LINE__ << std::endl;

    assert(record1.qName == record2.qName && "Be sure that the sam files have the same read order ");
    // if any of the reads is not mapped
    if (((hasFlagAllProper(record1) == 1) && 0x4 == 4) || ((hasFlagAllProper(record2) == 1) && 0x4 == 4))
    {
        one_mate_unmapped += 1;
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
        // maybe this->
        if (is_duplicated(record1.qName, record1.beginPos, // chromosome_refName
                          record2.qName, record2.beginPos))
        {
            duplicated_pairs += 1;
        }
    }
    std::cout << __LINE__ << std::endl;
    std::cout << duplicated_pairs << std::endl;


    // buffer_mate1.push_back(record1);
    // buffer_mate2.push_back(record2);
    // j += 1;
    // if (((all_data_read != 0) && (getAlignmentLengthInRef(record1) != 0)) && getAlignmentLengthInRef(record2) != 0)
    // {
    //     return buffer_mate1, buffer_mate2, 1, duplicated_pairs,
    //            one_mate_unmapped, one_mate_not_unique, one_mate_low_quality,
    //            iter_num - buffer_mate1.size();
    // }
    // if ((all_data_read == 0 && buffer_mate1.size() == 0) || (buffer_mate2.size() == 0))
    // {
    //     return 0, 0, 1, duplicated_pairs,
    //            one_mate_unmapped, one_mate_not_unique, one_mate_low_quality,
    //            iter_num - buffer_mate1.size();
    // }
    // else
    //     return buffer_mate1, buffer_mate2, 0, duplicated_pairs, one_mate_unmapped,
    //            one_mate_not_unique, one_mate_low_quality, iter_num - buffer_mate1.size();
    // while (!atEnd(bamStreamIn1) && !atEnd(bamStreamIn2))
    // {
    //     readRecord(record1, bamStreamIn1);
    //     seqan::CharString mate_ref1 = record1.qName; //get<record1.rname>(pRefId2name);
    //     readRecord(record2, bamStreamIn2);
    //     seqan::CharString mate_ref2 = record2.qName; //get<record2.rname>(pRefId2name);
    //     int read_middle = record1.beginPos + int(seqan::BamAlignmentRecord::getContigLength(record1, bamStreamIn1) / 2);
    //     int middle_pos = int((start + end) / 2);

    //     std::vector<int> middle_position_element;
    //     std::string x = std::pSharedBinIntvalTree[middle_pos];
    //     for (int i = 0; i < middle_position_element.size(); i++)
    //     {
    //         for (int j = 0; j < pSharedBinIntvalTree.size(); j++)
    //         {
    //             if (pSharedBinIntvalTree[j] == x)
    //             {
    //                 middle_position_element.insert(x);
    //             }
    //         }
    //     }
    //     while (!start > end)
    //     {
    //         if ((middle_position_element.begin() <= read_middle) && (read_middle <= middle_position_element.end()))
    //         {
    //             std::string mate_bin = pSharedBinIntvalTree[middle_pos];
    //             mate_is_unasigned = 0;
    //         }
    //         else if (middle_position_element.begin() <= read_middle)
    //         {
    //             end = middle_pos - 1;
    //             middle_pos = int((start + end) / 2);
    //             mate_is_unasigned = 1;
    //         }
    //         else
    //         {
    //             start = middle_pos + 1;
    //             middle_pos = int((start + end) / 2);
    //             mate_is_unasigned = 1;
    //         }
    //     }
    //     try
    //     {
    //         throw mate_is_unasigned = 1;
    //     }
    //     catch (...)
    //     {
    //         std::cout << "An exception occurred.";
    //     }
    //     if (mate_bins.empty())
    //     {
    //         mate_is_unasigned = 1;
    //         break;
    //     }

    return 0;
}

std::unordered_map<std::string, IntervalTree<size_t, size_t>> HiCBuildMatrix::createInitialStructures(seqan::BamStream pBamStream, int pBinSize)
{

    // if (!isGood(pBamStream))
    // {
    //     std::cerr << "ERROR: Could not open the first file!\n";
    //     return NULL;
    // }

    std::vector<std::string> chromosome_refName; // vector with chromosomes
    std::vector<size_t> chromosome_size;         // vector with sizes

    for (size_t i = 0; i < length(pBamStream.header.sequenceInfos); ++i)
    {

        chromosome_refName.push_back(seqan::toCString(pBamStream.header.sequenceInfos[i].i1));

        chromosome_size.push_back(pBamStream.header.sequenceInfos[i].i2);
        //}
    }
    for (size_t j = 0; j < chromosome_refName.size(); j++) // read names of chromosomes
    {
        std::cout << chromosome_refName[j] << ' ' << chromosome_size[j] << std::endl;
        //   } // read sizes of chromosomes
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
