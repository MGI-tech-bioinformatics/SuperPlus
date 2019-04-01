// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
//
// Parse special barcoded fastq files.  These are not regulation fastqs.
//
// FASTQ lines:
// - read1, qual1, read2, qual2, barcode, barcode_qual, sample index, sample index qual
// - no plus (+) signs.  Line ending signifies the end of each of the above entries in a
//   record.
// - PHRED offset TBD (+33, +64)

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "feudal/BinaryStream.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/VirtualMasterVec.h"
#include "10X/Barcode.h"
#include <fstream>
#include <string>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <list>
#include <zlib.h>

namespace {

struct PairError : std::runtime_error {
     PairError(String const& s) : std::runtime_error(s) {};
     PairError(string const& s) : std::runtime_error(s) {};
     PairError(char const * s) : std::runtime_error(s) {};
};
struct ParseError : std::runtime_error {
     ParseError(String const& s) : std::runtime_error(s) {};
     ParseError(string const& s) : std::runtime_error(s) {};
     ParseError(char const * s) : std::runtime_error(s) {};
};
class EndOfFile {};

void convertPhred( String const& buf, qualvector& qtmp, const uchar phred_offset=33 )
{
     qtmp.clear();
     for ( char const c : buf  ) {
          if ( c == '\n' || c == '\r' ) continue;
          qtmp.push_back( c - phred_offset );
     }
}


bool hasGemGroup( String const& s )
{
     // presence of a dash -- searching from end for efficiency
     auto itr=s.end();
     while ( itr-- != s.begin() ) if ( *itr == '-' ) return true;
     return false;
}

// parse standard fastq files to extract barcodes
void countBarcode(const std::pair<String, String> &fastq_pair, std::unordered_set<barcode_t> &bc_set) {

    // count bc number based on fastq for read pair 1
    auto &file = fastq_pair.first;
    if (!(file.EndsWith(".fastq.gz") || file.EndsWith(".fq.gz") || file.EndsWith(".gz"))) {
        std::cerr << "read pair fastq input has to be in gz format" << endl;
        Scram(1);
    }
    gzFile fastq_in_fp = gzopen(file.c_str(), "rb");
    int buf_len = 1024;
    char buf_c_str[buf_len];
    /*
    String command;
    if (file.EndsWith(".gz")) {
        command = "zcat " + file;
    } else {
        command = "cat " + file;
    }

    fast_pipe_ifstream input( command );
    */
    String buf;
    buf.reserve(1024 * 1024);
    std::size_t read_cnt = 0;
    while ( 1 ) {
        try {
            // name and barcode
            Bool fail = !gzgets(fastq_in_fp, buf_c_str, buf_len);
            if (fail) throw EndOfFile();
            buf = buf_c_str;
            if ( !buf.StartsWith("@") ) throw ParseError(buf);
            buf = buf.SafeAfterLast("#");
            buf = buf.SafeBefore("\t");
            bc_set.emplace(createBarcodeIndex(buf));
            ++read_cnt;
            
            // skip the next 3 lines
            gzgets(fastq_in_fp, buf_c_str, buf_len);
            gzgets(fastq_in_fp, buf_c_str, buf_len);
            gzgets(fastq_in_fp, buf_c_str, buf_len);
        } catch ( ParseError const& e ) {
            FatalErr( "out of sync reading line: " + string(e.what()) );
        } catch ( EndOfFile const& e ) {
            break;
        }
    }
    gzclose(fastq_in_fp);
    std::cerr << "total reads: " << (read_cnt * 2)<< std::endl;
    std::cerr << "total barcodes: " << bc_set.size() << std::endl;
}


void newUnpackBarcodeSortedFastq( const String& file, vecbasevector& b0_bases,
          VecPQVec& b0_quals, vecbasevector& bc_bases, VecPQVec& bc_quals,
          vec<uint32_t>& bcs )
{
     // TODO: work for non-gzipped
     ForceAssert(file.EndsWith(".fastq.gz") || file.EndsWith(".fasth.gz"));

     cout << Date() << ": " << file << endl;

     String buf;
     fast_pipe_ifstream input( "zcat " + file );
     enum { START, READ1, QUAL1, READ2, QUAL2, BARC, QUALBARC, INDEX, QUALINDEX, END };
     array<basevector,2> btmp;
     array<qualvector,2> qtmp;
     uint32_t barc = bcs.size() ? bcs.back() : 0u;
     String lastb;
     size_t line = 0;
     std::map<String,size_t> seen;

     while ( 1 ) {
          try {
               // name
               getline(input, buf); line++;
               if ( input.fail() ) throw EndOfFile();
               if ( !buf.StartsWith("@") ) throw ParseError(buf);

               // read a block
               for ( size_t blocki = START+1; blocki < END; ++blocki  ) {
                    getline(input, buf); line++;
                    if ( input.fail() ) throw ParseError(buf);

                    for ( auto& c : buf )
                         if ( c == 'n' || c == 'N' ) c = 'A';

                    switch( blocki ) {
                         case READ1:
                              btmp[0].SetFromStringWithNs( buf );
                              break;
                         case READ2:
                              btmp[1].SetFromStringWithNs( buf );
                              break;
                         case QUAL1:
                              convertPhred( buf, qtmp[0] );
                              break;
                         case QUAL2:
                              convertPhred( buf, qtmp[1] );
                              break;
                         case BARC:
                              // lack of a gem group indicates failing barcode
                              // presence of a gem group, but lack of a barcode indicates no barcode read
                              // AND no whitelist.  Treat them the same.
                              if ( hasGemGroup( buf ) && buf[0] != '-' ) {
                                   buf = buf.SafeBefore(",");              // new FASTH format has un-corrected barcode after the comma
                                   if ( lastb != buf ) {    // should trigger 1st time
                                        barc++;             // ...so 1st bc is 1
                                        lastb = buf;
                                        if ( seen.count(buf) > 0 ) {
                                             cout  << "barcode " << buf << " at line " << line << " already seen at line "  << seen[buf] << endl;
                                        }
                                        seen[buf]=line;
//                                        PRINT3( barc, seen.size(), buf );
                                   }
                                   bcs.push_back(barc);
                                   bcs.push_back(barc);
                                   bc_bases.push_back( btmp[0] );
                                   bc_bases.push_back( btmp[1] );
                                   bc_quals.push_back( PQVec( qtmp[0] ) );
                                   bc_quals.push_back( PQVec( qtmp[1] ) );
                              } else {
                                   b0_bases.push_back( btmp[0] );
                                   b0_bases.push_back( btmp[1] );
                                   b0_quals.push_back( PQVec( qtmp[0] ) );
                                   b0_quals.push_back( PQVec( qtmp[1] ) );
                              }

                              break;
                         case QUALBARC:
                         case INDEX:
                         case QUALINDEX:
                         default:
                              break;
                    }
               }

          } catch ( ParseError const& e ) {
               FatalErr( "out of sync reading line: " + string(e.what()) );
          } catch ( EndOfFile const& e ) {
               break;
          }
     }
}

void readIndirect( String const& filename, vec<String>& v )
{
     Ifstream(INPUT, filename);
     String buf;
     while (1)  {
          getline(INPUT, buf);
          if ( INPUT.fail() ) break;
          if ( buf.size() > 0 ) {
               buf.resize( buf.size()-1 );
               v.push_back(buf);
          }
     }
}

void mergeBarcodedReadFiles( vec<String> const& MERGE_HEADS, String const& OUT_HEAD )
{
     // parse input, including any indirect files
     vec<String> merge_heads;
     for ( auto const& merge_head : MERGE_HEADS ) {
          ForceAssertGt(merge_head.size(), 0u);
          if ( merge_head[0] == '@' ) readIndirect( merge_head.After("@"), merge_heads );
          else merge_heads.push_back( merge_head );
     }

     std::cerr << OUT_HEAD << std::endl;
     // create output files
     BinaryIteratingWriter<vec<int64_t>>     bci_out( OUT_HEAD + ".bci" );
     IncrementalWriter<basevector>           bases_out( OUT_HEAD + ".fastb" );
     IncrementalWriter<PQVec>                quals_out( OUT_HEAD + ".qualp" );

     // two passes -- do zero bc first because they must be aggregated.
     // all other files are assumed to have unique barcodes
     size_t count = 0;
     uint32_t bc_bias = 0;
     bci_out.write(0u);

     for (unsigned pass = 0; pass <=1; ++pass ) {
          for ( auto &merge_head : merge_heads ) {

               cout << Date() << ": pass " << pass+1 << " (of 2): " << merge_head << endl;

               VirtualMasterVec<basevector> bases( merge_head + ".fastb" );
               VirtualMasterVec<PQVec> quals( merge_head + ".qualp" );
               vec<int64_t> bci; BinaryReader::readFile( merge_head + ".bci", &bci );

               ForceAssertEq( bases.size(), quals.size() );
               ForceAssertEq( bases.size(), bci.back() );

               if ( bases.size() == 0 ) continue;

               for ( size_t i = 1; i < bci.size(); ++i )
                    if ( bci[i] < bci[i-1] ) {
                         cout << "barcode indices are out of order" << endl;
                         PRINT2(bci[i], bci[i-1]);
                    }

               if ( pass == 0 ) {
                    size_t this_count = bci[1];
                    for ( size_t i = 0; i < this_count; ++i ) {
                         bases_out.add( bases[i] );
                         quals_out.add( quals[i] );
                    }

                    count += this_count;
               } else {
                    size_t zero_count = bci[1];

                    // increment all of the offsets
                    ForceAssertGe( count, zero_count );

                    for ( size_t i = 1; i < bci.size(); ++i )
                         bci[i] = bci[i] - zero_count + count;

                    // write all reads and barcodes past bc0
                    for ( size_t i = zero_count; i < bases.size(); ++i )  {
                         bases_out.add( bases[i] );
                         quals_out.add( quals[i] );
                    }

                    for ( size_t i=1; i < bci.size()-1; ++i )
                         bci_out.write( bci[i] );

                    count += bases.size() - zero_count;
               }
          }

          if ( pass == 1 ) bci_out.write( count );
     }
     // delete all tmp files
     for ( auto &merge_head : merge_heads ) {
         Remove(merge_head + ".fastb");
         Remove(merge_head + ".qualp");
         Remove(merge_head + ".bci");
     }
}

void bucketAndSortFastqByBarcode(const std::pair<String, String> &fastq_pair,
                                 const String &OUT_HEAD,
                                 const std::unordered_set<barcode_t> &bc_set,
                                 std::size_t reads_per_bc = 0,
                                 std::size_t bucket_num = 256) {
    ForceAssert((bucket_num > 0) || (bucket_num <= 256));
    // put barcode into different buckets
    vec<vec<barcode_t>> bc_bucket_vec;
    if (bucket_num > bc_set.size()) {
        bucket_num = bc_set.size();
    }
    std::size_t bucket_size = bc_set.size() / bucket_num;
    vec<barcode_t> bc_vec;
    bc_vec.reserve(bucket_size);
    for (auto it = bc_set.cbegin(); it != bc_set.cend(); ++it) {
        if (bc_vec.size() < bucket_size) {
            bc_vec.emplace_back(*it);
        } else {
            bc_bucket_vec.emplace_back(bc_vec);
            bc_vec.assign(1, *it);
        }
    }
    // deal last part of barcode
    if (!bc_vec.empty()) {
        auto &last_vec = bc_bucket_vec.back();
        last_vec.insert(last_vec.end(), bc_vec.begin(), bc_vec.end());
    }
    vec<barcode_t> ().swap(bc_vec);
    
    /*
    for (auto &bc : bc_bucket_vec.front()) {
        std::cerr << bc << std::endl;
    }
    */

    // start to parse pair fastq file
    std::cerr << "There are " << bc_bucket_vec.size() << " buckets in total" << std::endl;
    Mkpath(OUT_HEAD.SafeBeforeLast("/"));
    // to record each tmp file name
    vec<String> tmp_out_heads;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < bc_bucket_vec.isize(); ++i) {
        auto &bc_vec = bc_bucket_vec.at(i);
        std::set<barcode_t> bc_set_partial(bc_vec.begin(), bc_vec.end());
        auto &file1 = fastq_pair.first, &file2 = fastq_pair.second;
        if (!(file1.EndsWith(".fastq.gz") || file1.EndsWith(".fq.gz") || file1.EndsWith(".gz"))) {
            std::cerr << "read pair fastq input has to be in gz format" << endl;
            Scram(1);
        }
        if (!(file2.EndsWith(".fastq.gz") || file2.EndsWith(".fq.gz") || file2.EndsWith(".gz"))) {
            std::cerr << "read pair fastq input has to be in gz format" << endl;
            Scram(1);
        }
        gzFile fastq_in_fp1 = gzopen(file1.c_str(), "rb");
        gzFile fastq_in_fp2 = gzopen(file2.c_str(), "rb");
        int buf_len = 1024;
        char buf_c_str1[buf_len], buf_c_str2[buf_len];
        String buf1, buf2;
        buf1.reserve(1024 * 1024);
        buf2.reserve(1024 * 1024);
        // std::size_t line1 = 0, line2 = 0;
        // pair <barcode, pair<bases, qual>>
        std::unordered_map<barcode_t, std::list<pair<basevector, basevector>>> barcode_base_map;
        std::unordered_map<barcode_t, std::list<pair<PQVec, PQVec>>> barcode_qual_map;
        qualvector qtmp1, qtmp2;
        while ( 1 ) {
            try {
                // name and barcode
                Bool fail1 = !gzgets(fastq_in_fp1, buf_c_str1, buf_len);
                Bool fail2 = !gzgets(fastq_in_fp2, buf_c_str2, buf_len);
                if (fail1 && fail2) throw EndOfFile();
                if ((fail1 && !fail2) || (!fail1 && fail2)) {
                    throw  PairError(file1 + " or " + file2);
                }
                buf1 = buf_c_str1;
                buf2 = buf_c_str2;
                if (!buf1.StartsWith("@")) throw ParseError(file1 + ": " + buf1);
                if (!buf2.StartsWith("@")) throw ParseError(file2 + ": " + buf2);

                buf1 = buf1.SafeAfterLast("#"); buf2 = buf2.SafeAfterLast("#");
                buf1 = buf1.SafeBefore("/"); buf2 = buf2.SafeBefore("/");
                barcode_t barcode_int1 = createBarcodeIndex(buf1);
                barcode_t barcode_int2 = createBarcodeIndex(buf2);
                if (barcode_int1 != barcode_int2) {
                    throw PairError(file1 + " or " + file2);
                }
                if (bc_set_partial.find(barcode_int1) == bc_set_partial.end()) {
                    gzgets(fastq_in_fp1, buf_c_str1, buf_len);
                    gzgets(fastq_in_fp1, buf_c_str1, buf_len);
                    gzgets(fastq_in_fp1, buf_c_str1, buf_len);
                    gzgets(fastq_in_fp2, buf_c_str2, buf_len);
                    gzgets(fastq_in_fp2, buf_c_str2, buf_len);
                    gzgets(fastq_in_fp2, buf_c_str2, buf_len);
                    continue;
                }

                // read 1 and 2
                gzgets(fastq_in_fp1, buf_c_str1, buf_len);
                gzgets(fastq_in_fp2, buf_c_str2, buf_len);
                buf1 = String(buf_c_str1).SafeBefore("\n");
                buf2 = String(buf_c_str2).SafeBefore("\n");
                for ( auto& c : buf1 ) {
                    if ( c == 'n' || c == 'N' ) c = 'A';
                }
                for ( auto& c : buf2 ) {
                    if ( c == 'n' || c == 'N' ) c = 'A';
                }
                auto base_pair = std::make_pair(basevector(buf1), basevector(buf2));
                // skip line '+'
                gzgets(fastq_in_fp1, buf_c_str1, buf_len);
                gzgets(fastq_in_fp2, buf_c_str2, buf_len);
                // quals
                gzgets(fastq_in_fp1, buf_c_str1, buf_len);
                gzgets(fastq_in_fp2, buf_c_str2, buf_len);
                buf1 = String(buf_c_str1).SafeBefore("\n");
                buf2 = String(buf_c_str2).SafeBefore("\n");
                convertPhred( buf1, qtmp1 );
                convertPhred( buf2, qtmp2 );
                auto &base_pairs = barcode_base_map[barcode_int1];
                auto &qual_pairs = barcode_qual_map[barcode_int1];
                // std::cerr << buf1 << std::endl;
                // std::cerr << buf2 << std::endl;
                // switch order if necessary
                if (barcode_int1 == 0) {
                    base_pairs.emplace_back(base_pair);
                    qual_pairs.emplace_back(PQVec(qtmp1), PQVec(qtmp2));
                    continue;
                }
                Bool is_insert = false;
                auto it_qual = qual_pairs.begin();
                for (auto it_base = base_pairs.begin();
                     it_base != base_pairs.end();
                     ++it_base, ++it_qual ) {
                    if (base_pair > *it_base) {
                        base_pairs.emplace(it_base, base_pair);
                        qual_pairs.emplace(it_qual, PQVec(qtmp1), PQVec(qtmp2));
                        is_insert = true;
                        break;
                    }
                }
                if (!is_insert) {
                    base_pairs.emplace_back(base_pair);
                    qual_pairs.emplace_back(PQVec(qtmp1), PQVec(qtmp2));
                }
                /*
                for (int j = base_pairs.isize() - 1; j >= 1; --j) {
                    auto &base_rhs = base_pairs.at(j);
                    auto &base_lhs = base_pairs.at(j - 1);
                    if (base_rhs > base_lhs) {
                        base_lhs.swap(base_rhs);
                        // switch qual too
                        auto &qual_rhs = qual_pairs.at(j);
                        auto &qual_lhs = qual_pairs.at(j - 1);
                        qual_lhs.swap(qual_rhs);
                    }
                }
                */
            } catch ( ParseError const& e ) {
                FatalErr( "out of sync reading line: " + string(e.what()) );
            } catch (PairError const & e) { 
                FatalErr( "something not match with pair file: " + string(e.what()) );
            } catch ( EndOfFile const& e ) {
                break;
            }
        }
        gzclose(fastq_in_fp1);
        gzclose(fastq_in_fp2);
        // write to tmp output files
        String tmp_file_prefix = OUT_HEAD + "_tmp_" + ToString(i);
        std::cerr << "Bucket " << i << " finished. Writing out tmp file: "
                  << tmp_file_prefix << ".{fastb, qualp, bci}" << std::endl;
        #pragma omp critical(file_prefix)
        {
            tmp_out_heads.emplace_back(tmp_file_prefix);
        }
        // create output files
        BinaryIteratingWriter<vec<int64_t>>     bci_out( tmp_file_prefix + ".bci" );
        IncrementalWriter<basevector>           bases_out( tmp_file_prefix + ".fastb" );
        IncrementalWriter<PQVec>                quals_out( tmp_file_prefix + ".qualp" );
        // first element in
        std::size_t bc_cnt = 0;
        bci_out.write(bc_cnt);

        // check if barcode 0 in this bucket
        if (barcode_base_map.find(0) != barcode_base_map.end()) {  // find barcode 0
            auto &base_pairs = barcode_base_map.at(0);
            auto &qual_pairs = barcode_qual_map.at(0);
            auto it_qual = qual_pairs.begin();
            for (auto it_base = base_pairs.cbegin();
                 it_base != base_pairs.cend();
                 ++it_base, ++it_qual ) {
                bases_out.add(it_base->first);
                bases_out.add(it_base->second);
                quals_out.add(it_qual->first);
                quals_out.add(it_qual->second);
            }

            bc_cnt += base_pairs.size() * 2;
        }
        bci_out.write(bc_cnt);
        for (auto it_bc = bc_set_partial.cbegin();
             it_bc != bc_set_partial.cend();
             ++it_bc) {
            auto &barcode_int = *it_bc;
            if (barcode_int == 0) {
                continue;
            }
            auto &base_pairs = barcode_base_map.at(barcode_int);
            if (reads_per_bc && ((base_pairs.size() * 2) >= reads_per_bc)) {
                continue;
            }
            auto &qual_pairs = barcode_qual_map.at(barcode_int);
            auto it_qual = qual_pairs.begin();
            for (auto it_base = base_pairs.cbegin();
                    it_base != base_pairs.cend();
                    ++it_base, ++it_qual ) {
                bases_out.add(it_base->first);
                bases_out.add(it_base->second);
                quals_out.add(it_qual->first);
                quals_out.add(it_qual->second);
            }
            bc_cnt += base_pairs.size() * 2;
            bci_out.write(bc_cnt);
        }

        std::unordered_map<barcode_t, std::list<pair<basevector, basevector>>> ().swap(barcode_base_map);
        std::unordered_map<barcode_t, std::list<pair<PQVec, PQVec>>> ().swap(barcode_qual_map);
        std::cerr << "finish writing tmp file for bucket " << i << std::endl;
    }  // finish bucket and sort fastq by barcode
    
    // merge tmp files
    std::cerr << "start to merge all tmp files" << std::endl;
    mergeBarcodedReadFiles(tmp_out_heads, OUT_HEAD);
}

};

int main(int argc, char *argv[])
{
     double clock = WallClockTime();

     RunTime( );

     BeginCommandArguments;
     CommandArgument_StringSet_OrDefault_Doc(FASTQS, "", "list of barcode-fastq files with reads\n"
     "will consider a read to be 'unbarcoded' if there is no appended gem group (e.g. BBBBBBBBBBBBBB-1)");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_BUCKETS, 256,                   
                         "Number of bucket. By default,  256");
     CommandArgument_UnsignedInt_OrDefault_Doc(READS_PER_BC, 0,                   
                         "Maximium number of reads per barcode.  By default, 0.");
     CommandArgument_String_Doc(OUT_HEAD, "basename of output files {.fastb, .qualp, .bc, .bci}");
     CommandArgument_StringSet_OrDefault_Doc(MERGE_HEADS, "",
               "list of heads for fastb/qualp/bc/bci files to merge, rather than produce" );
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,                   
                         "Number of threads.  By default, the number of processors online.");
     CommandArgument_Double_OrDefault_Doc(MAX_MEM_GB, 0,                        
             "if specified, maximum *suggested* RAM use in GB; in some cases may be "
             "exceeded by our code");             
     EndCommandArguments;
     SetThreads( NUM_THREADS, False );                                          
     SetMaxMemoryGBCheck(MAX_MEM_GB);  

     if (OUT_HEAD.EndsWith("/")) {
         std::cerr << "OUT_HEAD can not end with '/'" << std::endl;
         Scram(1);
         
     }
     // Define data structures.

     if ( MERGE_HEADS.size() ) {

          // MERGE path
          if ( FASTQS.size() ) FatalErr( "must not supply FASTQ/FASTQ0 with MERGE_HEADS" );

          mergeBarcodedReadFiles( MERGE_HEADS, OUT_HEAD );

     } else {
          // FASTQ->FASTB/QUALP path

          if ( FASTQS.size() != 2) FatalErr( ToString(FASTQS[0]) );
          std::pair<String, String> fastq_pair(FASTQS[0], FASTQS[1]);
          std::unordered_set<barcode_t> bc_set;
          countBarcode(fastq_pair, bc_set);
          bucketAndSortFastqByBarcode(fastq_pair, OUT_HEAD, bc_set, READS_PER_BC, NUM_BUCKETS);
#if 0
          // FASTQ->FASTB/QUALP path

          if ( !FASTQS.size() ) FatalErr( "no inputs specified FASTQ or FASTQ0" );

          vec<String> fastqs;

          for ( auto const& file : FASTQS )
               if ( file[0] == '@' ) readIndirect(file.After("@"), fastqs );
               else fastqs.push_back(file);


          // we're assuming that barcodes are not split across files
          vecbasevector bc_bases,b0_bases;
          VecPQVec bc_quals,b0_quals;
          vec<uint32_t> bcs;
          for ( auto const& file : FASTQS )
               newUnpackBarcodeSortedFastq( file, b0_bases, b0_quals, bc_bases, bc_quals, bcs );


          vec<int64_t> bci;
          bci.push_back(0);
          int64_t b0_size = b0_bases.size();
          int64_t cur = 0;
          for (size_t i = 0; i < bcs.size(); ++i )
               if ( bcs[i] != cur ) {        // bcs should start != 0, so we cover the end of 0 here
                    cur = bcs[i];
                    bci.push_back(i+b0_size);
               }
          bci.push_back(b0_size+bcs.size());


          cout << Date() << ": writing output to " + OUT_HEAD + " .fastb,.qualp,.bci " << endl;
          IncrementalWriter<basevector> bases_out( OUT_HEAD + ".fastb" );
          IncrementalWriter<PQVec>      quals_out( OUT_HEAD + ".qualp" );
          bases_out.add( b0_bases.begin(), b0_bases.end() );
          bases_out.add( bc_bases.begin(), bc_bases.end() );
          quals_out.add( b0_quals.begin(), b0_quals.end() );
          quals_out.add( bc_quals.begin(), bc_quals.end() );
          BinaryWriter::writeFile( OUT_HEAD + ".bci", bci ) ;

          // barcode stats
          size_t count = 0;
          for ( size_t i = 1; i < bci.size()-1; ++i )  {
               if ( bci[i+1] - bci[i] > 10 ) count++;
          }
          cout << Date() << ": " << count << " barcodes have more than 10 reads, out of a total of "
               << bci.size()-2 << " barcodes" << endl;

          // Done.
#endif
     }

     cout << TimeSince(clock) << " used, peak mem = " << PeakMemUsageGBString( )
          << endl;
     cout << Date( ) << ": done" << endl;
     Scram(0);
}
