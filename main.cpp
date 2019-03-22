#include <cstdlib>
#include <cmath>
#include <errno.h>

#include <memory>
#include <iostream>

#include <bitset>

#include "range_coder.h"
#include "frequency_model.h"
#include "vcf_reader.h"

#include <zstd.h>
#include <zstd_errors.h>

// #define bcf_gt_phased(idx) (((idx)+1)<<1|1)
// phase = & 1
// (idx-1) >> 1

// (allele+ 1)<<1|phased


const uint8_t BCF_UNPACK_TOMAHAWK[3] = {2, 0, 1};
const uint8_t BCF_UNPACK_TOMAHAWK_GENERAL[16] = {15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
#define BCF_UNPACK_GENOTYPE(A) BCF_UNPACK_TOMAHAWK[(A >> 1)]
#define BCF_UNPACK_GENOTYPE_GENERAL(A) BCF_UNPACK_TOMAHAWK_GENERAL[(A >> 1)]

template <int n_symbols = 2>
class PBWT {
public:
    //PBWT() : n_samples(0), prev(nullptr), ppa(nullptr), n_q1(0), n_q2(0), q1(nullptr), q2(nullptr){}
    PBWT() :
        n_samples(0), n_steps(0),
        prev(nullptr), ppa(nullptr)
    {
        static_assert(n_symbols > 1, "Number of symbols must be > 1");
        for(int i = 0; i < n_symbols; ++i) {
            queue[i] = nullptr;
        }
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);
    }

    PBWT(int64_t n_samples) :
        n_samples(n_samples),
        n_steps(0),
        prev(new uint8_t[n_samples]),
        ppa(new uint32_t[n_samples])
    {
        static_assert(n_symbols > 1, "Number of symbols must be > 1");
        for(int i = 0; i < n_symbols; ++i) {
            queue[i] = new uint32_t[n_samples];
        }
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);
        reset();
    }

    ~PBWT() {
        delete[] prev;
        delete[] ppa;
        for(int i = 0; i < n_symbols; ++i) {
            delete[] queue[i];
        }
    }

    void Initiate(int64_t n_samples) {
        delete[] prev;
        delete[] ppa;
        for(int i = 0; i < n_symbols; ++i) {
            delete[] queue[i];
        }
        this->n_samples = n_samples;
        prev = new uint8_t[n_samples];
        ppa  = new uint32_t[n_samples];
        for(int i = 0; i < n_symbols; ++i) {
            queue[i] = new uint32_t[n_samples];
        }
        reset();
    }

    void reset() {
        if(ppa != nullptr && n_samples) {
            for(int i = 0; i < n_samples; ++i) {
                ppa[i] = i;
            }
        }
        if(prev != nullptr && n_samples) {
            memset(prev, 0, sizeof(uint8_t)*n_samples);
        }
        //n_q1 = n_q2 = 0;
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);
    }

    int Update(const int* arr) {
        //n_q1 = n_q2 = 0;
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

        for(int i = 0; i < n_samples; ++i) {
            const uint8_t& gt = BCF_UNPACK_GENOTYPE(arr[ppa[i]]);
            for(int j = 0; j < n_symbols; ++j) {
                if(gt == j)
                    queue[j][n_queue[j]++] = ppa[i];
            }
            prev[i] = gt;
            //std::cerr << " " << (int)gt;
        }
        //std::cerr << std::endl;

        uint32_t of = 0;
        for(int j = 0; j < n_symbols; ++j) {
            for(int i = 0; i < n_queue[j]; ++i, ++of) ppa[of] = queue[j][i];
        }
        //for(int i = 0; i < n_q2; ++i, ++of) ppa[of] = q2[i];
        assert(of == n_samples);
        // debug should be sorted here
        //for(int i = 0; i < n_samples; ++i) {
         //   std::cerr << " " << (int)BCF_UNPACK_GENOTYPE(arr[ppa[i]]);
        //}
        //std::cerr << std::endl;
        ++n_steps;

        return(1);
    }

    template <class T>
    int Update(const T* arr, uint32_t stride = 1) {
        //n_q1 = n_q2 = 0;
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

        for(int i = 0; i < n_samples; ++i) {
            const uint8_t& gt = BCF_UNPACK_GENOTYPE(arr[ppa[i] * stride]);
            assert(gt < n_symbols);
            for(int j = 0; j < n_symbols; ++j) {
                if(gt == j)
                    queue[j][n_queue[j]++] = ppa[i];
            }
            prev[i] = gt;
            //std::cerr << " " << (int)gt;
        }
        //assert(p == n_samples);
        //std::cerr << std::endl;

        uint32_t of = 0;
        for(int j = 0; j < n_symbols; ++j) {
            for(uint32_t i = 0; i < n_queue[j]; ++i, ++of)
                ppa[of] = queue[j][i];
        }
        //for(int i = 0; i < n_q2; ++i, ++of) ppa[of] = q2[i];
        //std::cerr << "of=" << of << "/" << n_samples << std::endl;
        assert(of == n_samples);
        // debug should be sorted here
        //for(int i = 0; i < n_samples; ++i) {
         //   std::cerr << " " << (int)BCF_UNPACK_GENOTYPE(arr[ppa[i]]);
        //}
        //std::cerr << std::endl;
        ++n_steps;

        return(1);
    }

    std::string ToPrettyString() const {
        std::string ret = "n=" + std::to_string(n_samples) + "{";
        ret += std::to_string(ppa[0]);
        for(int i = 1; i < n_samples; ++i) {
            ret += "," + std::to_string(ppa[i]);
        }
        ret += "}";
        return(ret);
    }

public:
    int64_t n_samples;
    uint64_t n_steps;
    uint8_t* prev;
    uint32_t* ppa;
    uint32_t n_queue[n_symbols];
    uint32_t* queue[n_symbols];
};

int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
}

int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int ret = ZSTD_decompress(out, out_capacity, in, n_in);
    return(ret);
}

void ReadVcfGT (const std::string& filename) {
    std::unique_ptr<tachyon::io::VcfReader> reader = tachyon::io::VcfReader::FromFile(filename,8);
    if(reader.get() == nullptr) {
        std::cerr << "failed read" << std::endl;
        exit(1);
    }


    //std::cerr << hr->n[0] << "," << hr->n[1] << "," << hr->n[2] << std::endl;
    //assert(reader->n_samples_ == 2504);
    //bcf1_t *rec = bcf_init1();

    //const uint32_t n_samples = hr->n[2];
    //int mgt_arr = 0;
    //int* gt_arr = nullptr;
   // int* gt1_a = new int[uint32_t(2.5*reader->n_samples_)];
    //int* gt2_a = new int[uint32_t(2.5*reader->n_samples_)];

    uint64_t n_out = 0, n_in = 0, n_data = 0, n_out1 = 0, n_out2 = 0, n_out_gtpbwt = 0;
    //uint8_t* data_buffer = new uint8_t[2*10000000];
    uint8_t* gtpbwt_buffer  = new uint8_t[reader->n_samples_];

    uint8_t* out1_buffer = new uint8_t[10000000];
    uint8_t* out2_buffer = new uint8_t[10000000];
    uint8_t* out_gtpbwt  = new uint8_t[10000000];
    uint8_t* out_gtpbwt_24  = new uint8_t[10000000];

    pil::RangeCoder rc1, rc2, rc_pbwt, rc_pbwt4, rc_pbwt24;
    pil::FrequencyModel<2>* fmodel1 = new pil::FrequencyModel<2>[4096]; // order-12 model
    pil::FrequencyModel<2>* fmodel2 = new pil::FrequencyModel<2>[4096]; // order-12 model
    //pil::FrequencyModel<2>* fmodel1_hard = new pil::FrequencyModel<2>[4096];
    //pil::FrequencyModel<2>* fmodel2_hard = new pil::FrequencyModel<2>[4096];

    std::cerr << "N-samples=" << reader->n_samples_ << std::endl;

#define MODEL_SIZE 16556
    pil::FrequencyModel<4>* gtpbwt_model = new pil::FrequencyModel<4>[MODEL_SIZE]; // order-7 model
    pil::FrequencyModel<16>* gtpbwt_model_4 = new pil::FrequencyModel<16>[MODEL_SIZE]; // order-2 model
    // if alt > 4 then use individual pbwts OR raw + zstd

    rc1.SetOutput(out1_buffer);
    rc2.SetOutput(out2_buffer);
    rc_pbwt.SetOutput(out_gtpbwt);
    rc_pbwt4.SetOutput(out_gtpbwt);
    rc1.StartEncode();
    rc2.StartEncode();
    rc_pbwt.StartEncode();
    rc_pbwt4.StartEncode();

    PBWT<2>* pbwt = new PBWT<2>[2];
    pbwt[0].Initiate(reader->n_samples_);
    pbwt[1].Initiate(reader->n_samples_);

    PBWT<4> gt_pbwt(reader->n_samples_);
    PBWT<16> gt_pbwt4(reader->n_samples_);

    uint32_t n_processed = 0;
    uint32_t n_lines = 0;
    uint32_t last_flush_pos = 0;

    uint32_t n_under_10 = 0;

    //uint32_t n_bitmaps = std::ceil((float)reader->n_samples_/64);
    //uint64_t* bitmaps = new uint64_t[n_bitmaps];

    while (reader->Next()) {
        //const char* chrom = bcf_seqname(hr,line) ;
        //if (!p->chrom) p->chrom = strdup (chrom) ;
        //else if (strcmp (chrom, p->chrom)) break ;
        //int pos = line->pos; // bcf coordinates are 0-based
        //char *ref, *REF;
        //ref = REF = strdup(line->d.allele[0]);
        //while ( (*ref = toupper(*ref)) ) ++ref ;
        if(reader->bcf1_->n_allele > 4) {
//            if(reader->bcf1_->n_allele <= 4) {
//
//            } else
            std::cerr << "skip: n_allele=" << reader->bcf1_->n_allele << std::endl;
            continue;
        }


//        pbwt[0].Update(reader->bcf1_->d.fmt[0].p+0, 2);
//        pbwt[1].Update(reader->bcf1_->d.fmt[0].p+1, 2);
        ++n_processed;
        ++n_lines;
        // encode

        n_in += 2*reader->n_samples_;

        if (reader->bcf1_->n_allele != 2) {
            // copy to 4-bit
            uint8_t* gts = reader->bcf1_->d.fmt[0].p;
            uint32_t j = 0;
            for(int i = 0; i < 2*reader->n_samples_; i+=2, ++j) {
                gtpbwt_buffer[j] = (BCF_UNPACK_GENOTYPE_GENERAL(gts[i+0]) << 2) |  BCF_UNPACK_GENOTYPE_GENERAL(gts[i+1]);
                //std::cerr << " " << std::bitset<4>(gtpbwt_buffer[j]);
            }
            //std::cerr << std::endl;
            gt_pbwt4.Update(gtpbwt_buffer, 1);

            //size_t before = rc_pbwt.OutSize();
            {
               uint32_t s = 0;
               for(int i = 0; i < j; ++i) {
                   gtpbwt_model_4[s].EncodeSymbol(&rc_pbwt4, gt_pbwt4.prev[i]);
                   s <<= 2;
                   s |= gt_pbwt4.prev[i];
                   s &= (MODEL_SIZE-1);
                   _mm_prefetch((const char *)&gtpbwt_model_4[s], _MM_HINT_T0);
               }
            }
            //std::cerr << rc_pbwt.OutSize() - before << " b for " << reader->bcf1_->n_allele << std::endl;
        }

        else
        {
            // copy to two-bit
            uint8_t* gts = reader->bcf1_->d.fmt[0].p;
            int j = 0;
            for(int i = 0; i < 2*reader->n_samples_; i+=2, ++j) {
                gtpbwt_buffer[j] = (BCF_UNPACK_GENOTYPE(gts[i+0]) << 1) |  BCF_UNPACK_GENOTYPE(gts[i+1]);
                //std::cerr << " " << std::bitset<2>(gtpbwt_buffer[j]);
            }
            //std::cerr << std::endl;
            gt_pbwt.Update(gtpbwt_buffer, 1);
            {
               uint32_t s = 0;
               for(int i = 0; i < reader->n_samples_; ++i) {
                   gtpbwt_model[s].EncodeSymbol(&rc_pbwt, gt_pbwt.prev[i]);
                   s <<= 2;
                   s |= gt_pbwt.prev[i];
                   s &= (MODEL_SIZE-1);
                   _mm_prefetch((const char *)&gtpbwt_model[s], _MM_HINT_T0);
               }
            }
        }

        //if (reader->bcf1_->pos - last_flush_pos > 500000) {
        if (n_lines == 8196 || rc_pbwt.OutSize() > 9000000) {
            n_out_gtpbwt += rc_pbwt.OutSize();
            n_out_gtpbwt += rc_pbwt4.OutSize();
            rc_pbwt.FinishEncode();
            rc_pbwt.SetOutput(out_gtpbwt);
            rc_pbwt.StartEncode();
            rc_pbwt4.FinishEncode();
            rc_pbwt4.SetOutput(out_gtpbwt);
            rc_pbwt4.StartEncode();
            std::cerr << "gtPBWT " << n_lines << " Compressed=" << n_in << "->" << n_out_gtpbwt << "(" << (float)n_in/n_out_gtpbwt << "-fold)" << std::endl;
            delete[] gtpbwt_model;
            gtpbwt_model = new pil::FrequencyModel<4>[MODEL_SIZE];
            delete[] gtpbwt_model_4;
            gtpbwt_model_4 = new pil::FrequencyModel<16>[MODEL_SIZE];
            gt_pbwt.reset();
            gt_pbwt4.reset();
            n_lines = 0;
            last_flush_pos = reader->bcf1_->pos;
        }


//        if(n_lines % 10000 == 0) {
//            //std::cerr << "PBWT1: " << pbwt[0].ToPrettyString() << std::endl;
//            //std::cerr << "PBWT2: " << pbwt[1].ToPrettyString() << std::endl;
//            std::cerr << "gtpbwt: " << gt_pbwt.ToPrettyString() << std::endl;
//        }


        // No PBWT order example
        /*
        {
            uint16_t s = 0;
            uint8_t* a = reader->bcf1_->d.fmt[0].p;
            for(int i = 0; i < 2*reader->n_samples_; i += 2) {
               const uint8_t& gt = BCF_UNPACK_GENOTYPE(a[i]);
               fmodel1[s].EncodeSymbol(&rc1, gt);
               s <<= 1;
               s |= (gt == 1);
               s &= 2047;
               _mm_prefetch((const char *)&fmodel1[s], _MM_HINT_T0);
               //data_buffer[n_data++] = pbwt[0].prev[i];
           }
        }

        {
            uint16_t s = 0;
            uint8_t* a = reader->bcf1_->d.fmt[0].p;
            for(int i = 0; i < 2*reader->n_samples_; i += 2) {
                const uint8_t& gt = BCF_UNPACK_GENOTYPE(a[i+1]);
               fmodel2[s].EncodeSymbol(&rc1, gt);
               s <<= 1;
               s |= (gt == 1);
               s &= 2047;
               _mm_prefetch((const char *)&fmodel2[s], _MM_HINT_T0);
               //data_buffer[n_data++] = pbwt[0].prev[i];
            }
        }
        */

//        if (pbwt[0].n_queue[1] >= 100) {
//           uint16_t s = 0;
//           for(int i = 0; i < reader->n_samples_; ++i) {
//               fmodel1_hard[s].EncodeSymbol(&rc1, pbwt[0].prev[i]);
//               s <<= 1;
//               s |= (pbwt[0].prev[i] == 1);
//               s &= 4095;
//               _mm_prefetch((const char *)&fmodel1_hard[s], _MM_HINT_T0);
//               //data_buffer[n_data++] = pbwt[0].prev[i];
//           }
//        }
//
//        else if(pbwt[0].n_queue[1] <= 10) {
//            uint32_t* out = reinterpret_cast<uint32_t*>(&out1_buffer[n_out1]);
//            uint32_t n = 0;
//            for(int i = 0; i < reader->n_samples_; ++i) {
//                if(pbwt[0].prev[i] == 1) {
//                    out[n++] = pbwt[0].ppa[i];
//                    //std::cerr << out[n-1] << std::endl;
//                }
//            }
//            //std::cerr << "ndata=" << n_data << std::endl;
//            n_out1 += sizeof(uint32_t)*pbwt[0].n_queue[1];
//        }

//        else

//        {
//            uint16_t s = 0;
//            for(int i = 0; i < reader->n_samples_; ++i) {
//                fmodel1[s].EncodeSymbol(&rc1, pbwt[0].prev[i]);
//                s <<= 1;
//                s |= (pbwt[0].prev[i] == 1);
//                s &= 4095;
//                _mm_prefetch((const char *)&fmodel1[s], _MM_HINT_T0);
//                //data_buffer[n_data++] = pbwt[0].prev[i];
//            }
//
//            n_under_10 += (pbwt[0].n_queue[1] <= 10);
//        }


//        if(pbwt[1].n_queue[1] >= 100) {
//           uint16_t s = 0;
//           for(int j = 0; j < reader->n_samples_; ++j) {
//               //std::cerr << (int)pbwt[1].prev[j] << " ";
//               fmodel2_hard[s].EncodeSymbol(&rc2, pbwt[1].prev[j]);
//               s <<= 1;
//               s |= (pbwt[1].prev[j] == 1);
//               s &= 4095;
//               _mm_prefetch((const char *)&fmodel2_hard[s], _MM_HINT_T0);
//               //data_buffer[n_data++] = pbwt[1].prev[j];
//           }
//        }
//        else if(pbwt[1].n_queue[1] <= 10) {
//            uint32_t* out = reinterpret_cast<uint32_t*>(&out2_buffer[n_out2]);
//            uint32_t n = 0;
//            for(int i = 0; i < reader->n_samples_; ++i) {
//                if(pbwt[1].prev[i] == 1) {
//                    out[n++] = pbwt[1].ppa[i];
//                    //std::cerr << out[n-1] << std::endl;
//                }
//            }
//            //std::cerr << "ndata=" << n_data << std::endl;
//            n_out2 += sizeof(uint32_t)*pbwt[1].n_queue[1];
//        }
//        else

//        {
//            uint16_t s = 0;
//            for(int j = 0; j < reader->n_samples_; ++j) {
//                //std::cerr << (int)pbwt[1].prev[j] << " ";
//                fmodel2[s].EncodeSymbol(&rc2, pbwt[1].prev[j]);
//                s <<= 1;
//                s |= (pbwt[1].prev[j] == 1);
//                s &= 4095;
//                _mm_prefetch((const char *)&fmodel2[s], _MM_HINT_T0);
//                //data_buffer[n_data++] = pbwt[1].prev[j];
//            }
//            //std::cerr << std::endl;
//            n_under_10 += (pbwt[1].n_queue[1] <= 10);
//        }

        //std::cerr << n_lines << "," << n_data << "/" << 10000000 << std::endl;
        //std::cerr << rc1.OutSize() << " + " << rc2.OutSize() << std::endl;
        //if(rc1.OutSize() > 1000000 || rc2.OutSize() > 1000000) {


        /*
        if(n_data > 9000000) {
           // std::cerr << "before zstd" << std::endl;
        //if(n_data > 9000000) {
            int ret = ZstdCompress(data_buffer,n_data,zstd_buffer,10000000,6);
            n_out += ret;
            n_out_zstd += ret;
            std::cerr << "Zstd " << n_lines << " Compressed=" << n_data << "->" << ret << "(" << (float)n_data/ret << "-fold)" << std::endl;
            //std::cerr << "Zstd Compressed=" << n_in << "->" << n_out_zstd << "(" << (float)n_in/n_out_zstd << "-fold)" << std::endl;
            n_data = 0;
        }
*/

        //std::cerr << n_lines << " " << rc1.OutSize() << " " << rc2.OutSize() << std::endl;
        /*if(n_lines == 8196 || rc1.OutSize() > 9000000 || rc2.OutSize() > 9000000 || n_out1 > 9000000 || n_out2 > 9000000) {
        //if(n_lines == 50000) {
//            int ret = ZstdCompress(out1_buffer,n_out1,zstd_buffer,10000000,6);
//            n_out += ret;
//            std::cerr << "Zstd " << n_lines << " Compressed=" << n_out1 << "->" << ret << "(" << (float)n_out1/ret << "-fold)" << std::endl;
//            n_out1 = 0;
//            ret = ZstdCompress(out2_buffer,n_out2,zstd_buffer,10000000,6);
//            n_out += ret;
//            std::cerr << "Zstd " << n_lines << " Compressed=" << n_out2 << "->" << ret << "(" << (float)n_out2/ret << "-fold)" << std::endl;
//            n_out2 = 0;

            n_out += rc1.OutSize();
            rc1.FinishEncode();
            rc1.SetOutput(out1_buffer);
            rc1.StartEncode();
            //std::cerr << "RC1 " << n_lines << " Compressed=" << n_in << "->" << n_out << "(" << (float)n_in/n_out << "-fold)" << std::endl;
            delete[] fmodel1;
//            delete[] fmodel1_hard;
            fmodel1 = new pil::FrequencyModel<2>[4096];
//            fmodel1_hard = new pil::FrequencyModel<2>[4096];
            pbwt[0].reset();

            n_out += rc2.OutSize();
            rc2.FinishEncode();
            rc2.SetOutput(out2_buffer);
            rc2.StartEncode();
            std::cerr << "RC " << n_lines << " Compressed=" << n_in << "->" << n_out << "(" << (float)n_in/n_out << "-fold)" << std::endl;
            delete[] fmodel2;
//            delete[] fmodel2_hard;
            fmodel2 = new pil::FrequencyModel<2>[4096];
//            fmodel2_hard = new pil::FrequencyModel<2>[4096];
            pbwt[1].reset();
            n_lines = 0;
        }
        */
    }


//    if(n_out1 || n_out2) {
//        int ret = ZstdCompress(out1_buffer,n_out1,zstd_buffer,10000000,6);
//        n_out += ret;
//        std::cerr << "Zstd " << n_lines << " Compressed=" << n_out1 << "->" << ret << "(" << (float)n_out1/ret << "-fold)" << std::endl;
//        n_out1 = 0;
//        ret = ZstdCompress(out2_buffer,n_out2,zstd_buffer,10000000,6);
//        n_out += ret;
//        std::cerr << "Zstd " << n_lines << " Compressed=" << n_out2 << "->" << ret << "(" << (float)n_out2/ret << "-fold)" << std::endl;
//        n_out2 = 0;
//    }

    /*if(rc1.OutSize() || rc2.OutSize()) {
        n_out += rc1.OutSize();
        rc1.FinishEncode();
        rc1.SetOutput(out1_buffer);
        rc1.StartEncode();
        //std::cerr << "RC1 Compressed=" << n_in << "->" << n_out << "(" << (float)n_in/n_out << "-fold)" << std::endl;

        n_out += rc2.OutSize();
        rc2.FinishEncode();
        rc2.SetOutput(out2_buffer);
        rc2.StartEncode();
        std::cerr << "RC Compressed=" << n_in << "->" << n_out << "(" << (float)n_in/n_out << "-fold)" << std::endl;
    }*/

    if (rc_pbwt.OutSize()) {
        n_out_gtpbwt += rc_pbwt.OutSize();
        rc_pbwt.FinishEncode();
        rc_pbwt.SetOutput(out_gtpbwt);
        rc_pbwt.StartEncode();
        std::cerr << "gtPBWT " << n_lines << " Compressed=" << n_in << "->" << n_out_gtpbwt << "(" << (float)n_in/n_out_gtpbwt << "-fold)" << std::endl;
    }

    std::cerr << "under 10= " << n_under_10 << " -> " << n_under_10*sizeof(uint32_t) << std::endl;

    delete[] fmodel1; delete[] fmodel2;
//    delete[] fmodel1_hard; delete[] fmodel2_hard;
    //delete[] gt_arr;
    //delete[] data_buffer;
    delete[] out1_buffer; delete[] out2_buffer;
    //delete[] zstd_buffer;
    //delete[] gt1_a; delete[] gt2_a;
    delete[] pbwt;
    delete[] gtpbwt_buffer;
    delete[] gtpbwt_model;
}

int main(int argc, char** argv) {
    std::cerr << argc << std::endl;
    if(argc == 1) return(1);
    else {
        std::cerr << std::string(argv[1]) << std::endl;
        ReadVcfGT(std::string(argv[1]));
        return(0);
    }
}
