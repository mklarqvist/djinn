#include <cstdlib>
#include <cmath>
#include <errno.h>

#include <memory>
#include <iostream>

#include "range_coder.h"
#include "frequency_model.h"

#include <zstd.h>
#include <zstd_errors.h>
#include <htslib/synced_bcf_reader.h>

// #define bcf_gt_phased(idx) (((idx)+1)<<1|1)
// phase = & 1
// (idx-1) >> 1

// (allele+ 1)<<1|phased


const uint8_t BCF_UNPACK_TOMAHAWK[3] = {2, 0, 1};
#define BCF_UNPACK_GENOTYPE(A) BCF_UNPACK_TOMAHAWK[(A >> 1)]

class PBWT {
public:
    //PBWT() : n_samples(0), prev(nullptr), ppa(nullptr), n_q1(0), n_q2(0), q1(nullptr), q2(nullptr){}
    PBWT() :
        n_samples(0), n_steps(0),
        prev(nullptr), ppa(nullptr),
        n_q1(0), n_q2(0),
        q1(nullptr), q2(nullptr)
    {

    }

    PBWT(int32_t n_samples) :
        n_samples(n_samples),
        n_steps(0),
        prev(new uint8_t[n_samples]),
        ppa(new uint32_t[n_samples]),
        n_q1(0), n_q2(0),
        q1(new uint32_t[n_samples]), q2(new uint32_t[n_samples])
    {
        reset();
    }

    virtual ~PBWT() {
        delete[] prev;
        delete[] ppa;
        delete[] q1;
        delete[] q2;
    }

    void Initiate(int32_t n_samples) {
        delete[] prev;
        delete[] ppa;
        delete[] q1;
        delete[] q2;
        this->n_samples = n_samples;
        prev = new uint8_t[n_samples];
        ppa  = new uint32_t[n_samples];
        q1   = new uint32_t[n_samples];
        q2   = new uint32_t[n_samples];
        reset();
    }

    void reset() {
        if(ppa != nullptr && n_samples) {
            for(int i = 0; i < n_samples; ++i) {
                ppa[i] = i;
            }
        }
        if(prev != nullptr && n_samples) {
            memset(prev, 0, sizeof(uint32_t)*n_samples);
        }
        n_q1 = n_q2 = 0;
    }

    virtual int Update(const int* arr) {
        n_q1 = n_q2 = 0;

        for(int i = 0; i < n_samples; ++i) {
            uint8_t gt = BCF_UNPACK_GENOTYPE(arr[ppa[i]]);
            switch(gt) {
            case(0): q1[n_q1++] = ppa[i]; break;
            case(1): q2[n_q2++] = ppa[i]; break;
            default: std::cerr << "error!!!" << std::endl; exit(1);
            }
            prev[ppa[i]] = gt;
            //std::cerr << " " << (int)gt;
        }
        //std::cerr << std::endl;

        uint32_t of = 0;
        for(int i = 0; i < n_q1; ++i, ++of) ppa[of] = q1[i];
        for(int i = 0; i < n_q2; ++i, ++of) ppa[of] = q2[i];
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
    int32_t n_samples;
    uint64_t n_steps;
    uint8_t* prev;
    uint32_t* ppa;
    uint32_t n_q1, n_q2;
    uint32_t *q1, *q2;
};

class BitmapPBWT : public PBWT {
public:
    BitmapPBWT(int32_t n_samples) : PBWT(n_samples), n_bitmaps(0), bitmaps(nullptr) {
        n_bitmaps = std::ceil(n_samples / 64.0);
        bitmaps = new uint64_t[n_bitmaps];
        memset(bitmaps,0,sizeof(uint64_t)*n_bitmaps);
    }

    ~BitmapPBWT() {
        delete[] bitmaps;
    }

    void resetBitmaps() {
        memset(bitmaps,0,sizeof(uint64_t)*n_bitmaps);
    }

    int Update(const int* arr) override {
        resetBitmaps();

        for(int i = 0; i < n_samples; ++i) {
            uint64_t gt = BCF_UNPACK_GENOTYPE(arr[ppa[i]]);
            switch(gt) {
            case(0): q1[n_q1++] = ppa[i]; break;
            case(1): q2[n_q2++] = ppa[i]; break;
            default: std::cerr << "error!!!" << std::endl; exit(1);
            }
            bitmaps[ppa[i] / 64] |= ((gt == 1) << (ppa[i] / 64 & 4));
        }

        uint32_t of = 0;
        for(int i = 0; i < n_q1; ++i, ++of) ppa[of] = q1[i];
        for(int i = 0; i < n_q2; ++i, ++of) ppa[of] = q2[i];
        assert(of == n_samples);
        n_q1 = n_q2 = 0;
        ++n_steps;

        return(1);
    }

public:
    int32_t n_bitmaps;
    uint64_t* bitmaps;
};

int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
}

void ReadVcfGT (const std::string& filename) {
    bcf_srs_t *sr = bcf_sr_init();
    if (!bcf_sr_add_reader (sr, filename.c_str())) {
        std::cerr << "reader dead" << std::endl;
        exit(1) ;
    }

    //htsFile *fp    = hts_open(filename, "rb");
    bcf_hdr_t *hr = sr->readers[0].header;
    //std::cerr << hr->n[0] << "," << hr->n[1] << "," << hr->n[2] << std::endl;
    assert(hr->n[2] == 2504);
    //bcf1_t *rec = bcf_init1();

    const uint32_t n_samples = hr->n[2];
    int mgt_arr = 0, *gt_arr = NULL;
    int* gt1_a = new int[uint32_t(2.5*n_samples)];
    int* gt2_a = new int[uint32_t(2.5*n_samples)];

    uint64_t n_out = 0, n_in = 0, n_data = 0, n_out_zstd = 0;
    uint8_t* data_buffer = new uint8_t[10000000];
    uint8_t* out1_buffer = new uint8_t[10000000];
    uint8_t* out2_buffer = new uint8_t[10000000];
    uint8_t* zstd_buffer = new uint8_t[10000000];
    pil::RangeCoder rc1, rc2;
    pil::FrequencyModel<2>* fmodel1 = new pil::FrequencyModel<2>[1024];
    pil::FrequencyModel<2>* fmodel2 = new pil::FrequencyModel<2>[1024];
    rc1.SetOutput(out1_buffer);
    rc1.StartEncode();
    rc2.SetOutput(out2_buffer);
    rc2.StartEncode();

    PBWT* pbwt = new PBWT[2];
    pbwt[0].Initiate(n_samples);
    pbwt[1].Initiate(n_samples);

    uint32_t n_processed = 0;
    uint32_t n_lines = 0;

    while (bcf_sr_next_line (sr)) {
        bcf1_t *line = bcf_sr_get_line(sr,0) ;
        //const char* chrom = bcf_seqname(hr,line) ;
        //if (!p->chrom) p->chrom = strdup (chrom) ;
        //else if (strcmp (chrom, p->chrom)) break ;
        //int pos = line->pos; // bcf coordinates are 0-based
        //char *ref, *REF;
        //ref = REF = strdup(line->d.allele[0]);
        //while ( (*ref = toupper(*ref)) ) ++ref ;
        if(line->n_allele != 2) {
            //std::cerr << "skip: n_allele=" << line->n_allele << std::endl;
            continue;
        }

        // get a copy of GTs
        mgt_arr = 0;
        int ngt = bcf_get_genotypes(hr, line, &gt_arr, &mgt_arr);
        if (ngt <= 0) continue;  // it seems that -1 is used if GT is not in the FORMAT
        if (ngt != 2*n_samples) {
            std::cerr << "not multiple of 2*n_s" << std::endl;
            continue;
        }
        //memcpy(gt_a, gt_arr, ngt*sizeof(int));
        int j = 0;
        for(int i = 0; i < 2*n_samples; i += 2, ++j) {
            gt1_a[j] = gt_arr[i+0];
            gt2_a[j] = gt_arr[i+1];
        }

        pbwt[0].Update(gt1_a);
        pbwt[1].Update(gt2_a);
        ++n_processed;
        ++n_lines;
        // encode
        n_in += 2*n_samples;

        if(pbwt[0].n_q2 > 100) {
            for(int i = 0; i < n_samples; ++i) {
                data_buffer[n_data++] = pbwt[0].prev[i];
            }
        } else {
            uint16_t s = 0;
            for(int i = 0; i < n_samples; ++i) {
                fmodel1[s].EncodeSymbol(&rc1, pbwt[0].prev[i]);
                s <<= 1;
                s |= (pbwt[0].prev[i] == 1);
                s &= 1023;
                _mm_prefetch((const char *)&fmodel1[s], _MM_HINT_T0);
                //data_buffer[n_data++] = pbwt[0].prev[i];
            }
        }

        if(pbwt[1].n_q2 > 100) {
            for(int i = 0; i < n_samples; ++i) {
                data_buffer[n_data++] = pbwt[1].prev[i];
            }
        } else {
            uint16_t s = 0;
            for(int j = 0; j < n_samples; ++j) {
                fmodel2[s].EncodeSymbol(&rc2, pbwt[1].prev[j]);
                s <<= 1;
                s |= (pbwt[1].prev[j] == 1);
                s &= 1023;
                _mm_prefetch((const char *)&fmodel2[s], _MM_HINT_T0);
                //data_buffer[n_data++] = pbwt[1].prev[j];
            }
        }

        if(n_lines == 75000) {
        //if(n_data > 9000000) {
            int ret = ZstdCompress(data_buffer,n_data,zstd_buffer,10000000,22);
            n_out += ret;
            n_out_zstd += ret;
            std::cerr << "Zstd " << n_lines << " Compressed=" << n_in << "->" << n_out << "(" << (float)n_in/n_out << "-fold)" << std::endl;
            //std::cerr << "Zstd Compressed=" << n_in << "->" << n_out_zstd << "(" << (float)n_in/n_out_zstd << "-fold)" << std::endl;
            n_data = 0;
        //}

        //if(rc1.OutSize() > 1000000 || rc2.OutSize() > 1000000) {
        //if(n_lines == 50000) {
            n_out += rc1.OutSize();
            rc1.FinishEncode();
            rc1.SetOutput(out1_buffer);
            rc1.StartEncode();
            //std::cerr << "RC1 " << n_lines << " Compressed=" << n_in << "->" << n_out << "(" << (float)n_in/n_out << "-fold)" << std::endl;
            delete[] fmodel1;
            fmodel1 = new pil::FrequencyModel<2>[1024];
            pbwt[0].reset();

            n_out += rc2.OutSize();
            rc2.FinishEncode();
            rc2.SetOutput(out2_buffer);
            rc2.StartEncode();
            std::cerr << "RC " << n_lines << " Compressed=" << n_in << "->" << n_out << "(" << (float)n_in/n_out << "-fold)" << std::endl;
            delete[] fmodel2;
            fmodel2 = new pil::FrequencyModel<2>[1024];
            pbwt[1].reset();
            n_lines = 0;
        }
    }

    if(n_data) {
        int ret = ZstdCompress(data_buffer,n_data,zstd_buffer,10000000,1);
        n_out_zstd += ret;
        n_out += ret;
        std::cerr << "Zstd " << n_lines << " Compressed=" << n_in << "->" << n_out << "(" << (float)n_in/n_out << "-fold)" << std::endl;
        //std::cerr << "Zstd Compressed=" << n_in << "->" << n_out_zstd << "(" << (float)n_in/n_out_zstd << "-fold)" << std::endl;
        n_data = 0;
    }

    if(rc1.OutSize() || rc2.OutSize()) {
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
    }

    delete[] fmodel1; delete[] fmodel2;
    if (gt_arr) free(gt_arr);
    delete[] data_buffer;
    delete[] out1_buffer; delete[] out2_buffer;
    delete[] zstd_buffer;
    delete[] gt1_a; delete[] gt2_a;
    delete[] pbwt;
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
