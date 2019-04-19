#include <cstdlib>
#include <cmath>
#include <errno.h>

#include <memory>
#include <iostream>

#include <bitset>
#include <vector>

#include "range_coder.h"
#include "frequency_model.h"
#include "vcf_reader.h"

#include "model_ad.h"

#include <zstd.h>
#include <zstd_errors.h>

int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
}

int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int ret = ZSTD_decompress(out, out_capacity, in, n_in);
    return(ret);
}

const uint8_t TWK_BCF_GT_UNPACK[3] = {2, 0, 1};
const uint8_t TWK_BCF_GT_UNPACK_GENERAL[16] = {15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
#define BCF_UNPACK_GENOTYPE(A) TWK_BCF_GT_UNPACK[(A >> 1)]
#define BCF_UNPACK_GENOTYPE_GENERAL(A) TWK_BCF_GT_UNPACK_GENERAL[(A >> 1)]

class PBWT {
public:
    PBWT(int64_t n_samples, int n_symbols) :
        n_symbols(n_symbols),
        n_samples(n_samples),
        n_steps(0),
        prev(new uint8_t[n_samples]),
        ppa(new uint32_t[n_samples]),
        n_queue(new uint32_t[n_symbols]),
        queue(new uint32_t*[n_symbols])
    {
        assert(n_symbols > 1);

        for (int i = 0; i < n_symbols; ++i)
            queue[i] = new uint32_t[n_samples];

        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);
        reset();
    }

    ~PBWT() {
        delete[] prev;
        delete[] ppa;

        for (int i = 0; i < n_symbols; ++i)
            delete[] queue[i];

        delete[] queue;
        delete[] n_queue;
    }

    void Initiate(int64_t n_samples) {
        exit(1);
        delete[] prev;
        delete[] ppa;
        for (int i = 0; i < n_symbols; ++i)
            delete[] queue[i];

        this->n_samples = n_samples;
        prev = new uint8_t[n_samples];
        ppa  = new uint32_t[n_samples];
        for (int i = 0; i < n_symbols; ++i)
            queue[i] = new uint32_t[n_samples];

        reset();
    }

    void reset() {
        if (ppa != nullptr) {
            for (int i = 0; i < n_samples; ++i)
                ppa[i] = i;
        }
        if (prev != nullptr) {
            memset(prev, 0, sizeof(uint8_t)*n_samples);
        }
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);
    }

    int Update(const int* arr) {
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

        for (int i = 0; i < n_samples; ++i) {
            const uint8_t& gt = BCF_UNPACK_GENOTYPE(arr[ppa[i]]);
            for (int j = 0; j < n_symbols; ++j) {
                if (gt == j)
                    queue[j][n_queue[j]++] = ppa[i];
            }
            prev[i] = gt;
            //std::cerr << " " << (int)gt;
        }
        //std::cerr << std::endl;

        uint32_t of = 0;
        for (int j = 0; j < n_symbols; ++j) {
            for (int i = 0; i < n_queue[j]; ++i, ++of) {
                ppa[of] = queue[j][i];
            }
        }
        //for (int i = 0; i < n_q2; ++i, ++of) ppa[of] = q2[i];
        assert(of == n_samples);
        // debug should be sorted here
        //for (int i = 0; i < n_samples; ++i) {
         //   std::cerr << " " << (int)BCF_UNPACK_GENOTYPE(arr[ppa[i]]);
        //}
        //std::cerr << std::endl;
        ++n_steps;

        return(1);
    }

    int Update(const uint8_t* arr, uint32_t stride = 1) {
        //n_q1 = n_q2 = 0;
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

        for (int i = 0; i < n_samples; ++i) {
            const uint8_t& gt = BCF_UNPACK_GENOTYPE(arr[ppa[i] * stride]);
            assert(gt < n_symbols);
            for (int j = 0; j < n_symbols; ++j) {
                if (gt == j)
                    queue[j][n_queue[j]++] = ppa[i];
            }
            prev[i] = gt;
            //std::cerr << " " << (int)gt;
        }
        //std::cerr << std::endl;

        uint32_t of = 0;
        for (int j = 0; j < n_symbols; ++j) {
            for (uint32_t i = 0; i < n_queue[j]; ++i, ++of)
                ppa[of] = queue[j][i];
        }
        //for (int i = 0; i < n_q2; ++i, ++of) ppa[of] = q2[i];
        //std::cerr << "of=" << of << "/" << n_samples << std::endl;
        assert(of == n_samples);
        // debug should be sorted here
        //std::cerr << ToPrettyString() << std::endl;
        //std::cerr << "Alts=" << n_queue[1] << std::endl;
        //for (int i = 0; i < n_samples; ++i) {
        //   std::cerr << " " << (int)BCF_UNPACK_GENOTYPE(arr[ppa[i]*stride]);
        //}
        //std::cerr << std::endl;
        ++n_steps;

        return(1);
    }

    int UpdateGeneral(const uint8_t* arr, uint32_t stride = 1) {
        // Reset queues.
        memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

        for (int i = 0; i < n_samples; ++i) {
            const uint8_t& gt = (arr[ppa[i] * stride]);
            assert(gt < n_symbols);
            for (int j = 0; j < n_symbols; ++j) {
                if (gt == j)
                    queue[j][n_queue[j]++] = ppa[i];
            }
            prev[i] = gt;
        }

        uint32_t of = 0;
        for (int j = 0; j < n_symbols; ++j) {
            for (uint32_t i = 0; i < n_queue[j]; ++i, ++of)
                ppa[of] = queue[j][i];
        }
        assert(of == n_samples);
        ++n_steps;

        return(1);
    }

    std::string ToPrettyString() const {
        std::string ret = "n=" + std::to_string(n_samples) + " {";
        ret += std::to_string(ppa[0]);
        for (int i = 1; i < n_samples; ++i) {
            ret += "," + std::to_string(ppa[i]);
        }
        ret += "}";
        return(ret);
    }

public:
    const int  n_symbols; // universe of symbols (number of unique symbols)
    int64_t    n_samples; // number of samples (free interpretation)
    uint64_t   n_steps; // number of updates made (debugging)
    uint8_t*   prev; // previous output array
    uint32_t*  ppa; // current PPA
    uint32_t*  n_queue; // number of elements in each positional queue
    uint32_t** queue; // the positional queues themselves
};

#define MODEL_SIZE 65536

class GeneralPBWTModel {
public:
    GeneralPBWTModel() noexcept :
        max_model_symbols(0),
        model_context_shift(0),
        model_context(0),
        n_buffer(10000000), buffer(new uint8_t[n_buffer])
    {}

    GeneralPBWTModel(int64_t n_samples, int n_symbols) :
        max_model_symbols(n_symbols),
        model_context_shift(ceil(log2(n_symbols))),
        model_context(0),
        pbwt(std::make_shared<PBWT>(n_samples, n_symbols)),
        range_coder(std::make_shared<pil::RangeCoder>()),
        n_buffer(10000000),
        buffer(new uint8_t[n_buffer])
    {
        assert(n_symbols > 1);
        models.resize(MODEL_SIZE);
        for (int i = 0; i < MODEL_SIZE; ++i) models[i] = std::make_shared<pil::FrequencyModel>();

        Reset();
    }

    // Constructor outside constructor.
    void Construct(int64_t n_samples, int n_symbols) {
        max_model_symbols = n_symbols;
        model_context_shift = ceil(log2(n_symbols));
        model_context = 0;
        pbwt = std::make_shared<PBWT>(n_samples, n_symbols);
        range_coder = std::make_shared<pil::RangeCoder>();

        assert(n_symbols > 1);
        models.resize(MODEL_SIZE);
        for (int i = 0; i < MODEL_SIZE; ++i) models[i] = std::make_shared<pil::FrequencyModel>();

        Reset();
    }

    void ResetModels() {
        for (int i = 0; i < MODEL_SIZE; ++i) {
            models[i]->Initiate(max_model_symbols, max_model_symbols);
            //models[i]->SetShift(11);
        }
    }

    void ResetPBWT() {
        if (pbwt.get() != nullptr)
            pbwt->reset();
    }

    void ResetContext() { model_context = 0; }

    void Reset() {
        ResetModels();
        ResetPBWT();
        ResetContext();
    }

    int FinishEncoding() {
        range_coder->FinishEncode();
        int ret = range_coder->OutSize();
        return(ret);
    }

    void StartEncoding() {
        range_coder->SetOutput(buffer);
        range_coder->StartEncode();
    }

    void EncodeSymbol(const uint16_t symbol) {
        models[model_context]->EncodeSymbol(range_coder.get(), symbol);
        model_context <<= model_context_shift;
        model_context |= symbol;
        model_context &= (MODEL_SIZE-1);
        _mm_prefetch((const char *)&(models[model_context]), _MM_HINT_T0);
    }

    ~GeneralPBWTModel() {
        delete[] buffer;
    }

public:
    int max_model_symbols;
    int model_context_shift;
    uint32_t model_context;
    std::shared_ptr<PBWT> pbwt;
    std::shared_ptr<pil::RangeCoder> range_coder;
    std::vector < std::shared_ptr<pil::FrequencyModel> > models;
    size_t n_buffer;
    uint8_t* buffer; // fix
};

struct Buffer {
    Buffer() noexcept : len(0), cap(0), data(nullptr){}
    Buffer(const size_t size) noexcept : len(0), cap(size), data(new uint8_t[size]){}
    ~Buffer() { delete[] data; }

    const size_t& size() const { return len; }
    const size_t& capacity() const { return cap; }

    int resize() {
        uint8_t* old = data;
        size_t new_cap = len * 1.2 - len < 65536 ? 65536 : len * 1.2;
        data = new uint8_t[new_cap];
        cap = new_cap;
        memcpy(data,old,len);
        delete[] old;
        
        return 1;
    }

    int resize(const size_t desired_size) {
        if (desired_size < cap) {
            len = desired_size < len ? desired_size : len;
            return 1;
        }

        uint8_t* old = data;
        data = new uint8_t[desired_size];
        cap = desired_size;
        memcpy(data,old,len);
        delete[] old;
        
        return(0);
    }

    void reset() { len = 0; }

    size_t len, cap;
    uint8_t* data;
};

class GenotypePermuter {
public:
    GenotypePermuter(int64_t n_s) :
        n_samples(n_s),
        block_size(8192),
        processed_lines(0),
        processed_lines_local(0),
        bytes_in(0), bytes_out(0),
        models(nullptr)
    {
        base_models[0].Construct(n_samples, 2);
        base_models[1].Construct(n_samples, 2);
        base_models[2].Construct(n_samples, 3);
        base_models[3].Construct(n_samples, 3);
        base_models_complex[0].Construct(n_samples, 16);
        base_models_complex[1].Construct(n_samples, 16);

        buf_general[0].resize(10000000);
        buf_general[1].resize(10000000);
        buf_raw.resize(10000000);

        base_models[0].StartEncoding();
        base_models[1].StartEncoding();
        base_models[2].StartEncoding();
        base_models[3].StartEncoding();
        base_models_complex[0].StartEncoding();
        base_models_complex[1].StartEncoding();
    }

    ~GenotypePermuter() {
        delete[] models;
    }

    inline bool CheckLimit() const {
        return (processed_lines_local == 8196 || 
                base_models[0].range_coder->OutSize() > 9000000 || 
                base_models[1].range_coder->OutSize() > 9000000 || 
                base_models[2].range_coder->OutSize() > 9000000 || 
                base_models[3].range_coder->OutSize() > 9000000 ||
                base_models_complex[0].range_coder->OutSize() > 9000000 ||
                base_models_complex[1].range_coder->OutSize() > 9000000 ||
                buf_raw.len > 9000000);
    }

    //
    int Encode(const bcf1_t* bcf) {
        if (bcf == NULL) return 0;

        bytes_in += bcf->d.fmt[0].p_len;
        
        return(Encode2N(bcf));
    }

private:
    // Diploid wrapper.
    int Encode2N(const bcf1_t* bcf) {
        if (bcf->n_allele == 2) return(Encode2N2M(bcf->d.fmt));
        else if (bcf->n_allele < 4) return(Encode2NXM(bcf->d.fmt)); 
        else {
            // std::cerr << "alleles=" << bcf->n_allele << std::endl;
            uint8_t* gts = bcf->d.fmt[0].p;
            for (int i = 0; i < 2*n_samples; ++i) {
                buf_raw.data[buf_raw.len++] = gts[i];
            }

            ++processed_lines_local;
            ++processed_lines;
        }
    }

    // Wrapper for 2N2M
    int Encode2N2M(const bcf_fmt_t* fmt) {
        if (CheckLimit()) {
            Compress();
        }

        // Todo: assert genotypes are set for this variant.
        const uint8_t* gts = fmt[0].p; // data pointer
        
        // Todo: Assert that total number of alleles < 15.
        uint32_t alts[256] = {0};
        for (int i = 0; i < 2*n_samples; ++i) {
            ++alts[BCF_UNPACK_GENOTYPE_GENERAL(gts[i])];
        }

        if (alts[15] == 0) { // No missing values.
            return Encode2N2MC(fmt);

        } else { // Having missing values.
            // std::cerr << "using extended model" << std::endl;
            return Encode2N2MM(fmt);
        }

        return 1;
    }

    // 2N2M complete
    int Encode2N2MC(const bcf_fmt_t* fmt) {
        const uint8_t* gts = fmt[0].p; // data pointer
        
        base_models[0].pbwt->Update(&gts[0], 2);
        base_models[1].pbwt->Update(&gts[1], 2);

        base_models[0].ResetContext();
        for (int j = 0; j < n_samples; ++j) {
            assert(base_models[0].pbwt->prev[j] < 2);
            base_models[0].EncodeSymbol(base_models[0].pbwt->prev[j]);
        }
    
        base_models[1].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models[1].pbwt->prev[i] < 2);
            base_models[1].EncodeSymbol(base_models[1].pbwt->prev[i]);   
        }

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }

    // 2N2M with missing
    int Encode2N2MM(const bcf_fmt_t* fmt) {
        const uint8_t* gts = fmt[0].p; // data pointer
        
        base_models[3].pbwt->Update(&gts[0], 2);
        base_models[4].pbwt->Update(&gts[1], 2);

        base_models[3].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models[3].pbwt->prev[i] < 3);
            base_models[3].EncodeSymbol(base_models[3].pbwt->prev[i]);
        }

        base_models[4].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models[4].pbwt->prev[i] < 3);
            base_models[4].EncodeSymbol(base_models[4].pbwt->prev[i]);
        }

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }

    // 2N any M (up to 16)
    int Encode2NXM(const bcf_fmt_t* fmt) {
        if (CheckLimit()) {
            Compress();
        }

        // Todo: assert genotypes are set for this variant.
        const uint8_t* gts = fmt[0].p; // data pointer
        base_models_complex[0].pbwt->UpdateGeneral(&gts[0], 2);
        
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models_complex[0].pbwt->prev[i] < 16);
            base_models_complex[0].EncodeSymbol(base_models_complex[0].pbwt->prev[i]);
        }

        base_models_complex[1].pbwt->UpdateGeneral(&gts[1], 2);

        for (int i = 0; i < n_samples; ++i) {
            assert(base_models_complex[1].pbwt->prev[i] < 16);
            base_models_complex[1].EncodeSymbol(base_models_complex[1].pbwt->prev[i]);
        }

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }

    int Compress() {
        // flush: temp
        int p1  = base_models[0].FinishEncoding();
        int p2  = base_models[1].FinishEncoding();
        int p1E = base_models[2].FinishEncoding();
        int p2E = base_models[3].FinishEncoding();
        int p2X = base_models_complex[0].FinishEncoding();
        int p2X2 = base_models_complex[1].FinishEncoding();
        int praw = ZstdCompress(buf_raw.data, buf_raw.len,
                                buf_general[0].data, buf_general[0].capacity(),
                                10);

        base_models[0].Reset();
        base_models[0].StartEncoding();
        base_models[1].Reset();
        base_models[1].StartEncoding();
        base_models[2].Reset();
        base_models[2].StartEncoding();
        base_models[3].Reset();
        base_models[3].StartEncoding();
        base_models_complex[0].Reset();
        base_models_complex[0].StartEncoding();
        base_models_complex[1].Reset();
        base_models_complex[1].StartEncoding();

        processed_lines_local = 0;
        std::cerr << "Flushed: " << p1 << "," << p2 << "," << p1E << "," << p2E << "," << p2X << "," << p2X2 << "," << praw << std::endl;
        buf_raw.reset();
        bytes_out += p1+p2+p1E+p2E+p2X+p2X2+praw;
        std::cerr << "[PROGRESS] " << bytes_in << "->" << bytes_out << " (" << (double)bytes_in/bytes_out << "-fold)" << std::endl;
    }

public:
    int64_t n_samples;
    uint32_t block_size;
    uint32_t processed_lines;
    uint32_t processed_lines_local;
    uint64_t bytes_in, bytes_out;
    GeneralPBWTModel base_models[4];
    GeneralPBWTModel base_models_complex[2];
    Buffer buf_general[2];
    Buffer buf_raw;
    GeneralPBWTModel* models;
};

void ReadVcfGT (const std::string& filename) {
    std::unique_ptr<tachyon::io::VcfReader> reader = tachyon::io::VcfReader::FromFile(filename,8);
    if (reader.get() == nullptr) {
        std::cerr << "failed read" << std::endl;
        exit(1);
    }

    uint64_t n_out = 0, n_in = 0, n_data = 0, n_out1 = 0, n_out2 = 0;
    uint64_t n_out_gtpbwt = 0;
    //uint8_t* data_buffer = new uint8_t[2*10000000];
    uint8_t* gtpbwt_buffer  = new uint8_t[reader->n_samples_];

    uint8_t* out1_buffer = new uint8_t[10000000];
    uint8_t* out2_buffer = new uint8_t[10000000];
    uint8_t* out3_buffer = new uint8_t[10000000];
    uint8_t* out_gtpbwt  = new uint8_t[10000000];

    pil::RangeCoder rc_pbwt, rc_pbwt4;
    pil::FrequencyModel* gtpbwt_model = new pil::FrequencyModel[MODEL_SIZE];
    for (int i = 0; i < MODEL_SIZE; ++i) gtpbwt_model[i].Initiate(4,4);

    //pil::FrequencyModel<16>* gtpbwt_model_4 = new pil::FrequencyModel<16>[MODEL_SIZE]; // order-7 model
    pil::FrequencyModel* gtpbwt_model_4 = new pil::FrequencyModel[MODEL_SIZE];
    for (int i = 0; i < MODEL_SIZE; ++i) gtpbwt_model_4[i].Initiate(16,16);

    rc_pbwt.SetOutput(out_gtpbwt);
    rc_pbwt4.SetOutput(out3_buffer);
    rc_pbwt.StartEncode();
    rc_pbwt4.StartEncode();

    //PBWT pbwt[2] = {{reader->n_samples_, 2}, {reader->n_samples_, 2}};
    PBWT gt_pbwt(reader->n_samples_,   4);
    PBWT gt_pbwt4(reader->n_samples_, 16);

    PBWT pbwt_4[2] = {{reader->n_samples_, 4}, {reader->n_samples_, 4}};

    uint32_t n_processed = 0;
    uint32_t n_lines = 0, n_lines_total = 0;
    //uint32_t last_flush_pos = 0;

    uint32_t n_bstep = reader->n_samples_/32;
    size_t n_bitmaps = std::ceil((float)reader->n_samples_/n_bstep);
    std::cerr << "n_samples=" << reader->n_samples_ << std::endl;
    std::cerr << "n_bitmaps=" << n_bitmaps << std::endl;
    uint32_t* bitmaps1 = new uint32_t[n_bitmaps];
    uint32_t* bitmaps2 = new uint32_t[n_bitmaps];

    // test pbwt model
    GeneralPBWTModel pmodel1(reader->n_samples_, 2);
    GeneralPBWTModel pmodel2(reader->n_samples_, 2);
    pmodel1.StartEncoding();
    pmodel2.StartEncoding();

    GeneralPBWTModel pmodel1E(reader->n_samples_, 3);
    GeneralPBWTModel pmodel2E(reader->n_samples_, 3);
    pmodel1E.StartEncoding();
    pmodel2E.StartEncoding();

    // AD: temp
    // FormatAlelleDepth fmt_ad(10000,reader->n_samples_,true);

    GenotypePermuter gtperm(reader->n_samples_);
    
    // While there are bcf records available.
    while (reader->Next()) {
        //const char* chrom = bcf_seqname(hr,line) ;
        //if (!p->chrom) p->chrom = strdup (chrom) ;
        //else if (strcmp (chrom, p->chrom)) break ;
        //int pos = line->pos; // bcf coordinates are 0-based
        //char *ref, *REF;
        //ref = REF = strdup(line->d.allele[0]);
        //while ( (*ref = toupper(*ref)) ) ++ref ;

        // const bcf_fmt_t* fmt = bcf_get_fmt(reader->header_, reader->bcf1_, "AD");
        // if (reader->bcf1_->n_allele == 2) fmt_ad.AddRecord(fmt);
         gtperm.Encode(reader->bcf1_);

        // 
        if (reader->bcf1_->n_allele > 4) {
            uint8_t* gts = reader->bcf1_->d.fmt[0].p;
            for (int i = 0; i < 2*reader->n_samples_; ++i) {
                out1_buffer[n_out1++] = gts[i];
            }

            continue;
        }

        ++n_processed;
        ++n_lines;
        ++n_lines_total;

        n_in += reader->bcf1_->d.fmt[0].p_len;

       
        
        if (reader->bcf1_->n_allele != 2) {
            // Todo: fix
            // continue;
            //std::cerr << "not2" << std::endl;
            // copy to 4-bit
            uint8_t* gts = reader->bcf1_->d.fmt[0].p;
            //uint32_t j = 0;
            //for (int i = 0; i < 2*reader->n_samples_; ++j) {
                //gtpbwt_buffer[j] = (BCF_UNPACK_GENOTYPE_GENERAL(gts[i+0]) << 2) |  BCF_UNPACK_GENOTYPE_GENERAL(gts[i+1]);
            //}
            pbwt_4[0].Update(&gts[0], 2);
            pbwt_4[1].Update(&gts[1], 2);

            //size_t before = rc_pbwt.OutSize();
            {
               uint32_t s = 0;
               for (int i = 0; i < reader->n_samples_; ++i) {
                   gtpbwt_model_4[s].EncodeSymbol(&rc_pbwt4, pbwt_4[0].prev[i]);
                   s <<= 4;
                   s |= pbwt_4[0].prev[i];
                   s &= (MODEL_SIZE-1);
                   _mm_prefetch((const char *)&gtpbwt_model_4[s], _MM_HINT_T0);
               }
            }

            {
               uint32_t s = 0;
               for (int i = 0; i < reader->n_samples_; ++i) {
                   gtpbwt_model_4[s].EncodeSymbol(&rc_pbwt4, pbwt_4[1].prev[i]);
                   s <<= 4;
                   s |= pbwt_4[1].prev[i];
                   s &= (MODEL_SIZE-1);
                   _mm_prefetch((const char *)&gtpbwt_model_4[s], _MM_HINT_T0);
               }
            }
            //std::cerr << rc_pbwt.OutSize() - before << " b for " << reader->bcf1_->n_allele << std::endl;*/
        }

        else
        {
            
            
            uint8_t* gts = reader->bcf1_->d.fmt[0].p;
            //pbwt[0].Update(&gts[0], 2);
            //pbwt[1].Update(&gts[1], 2);

            uint32_t alts[256] = {0};
            for (int i = 0; i < 2*reader->n_samples_; ++i) {
                ++alts[BCF_UNPACK_GENOTYPE_GENERAL(gts[i])];
            }

            /*for (int i = 0; i < 256; ++i) {
                if (alts[i] != 0) {
                    std::cerr << " " << i << ":" << alts[i];
                }
            }
            std::cerr << std::endl;*/

            if (alts[15] == 0) {

                pmodel1.pbwt->Update(&gts[0], 2);
                pmodel2.pbwt->Update(&gts[1], 2);

            {
               uint32_t n_bused = 0;
               memset(bitmaps1, 0, n_bitmaps*sizeof(uint32_t));
               for (int i = 0; i < reader->n_samples_; ++i) {
                   //bitmaps1[i/n_bstep] += (pbwt[0].prev[i] == 1);
                   bitmaps1[i/n_bstep] += (pmodel1.pbwt->prev[i] == 1);
               }

               pmodel1.ResetContext();

               //size_t before = rc1.OutSize();
                  uint32_t n_seen = 0;
                  uint32_t p = 0, n_p = n_bstep;
                  for (int i = 0; i < n_bitmaps; ++i) {

                          ++n_bused;
                          for (int j = 0; j < n_p; ++j) {
                              n_seen += (pmodel1.pbwt->prev[p+j] == 1);
                              pmodel1.EncodeSymbol(pmodel1.pbwt->prev[p+j]);
                          }
                        
                      p += n_bstep;
                      n_p = (p + n_bstep > reader->n_samples_ ? reader->n_samples_ - p : n_bstep);
                  }
                  assert(p + n_p == reader->n_samples_);
                  //if ((rc1.OutSize() - before) > pbwt[0].n_queue[1]*sizeof(uint16_t))
                  //std::cerr << pbwt[1].n_queue[1] << "\t" <<  n_bused << "\t" << rc1.OutSize() - before << "\t" << (rc1.OutSize() - before) / ((float)pbwt[0].n_queue[1]*sizeof(uint16_t)) << std::endl;
                  //std::cerr << "Seen=" << n_seen << "/" << pbwt[1].n_queue[1] << std::endl;
                  //assert(n_seen == pbwt[0].n_queue[1]);
                  assert(n_seen == pmodel1.pbwt->n_queue[1]);
            }


            {
               
               uint32_t n_seen = 0;
               uint32_t p = 0, n_p = n_bstep, s = 0;

               pmodel2.ResetContext();

               for (int i = 0; i < n_bitmaps; ++i) {
                       //std::cerr << "i=" << i << "/" << n_bitmaps << " p/np=" << p << ":" << n_p << "/" << reader->n_samples_  << std::endl;
                       for (int j = 0; j < n_p; ++j) {
                           //std::cerr << p+j << "->" << p+n_p << "/" << reader->n_samples_ << ": " << (int)pbwt[1].prev[p+j] << std::endl;
                           //assert(p+j < reader->n_samples_);
                           //assert(pmodel2.model_context == s); // assert parity in context
                           pmodel2.EncodeSymbol(pmodel2.pbwt->prev[p+j]);
                           n_seen += (pmodel2.pbwt->prev[p+j] == 1);

                           /*fmodel2[s].EncodeSymbol(&rc2, pbwt[1].prev[p+j]);
                           n_seen += (pbwt[1].prev[p+j] == 1);
                           s <<= 1;
                           s |= pbwt[1].prev[p+j];
                           s &= (MODEL_SIZE - 1);
                           _mm_prefetch((const char *)&fmodel2[s], _MM_HINT_T0);*/
                       }
                 
                   p += n_bstep;
                   n_p = (p + n_bstep > reader->n_samples_ ? reader->n_samples_ - p : n_bstep);
               }
               //std::cerr << pbwt[1].n_queue[1] << "\t" <<  n_bused << "\t" << rc2.OutSize() - before << "\t" << (rc2.OutSize() - before) / ((float)pbwt[0].n_queue[1]) << std::endl;
               //std::cerr << "Seen=" << n_seen << "/" << pbwt[1].n_queue[1] << std::endl;
               //assert(n_seen == pbwt[1].n_queue[1]);
               assert(n_seen == pmodel2.pbwt->n_queue[1]);



              /* uint32_t s = 0;
               for (int i = 0; i < reader->n_samples_; ++i) {
                   fmodel2[s].EncodeSymbol(&rc2, pbwt[1].prev[i]);
                   s <<= 1;
                   s |= pbwt[1].prev[i];
                   s &= (MODEL_SIZE-1);
                   _mm_prefetch((const char *)&fmodel2[s], _MM_HINT_T0);
               }*/
            }

            } else {
                // std::cerr << "using extended model with " << reader->bcf1_->n_allele << " alleles" << std::endl;
                pmodel1E.pbwt->Update(&gts[0], 2);
                pmodel2E.pbwt->Update(&gts[1], 2);
                //std::cerr << "done" << std::endl;

                pmodel1E.ResetContext();
                for (int i = 0; i < reader->n_samples_; ++i) {
                    pmodel1E.EncodeSymbol(pmodel1E.pbwt->prev[i]);
                }

                pmodel2E.ResetContext();
                for (int i = 0; i < reader->n_samples_; ++i) {
                    pmodel2E.EncodeSymbol(pmodel2E.pbwt->prev[i]);
                }
            }
        }

        if (n_lines == 8196 || pmodel1.range_coder->OutSize() > 9000000 || pmodel2.range_coder->OutSize() > 9000000 || pmodel1E.range_coder->OutSize() > 9000000 || pmodel2E.range_coder->OutSize() > 9000000 || rc_pbwt4.OutSize() > 9000000 || n_out1 > 9000000) {
            int zstd_ret = 0;
            if (n_out1) {
                int ret = ZstdCompress(out1_buffer,n_out1,out2_buffer,10000000,6);
                n_out += ret;
                std::cerr << "Zstd " << n_lines << " Compressed=" << n_out1 << "->" << ret << "(" << (float)n_out1/ret << "-fold)" << std::endl;
                n_out1 = 0;
                n_out_gtpbwt += ret;
                zstd_ret = ret;
            }

            //
            //std::cerr << "done here" << std::endl;
            //std::cerr << "test rc=" << pmodel1.FinishEncoding() << "," << pmodel2.FinishEncoding() << std::endl;
            int p1 = pmodel1.FinishEncoding();
            int p2 = pmodel2.FinishEncoding();
            int p1E = pmodel1E.FinishEncoding();
            int p2E = pmodel2E.FinishEncoding();
            n_out_gtpbwt += p1;
            n_out_gtpbwt += p2;
            n_out_gtpbwt += p1E;
            n_out_gtpbwt += p2E;
            pmodel1.Reset();
            pmodel1.StartEncoding();
            pmodel2.Reset();
            pmodel2.StartEncoding();
            pmodel1E.Reset();
            pmodel1E.StartEncoding();
            pmodel2E.Reset();
            pmodel2E.StartEncoding();

            //

            //n_out_gtpbwt += rc_pbwt.OutSize();
            //rc1.FinishEncode();
            //rc2.FinishEncode();
            rc_pbwt4.FinishEncode();

            //std::cerr << "RC: " << rc1.OutSize() << "," << rc2.OutSize() << "," << rc_pbwt4.OutSize() << std::endl;
            n_out_gtpbwt += rc_pbwt4.OutSize();

            std::cerr << "Results " << p1 << " " << p2 << " " << p1E << " " << p2E << " " << rc_pbwt4.OutSize() << " " << zstd_ret << std::endl;
            std::cerr << "gtPBWT " << n_lines << " Compressed=" << n_in << "->" << n_out_gtpbwt << "(" << (float)n_in/n_out_gtpbwt << "-fold)" << std::endl;


            //n_out_gtpbwt += rc1.OutSize();
            //n_out_gtpbwt += rc2.OutSize();


            //rc1.SetOutput(out1_buffer);
            //rc2.SetOutput(out2_buffer);

            rc_pbwt4.SetOutput(out3_buffer);
            rc_pbwt4.StartEncode();


            for (int i = 0; i < MODEL_SIZE; ++i) {
                //fmodel1[i].Initiate(2,2);
                //fmodel2[i].Initiate(2,2);
                gtpbwt_model_4[i].Initiate(16,16);
                //fmodel2[i].SetShift(11);
                //fmodel2[i].SetShift(11);
            }

            n_lines = 0;

            //pbwt[0].reset();
            //pbwt[1].reset();
            pbwt_4[0].reset();
            pbwt_4[1].reset();
            gt_pbwt.reset();
            gt_pbwt4.reset();
        }
    }

    if (pmodel1.range_coder->OutSize() || pmodel2.range_coder->OutSize() || rc_pbwt4.OutSize() || n_out1) {
        if (n_out1) {
            int ret = ZstdCompress(out1_buffer,n_out1,out2_buffer,10000000,6);
            n_out += ret;
            std::cerr << "Zstd " << n_lines << " Compressed=" << n_out1 << "->" << ret << "(" << (float)n_out1/ret << "-fold)" << std::endl;
            n_out1 = 0;
        }

        n_out_gtpbwt += rc_pbwt4.OutSize();
        //n_out_gtpbwt += rc1.OutSize();
        //n_out_gtpbwt += rc2.OutSize();
        n_out_gtpbwt += pmodel1.FinishEncoding();
        n_out_gtpbwt += pmodel2.FinishEncoding();

        std::cerr << "gtPBWT " << n_lines << " Compressed=" << n_in << "->" << n_out_gtpbwt << "(" << (float)n_in/n_out_gtpbwt << "-fold)" << std::endl;
    }

    std::cerr << "gtPBWT " << n_lines_total << " Compressed=" << n_in << "->" << n_out_gtpbwt << "(" << (float)n_in/n_out_gtpbwt << "-fold)" << std::endl;

    std::cerr << "here" << std::endl;
    gtperm.Compress(); // final
    std::cerr << "Final=" << gtperm.bytes_in << "->" << gtperm.bytes_out << " (" << (double)gtperm.bytes_in/gtperm.bytes_out << ")" << std::endl;

    // std::cerr << "AD cost=" << fmt_ad.tot_in << "->" << fmt_ad.tot_out << std::endl;

    //delete[] fmodel1; delete[] fmodel2;
    delete[] out1_buffer; delete[] out2_buffer; delete[] out3_buffer;
    delete[] gtpbwt_buffer;
    delete[] gtpbwt_model;

    // delete[] out_bitmap;
    // delete[] bitmap_model;
}

int main(int argc, char** argv) {
    std::cerr << argc << std::endl;
    if (argc == 1) return(1);
    else {
        std::cerr << std::string(argv[1]) << std::endl;
        ReadVcfGT(std::string(argv[1]));
        return(0);
    }
}
