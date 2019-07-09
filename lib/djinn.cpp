#include "djinn.h"

namespace djinn {

djinn_variant_t::djinn_variant_t() : 
    ploidy(0), n_allele(0), data(nullptr), data_len(0), 
    data_alloc(0), data_free(false), errcode(0), 
    unpacked(0), d(nullptr) 
{}

djinn_variant_t::~djinn_variant_t() {
    if (data_free) delete[] data;
    delete d;
}

int djinn_variant_t::ToVcfDebug(char* out, const char phasing) const {
    int ret = ToVcf(out, phasing);
    if (ret <= 1) return ret;
    if (out[ret-1] == '\n' && out[ret-2] != '\t') {
        // std::cerr << "not true=" << (int)out[ret-1] << "," << (int)out[ret-2] << std::endl;
        out[ret-1] = '\t';
        out[ret++] = '\n';
    }
    return ret;
}

int djinn_variant_t::ToVcf(char* out, const char phasing) const {
    if (out == nullptr) return -1;

    uint32_t len_vcf = 0;
    if (unpacked == DJN_UN_IND) {
        if (n_allele < 10) {
            for (int j = 0; j < data_len; j += ploidy) {
                out[len_vcf++] = (char)data[j+0] + '0';
                for (int i = 1; i < ploidy; ++i) {
                    out[len_vcf++] = phasing;
                    out[len_vcf++] = (char)data[j+1] + '0';
                }
                out[len_vcf++] = '\t';
            }
            out[len_vcf++] = '\n';
        } else {
            std::cerr << "not implemented" << std::endl;
        }
    } else if (unpacked == DJN_UN_EWAH) {
        assert(d != nullptr);
        if (n_allele < 10) {
            uint32_t n_out = 0;
            uint32_t diff  = 0;
            // mul: number of packed items in 32-bit dirty bitvectors.
            // mask: bitmap for dirty bitvectors (either 1-bit or 4-bit selector).
            // shift: bit-shift width in bits for unpacking dirty words.
            const uint32_t mul   = (d->dirty_type == DJN_DIRTY_2MC ? 32 :  8);
            const uint8_t  mask  = (d->dirty_type == DJN_DIRTY_2MC ?  1 : 15);
            const uint8_t  shift = (d->dirty_type == DJN_DIRTY_2MC ?  1 :  4);

            uint32_t to = 0;
            for (int i = 0; i < d->n_ewah; ++i) {
                // Emit clean words
                for (int j = 0; j < d->ewah[i]->clean; ++j) {
                    to = n_out + mul > d->n_samples ? d->n_samples - n_out : mul;
                    for (int k = 0; k < to; ++k, ++n_out) {
                        diff = n_out % ploidy;
                        if (diff == 0 && n_out != 0) out[len_vcf++] = '\t';
                        if (diff != 0) {
                            out[len_vcf++] = phasing;
                        }
                        out[len_vcf++] = (d->ewah[i]->ref & mask) + '0';
                    }
                }
                
                // Emit dirty words
                for (int j = 0; j < d->ewah[i]->dirty; ++j) {
                    uint32_t ref = *(d->dirty[i] + j); // copy
                    to = n_out + mul > d->n_samples ? d->n_samples - n_out : mul;
                    for (int k = 0; k < to; ++k, ++n_out) {
                        diff = n_out % ploidy;
                        if (diff == 0 && n_out != 0) out[len_vcf++] = '\t';
                        if (diff != 0) {
                            out[len_vcf++] = phasing;
                        }
                        out[len_vcf++] = (ref & mask) + '0';
                        ref >>= shift;    
                    }
                }
            }
            assert(n_out == d->n_samples);
            // out[len_vcf++] = '\t';
            out[len_vcf++] = '\n';
        } else {
            std::cerr << "not implemented" << std::endl;
            return -1;
        }

    } else {
        return -1;
    }
    
    return len_vcf;
}

}