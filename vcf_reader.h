/*
* Copyright (c) 2019 Marcus D. R. Klarqvist
* Author(s): Marcus D. R. Klarqvist
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*   http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing,
* software distributed under the License is distributed on an
* "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
* KIND, either express or implied.  See the License for the
* specific language governing permissions and limitations
* under the License.
*/
#ifndef VCF_READER_H_
#define VCF_READER_H_

#include <iostream>
#include <string>
#include <memory>

#include <htslib/vcf.h>

namespace tachyon {
namespace io {

class VcfReader {
public:
    typedef VcfReader self_type;

public:
    // Singleton design pattern for retrieving a guaranteed unique pointer
    // to a VcfReader. This choice is to prevent inadvertent writes to the
    // target file as another file-handle is accessing it.
    static std::unique_ptr<self_type> FromFile(const std::string& variants_path, uint32_t n_extra_threads = 0){
        htsFile* fp = hts_open(variants_path.c_str(), "r");
        if(n_extra_threads){
            int ret = hts_set_threads(fp, n_extra_threads);
            if(ret < 0){
                std::cerr << "Failed to open multiple handles!" << std::endl;
                return nullptr;
            }
        }

        if (fp == nullptr) {
            std::cerr << "Could not open " << variants_path << std::endl;
            return nullptr;
        }

        bcf_hdr_t* header = bcf_hdr_read(fp);
        if (header == nullptr){
            std::cerr << "Couldn't parse header for " << fp->fn << std::endl;
            return nullptr;
        }

        return std::unique_ptr<self_type>(new self_type(variants_path, fp, header));
    }

    bool Next(const int unpack_level = BCF_UN_ALL){
        if (bcf_read(this->fp_, this->header_, this->bcf1_) < 0) {
            if (bcf1_->errcode) {
                std::cerr << "Failed to parse VCF record: " << bcf1_->errcode << std::endl;
                return false;
            } else {
                return false;
            }
        }

        bcf_unpack(this->bcf1_, unpack_level);
        return true;
    }

    bool Next(bcf1_t* bcf_entry, const int unpack_level = BCF_UN_ALL){
        if (bcf_read(this->fp_, this->header_, bcf_entry) < 0) {
            if (bcf_entry->errcode) {
                std::cerr << "Failed to parse VCF record: " << bcf1_->errcode << std::endl;
                return false;
            } else {
                //std::cerr << utility::timestamp("ERROR") << "Failed to retrieve a htslib bcf1_t record!" << std::endl;
                return false;
            }
        }

        bcf_unpack(bcf_entry, unpack_level);
        return true;
    }

private:
    // Private constructor.
    VcfReader(const std::string& variants_path,
              htsFile* fp,
              bcf_hdr_t* header) :
        n_samples_(0),
        fp_(fp),
        header_(header),
        bcf1_(bcf_init())
{
    if (this->header_->nhrec < 1) {
        std::cerr << "Empty header, not a valid VCF." << std::endl;
        return;
    }

    // Store the file-format header string
    if (std::string(this->header_->hrec[0]->key) != "fileformat") {
        std::cerr << "Not a valid VCF, fileformat needed: " << variants_path << std::endl;
    }

    // Fill in the contig info for each contig in the VCF header. Directly
    // accesses the low-level C struct because there are no indirection
    // macros/functions by htslib API.
    // BCF_DT_CTG: offset for contig (CTG) information in BCF dictionary (DT).
    const int n_contigs = this->header_->n[BCF_DT_CTG];
    for (int i = 0; i < n_contigs; ++i) {
        const bcf_idpair_t& idPair = this->header_->id[BCF_DT_CTG][i];
        //this->vcf_header_.AddContigInfo(idPair);
    }

    // Populate samples info.
    n_samples_ = bcf_hdr_nsamples(this->header_);
    //for (int i = 0; i < n_samples; i++)
    //    this->vcf_header_.AddSample(std::string(this->header_->samples[i]));

   // this->vcf_header_.BuildReverseMaps();

    // Build literal VCF header string for storage.
    kstring_t htxt = {0,0,0};
    bcf_hdr_format(this->header_, 0, &htxt);
    while (htxt.l && htxt.s[htxt.l-1] == '\0') --htxt.l; // kill trailing zeros
    std::string temp = std::string(htxt.s, htxt.l);
    size_t pos = temp.find("#CHROM"); // search for start of column header line
    temp = temp.substr(0, pos);
    //this->vcf_header_.literals_ = temp;
    free(htxt.s);
}

public:
    // Public destructor.
    ~VcfReader() {
        bcf_destroy(this->bcf1_);
        bcf_hdr_destroy(this->header_);
        hts_close(this->fp_);
    }

public:
    // Number of samples.
    int64_t n_samples_;

    // A pointer to the htslib file used to access the VCF data.
    htsFile * fp_;

    // A htslib header data structure obtained by parsing the header of this VCF.
    bcf_hdr_t * header_;

    // htslib representation of a parsed vcf line.
    bcf1_t* bcf1_;
};

}
}



#endif /* VCF_READER_H_ */
