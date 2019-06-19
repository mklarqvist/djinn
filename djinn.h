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
#ifndef DJINN_H_
#define DJINN_H_

#include <cstddef>//size_t
#include <cstdint>//uint

#include <cstdlib>//malloc
#include <cstring>//memcpy
#include <cassert>//assert

#include <iostream>//debug

#include <bitset>//debug

namespace djinn {

// Map missing to 2, 1->0, 2->1, and EOV -> 3.
const uint8_t TWK_BCF_GT_UNPACK[65] = {2,0,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3};
const uint8_t TWK_BCF_GT_PACK[3]   = {1, 2, 0};

// MISSING -> 14, EOV -> 15, other values as normal
const uint8_t TWK_BCF_GT_UNPACK_GENERAL[131] = 
{14,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,
47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,
62,63,15,65,66,67,68,69,70,71,72,73,74,75,76,
77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,
92,93,94,95,96,97,98,99,100,101,102,103,104,
105,106,107,108,109,110,111,112,113,114,115,
116,117,118,119,120,121,122,123,124,125,126,
127,128,129};

const uint8_t TWK_BCF_GT_UNPACK_GENERAL_REV[131] = 
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,15,16,17,18,19,20,21,22,
23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,
43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,
63,15,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,
83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,
103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,
118,119,120,121,122,123,124,125,126,127,128,129};


#define BCF_UNPACK_GENOTYPE(A) TWK_BCF_GT_UNPACK[(A) >> 1]
#define BCF_UNPACK_GENOTYPE_GENERAL(A) TWK_BCF_GT_UNPACK_GENERAL[(A) >> 1]

// EWAH structure
#pragma pack(push, 1)
struct djinn_ewah_t {
    djinn_ewah_t() : ref(0), clean(0), dirty(0){}
    void reset() { ref = clean = dirty = 0; }

    uint64_t ref: 4, clean: 30, dirty: 30;
};
#pragma pack(pop)

// Header structure
struct djinn_hdr_t {
    uint8_t version[3];
    uint8_t base_ploidy;
    int64_t n_samples; // number of samples
};

// Controller contexts
struct djinn_ctx_ctrl_t {
    uint16_t pbwt: 1, 
             dimc: 1, 
             dim: 1, 
             nm: 1, 
             unused: 12;
};

struct djinn_wah_ctrl_t {
    uint16_t pbwt:   1, 
             di2mc1: 1, 
             di2mc2: 1, 
             di2mm1: 1, 
             di2mm2: 1, 
             di2x1:  1, 
             di2x2:  1, 
             unused: 9;
};

#define DJINN_CTX_MODEL_ARC  0
#define DJINN_CTX_MODEL_DIMC 1
#define DJINN_CTX_MODEL_DIM  2
#define DJINN_CTX_MODEL_NM   3

#define DJINN_WAH_MODEL_DI2MC1      0
#define DJINN_WAH_MODEL_DI2MC2      1
#define DJINN_WAH_MODEL_DI2MM1      2
#define DJINN_WAH_MODEL_DI2MM2      3
#define DJINN_WAH_MODEL_DI2X1       4
#define DJINN_WAH_MODEL_DI2X2       5

struct djinn_data_container_t {
    djinn_data_container_t() : n(0), n_c(0), n_v(0), vptr(nullptr), vptr_len(0), vptr_off(0), vptr_free(0){}

    ~djinn_data_container_t() {
        if (vptr_free) {
            delete[] vptr;
        }
    }

    // Owner of this data. Copy it.
    int SetData(uint8_t* data, uint32_t len) {
        if (vptr_off == 0) {
            // std::cerr << "Allocation: " << len+65536 << std::endl;
            vptr = new uint8_t[len+65536];
            vptr_off = len+65536;
            vptr_free = true;
            vptr_len = 0;
        }

        if (len >= vptr_off) {
            // std::cerr << "resizing: " << len+65536 << std::endl;
            if (vptr_free) delete[] vptr;
            vptr = new uint8_t[len+65536];
            vptr_off = len+65536;
            vptr_free = true;
            vptr_len = 0;
        }

        memcpy(vptr, data, len);
        vptr_len = len;
        n = len;
        return 1;
    }

    // Do not own this data.
    int SetDataReference(uint8_t* data, uint32_t len) {
        if (vptr_off != 0) {
            delete[] vptr;
            vptr_free = false;
            vptr_off = 0;
        }
        vptr = data;
        vptr_len = len;
        n = len;
        return 1;
    }

    void reset() {
        n = 0; n_c = 0;
        vptr_len = 0;
    }

    int Serialize(std::ostream& stream) const {
        if (vptr_len) {
            stream.write((char*)&n,   sizeof(int));
            stream.write((char*)&n_c, sizeof(int));
            stream.write((char*)&n_v, sizeof(int));
            stream.write((char*)&vptr_len, sizeof(uint32_t));
            stream.write((char*)vptr, vptr_len);

            return 3*sizeof(int) + sizeof(uint32_t) + vptr_len;
        } else {
            return 0;
        }
    }

    int Serialize(uint8_t* dst) const {
        if (dst == nullptr) return -1;
        int dst_ptr = 0;
        if (vptr_len) {
            memcpy(&dst[dst_ptr], &n,   sizeof(int)); dst_ptr += sizeof(int);
            memcpy(&dst[dst_ptr], &n_c, sizeof(int)); dst_ptr += sizeof(int);
            memcpy(&dst[dst_ptr], &n_v, sizeof(int)); dst_ptr += sizeof(int);
            memcpy(&dst[dst_ptr], &vptr_len, sizeof(uint32_t)); dst_ptr += sizeof(uint32_t);
            memcpy(&dst[dst_ptr], vptr, vptr_len); dst_ptr += vptr_len;
            // std::cerr << n << "," << n_c << "," << n_v << "," << vptr_len << std::endl;

            return 3*sizeof(int) + sizeof(uint32_t) + vptr_len;
        } else {
            return 0;
        }
    }

    int Deserialize(std::istream& stream) {
        if (stream.good() == false) return -1; 
        stream.read((char*)&n,   sizeof(int));
        stream.read((char*)&n_c, sizeof(int));
        stream.read((char*)&n_v, sizeof(int));
        stream.read((char*)&vptr_len, sizeof(uint32_t));
        // std::cerr << n << "," << n_c << "," << n_v << "," << vptr_len << std::endl;

        if (vptr != nullptr) {
            vptr = new uint8_t[vptr_len+65536];
            vptr_free = true;
            vptr_off = vptr_len+65536;
        }

        if (vptr_off <= vptr_len) {
            if (vptr_free) delete[] vptr;
            vptr = new uint8_t[vptr_len+65536];
            vptr_free = true;
            vptr_off = vptr_len+65536;
        }
        
        stream.read((char*)vptr, vptr_len);

        return 3*sizeof(int) + sizeof(uint32_t) + vptr_len;
    }

    int n, n_c, n_v;     // n: size of decompressed data, n_c: size of compressed data, n_v: number of variants
    uint8_t* vptr;       // pointer to data array in bcf1_t->shared.s, excluding the size+type and tag id bytes
    uint32_t vptr_len;   // length of the vptr block or, when set, of the vptr_mod block, excluding offset
    uint32_t vptr_off:31,// vptr offset, i.e., the size of the INFO key plus size+type bytes
             vptr_free:1;// indicates that vptr-vptr_off must be freed; set only when modified and the new
                         //    data block is bigger than the original
};

struct djinn_data_t { // pure virtual base struct
    virtual ~djinn_data_t(){}
    
    virtual uint32_t size() const =0;
    virtual uint16_t GetController() const =0;
    virtual void reset() =0;
    virtual int Serialize(std::ostream& stream) const =0;
    virtual int Serialize(uint8_t* dst) const =0;
    virtual uint32_t SerializedSize() const =0;
    virtual int Deserialize(std::istream& stream, uint16_t ctrl) =0;
};

struct djinn_ctx_t : public djinn_data_t {
public:
    uint32_t size() const override {
        uint32_t tot = 0;

        for (int i = 0; i < 4; ++i) tot += ctx_models[i].vptr_len;

        return tot;
    }

    uint32_t SerializedSize() const override {
        uint32_t tot = 0;

        for (int i = 0; i < 4; ++i) {
            if (ctx_models[i].vptr_len) {
                tot += 3*sizeof(int) + sizeof(uint32_t) + ctx_models[i].vptr_len;
            }
        }

        return tot;
    }

    uint16_t GetController() const override {
        uint16_t raw = 0;
        djinn_ctx_ctrl_t* ctrl = (djinn_ctx_ctrl_t*)&raw;
        ctrl->dimc = ctx_models[0].vptr_len > 0;
        ctrl->dim  = ctx_models[1].vptr_len > 0;
        ctrl->nm   = ctx_models[2].vptr_len > 0;

        return raw;
    }

    void reset() override {
        for (int i = 0; i < 4; ++i) ctx_models[i].reset();
    }

    int Serialize(std::ostream& stream) const override {
        int ret = 0;
        for (int i = 0; i < 4; ++i) ret += ctx_models[i].Serialize(stream);

        return ret;
    }

    int Serialize(uint8_t* dst) const override {
        if (dst == nullptr) return -1;
        int ret = 0;
        for (int i = 0; i < 4; ++i) ret += ctx_models[i].Serialize(&dst[ret]);
        return ret;
    }

    int Deserialize(std::istream& stream, uint16_t ctrl) override {
        if (stream.good() == false) return -1;
        djinn_ctx_ctrl_t* controller = (djinn_ctx_ctrl_t*)&ctrl;

        int ret = 0;
        ctrl >>= 1;
        for (int i = 0; i < 4; ++i) {
            if (ctrl & 1) {
                ret += ctx_models[i].Deserialize(stream);
            }
            ctrl >>= 1;
        }
        return ret;
    }

public:
    djinn_data_container_t ctx_models[4]; // arcetype, 2mc, 2m, nm
};

// TODOOOOO
struct djinn_wah_t : public djinn_data_t {
    uint32_t size() const override {
        uint32_t tot = 0;

        for (int i = 0; i < 6; ++i) tot += wah_models[i].vptr_len;

        return tot;
    }

    uint32_t SerializedSize() const override {
        uint32_t tot = 0;

        for (int i = 0; i < 6; ++i) {
            if (wah_models[i].vptr_len) {
                tot += 3*sizeof(int) + sizeof(uint32_t) + wah_models[i].vptr_len;
            }
        }

        return tot;
    }

    uint16_t GetController() const override {
        uint16_t raw = 0;
        djinn_wah_ctrl_t* ctrl = (djinn_wah_ctrl_t*)&raw;
        ctrl->di2mc1 = wah_models[0].vptr_len > 0;
        ctrl->di2mc2 = wah_models[1].vptr_len > 0;
        ctrl->di2mm1 = wah_models[2].vptr_len > 0;
        ctrl->di2mm2 = wah_models[3].vptr_len > 0;
        ctrl->di2x1  = wah_models[4].vptr_len > 0;
        ctrl->di2x2  = wah_models[5].vptr_len > 0;

        return raw;
    }
    
    void reset() override {
        for (int i = 0; i < 6; ++i) wah_models[i].reset();
    }

    int Serialize(std::ostream& stream) const override { 
        int ret = 0;
        for (int i = 0; i < 6; ++i) ret += wah_models[i].Serialize(stream);

        return ret;
    }

    int Serialize(uint8_t* dst) const override {
        if (dst == nullptr) return -1;
        int ret = 0;
        for (int i = 0; i < 6; ++i) ret += wah_models[i].Serialize(&dst[ret]);
        return ret;
    }

    int Deserialize(std::istream& stream, uint16_t ctrl) override {
        if (stream.good() == false) return -1;
        djinn_wah_ctrl_t* controller = (djinn_wah_ctrl_t*)&ctrl;

        int ret = 0;
        ctrl >>= 1;
        for (int i = 0; i < 6; ++i) {
            if (ctrl & 1) {
                ret += wah_models[i].Deserialize(stream);
            }
            ctrl >>= 1;
        }
        return ret;
    }

public:
    djinn_data_container_t wah_models[6];
};

class djinn_block_t {
public:
    enum class BlockType : uint32_t { UNKNOWN = 0, WAH = 1, CONTEXT = 2 };

    djinn_block_t() : type(BlockType::UNKNOWN), n_rcds(0), p_len(0), ctrl(0), data(nullptr) {}
    virtual ~djinn_block_t() { delete data; }

    virtual int Serialize(std::ostream& stream) const { return -1; }
    virtual int Serialize(uint8_t* dst) const { return -1; }
    virtual int size() const { return -1; }
    virtual int Deserialize(std::istream& stream);

    BlockType type; // 
    size_t n_rcds; // number of records
    size_t p_len; // total compressed length
    uint16_t ctrl; // controller sequence -> recast as djinn_ctx_ctrl_t or djinn_wah_ctrl_t
    djinn_data_t* data; // recast as djinn_ctx_t or djinn_wah_t
};

class djinn_ctx_block_t : public djinn_block_t {
public:
    int Serialize(std::ostream& stream) const override {
        int out_type = (int)type;
        stream.write((char*)&out_type, sizeof(int));
        stream.write((char*)&n_rcds,   sizeof(size_t));
        stream.write((char*)&p_len,    sizeof(size_t));
        stream.write((char*)&ctrl,     sizeof(uint16_t));

        djinn_ctx_t* d = (djinn_ctx_t*)data;
        int ret = d->Serialize(stream);
        assert(ret == d->SerializedSize());
        assert(p_len == d->SerializedSize());

        return sizeof(int) + 2*sizeof(size_t) + sizeof(uint16_t) + ret;
    };

    int Serialize(uint8_t* dst) const override {
        if (dst == nullptr) return -1;
        int dst_ptr = 0;
        int out_type = (int)type;
        memcpy(&dst[dst_ptr], &out_type, sizeof(int));    dst_ptr += sizeof(int);
        memcpy(&dst[dst_ptr], &n_rcds,   sizeof(size_t)); dst_ptr += sizeof(size_t);
        memcpy(&dst[dst_ptr], &p_len,    sizeof(size_t)); dst_ptr += sizeof(size_t);
        memcpy(&dst[dst_ptr], &ctrl,     sizeof(uint16_t)); dst_ptr += sizeof(uint16_t);

        djinn_ctx_t* d = (djinn_ctx_t*)data;
        int ret = d->Serialize(&dst[dst_ptr]);
        dst_ptr += ret;
        assert(ret == d->SerializedSize());
        assert(p_len == d->SerializedSize());

        return dst_ptr;
    };

    int size() const override { 
         int ret = 0;
         djinn_ctx_t* d = (djinn_ctx_t*)data;
         ret += d->size();
         ret += sizeof(int) + 2*sizeof(size_t) + sizeof(uint16_t);
         return ret;
    }

    int Deserialize(std::istream& stream) override;
};

class djinn_wah_block_t : public djinn_block_t {
public:
    int Serialize(std::ostream& stream) const override {
        int out_type = (int)type;
        stream.write((char*)&out_type, sizeof(int));
        stream.write((char*)&n_rcds,   sizeof(size_t));
        stream.write((char*)&p_len,    sizeof(size_t));
        stream.write((char*)&ctrl,     sizeof(uint16_t));

        djinn_wah_t* d = (djinn_wah_t*)data;
        int ret = d->Serialize(stream);
        assert(ret == d->SerializedSize());
        assert(p_len == d->SerializedSize());

        return sizeof(int) + 2*sizeof(size_t) + sizeof(uint16_t) + ret;
    };

    int Serialize(uint8_t* dst) const override {
        if (dst == nullptr) return -1;
        int dst_ptr = 0;
        int out_type = (int)type;
        memcpy(&dst[dst_ptr], &out_type, sizeof(int));    dst_ptr += sizeof(int);
        memcpy(&dst[dst_ptr], &n_rcds,   sizeof(size_t)); dst_ptr += sizeof(size_t);
        memcpy(&dst[dst_ptr], &p_len,    sizeof(size_t)); dst_ptr += sizeof(size_t);
        memcpy(&dst[dst_ptr], &ctrl,     sizeof(uint16_t)); dst_ptr += sizeof(uint16_t);

        djinn_wah_t* d = (djinn_wah_t*)data;
        int ret = d->Serialize(&dst[dst_ptr]);
        dst_ptr += ret;
        assert(ret == d->SerializedSize());
        assert(p_len == d->SerializedSize());

        return dst_ptr;
    };

    int size() const override { 
         int ret = 0;
         djinn_wah_t* d = (djinn_wah_t*)data;
         ret += d->size();
         ret += sizeof(int) + 2*sizeof(size_t) + sizeof(uint16_t);
         return ret;
    }

    int Deserialize(std::istream& stream) override;
};

// djinn_block_t* djinn_create_context_block() {
//     djinn_block_t *djn = (djinn_block_t *)malloc(sizeof(djinn_block_t));
//     if (!djn) return NULL;
//     return djn;
// }

}

#endif