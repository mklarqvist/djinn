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

struct djinn_hdr_t {
    uint8_t version[3];
    uint8_t base_ploidy;
    int64_t n_samples; // number of samples
};

struct djinn_ctx_ctrl_t {
    uint16_t pbwt: 1, 
             di2mc1: 1, 
             di2mc2: 1, 
             di2mm1: 1, 
             di2mm2: 1, 
             di2x1:  1, 
             di2x2:  1, 
             di2mc_bit1: 1, 
             di2mc_bit2: 1, 
             unused: 7;
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

#define DJINN_CTX_MODEL_DI2MC1      0
#define DJINN_CTX_MODEL_DI2MC2      1
#define DJINN_CTX_MODEL_DI2MM1      2
#define DJINN_CTX_MODEL_DI2MM2      3
#define DJINN_CTX_MODEL_DI2X1       4
#define DJINN_CTX_MODEL_DI2X2       5
#define DJINN_CTX_MODEL_DI2MC_PART1 0
#define DJINN_CTX_MODEL_DI2MC_PART2 1

#define DJINN_WAH_MODEL_DI2MC1      0
#define DJINN_WAH_MODEL_DI2MC2      1
#define DJINN_WAH_MODEL_DI2MM1      2
#define DJINN_WAH_MODEL_DI2MM2      3
#define DJINN_WAH_MODEL_DI2X1       4
#define DJINN_WAH_MODEL_DI2X2       5

struct djinn_data_t { // pure virtual base struct
    virtual uint32_t size() const =0;
    virtual uint16_t GetController() const =0;
    virtual void reset() =0;
    virtual int Serialize(std::ostream& stream) const =0;
    virtual int Serialize(uint8_t* dst) const =0;
    virtual uint32_t SerializedSize() const =0;
    virtual int Deserialize(std::istream& stream, uint16_t ctrl) =0;
};

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

    int n, n_c, n_v;          // n: size of decompressed data, n_c: size of compressed data, n_v: number of variants
    uint8_t* vptr;       // pointer to data array in bcf1_t->shared.s, excluding the size+type and tag id bytes
    uint32_t vptr_len;   // length of the vptr block or, when set, of the vptr_mod block, excluding offset
    uint32_t vptr_off:31,// vptr offset, i.e., the size of the INFO key plus size+type bytes
             vptr_free:1;// indicates that vptr-vptr_off must be freed; set only when modified and the new
                         //    data block is bigger than the original
};

struct djinn_ctx_t : public djinn_data_t {
public:
    uint32_t size() const override {
        uint32_t tot = 0;

        for (int i = 0; i < 6; ++i) tot += ctx_models[i].vptr_len;
        for (int i = 0; i < 2; ++i) tot += ctx_partitions[i].vptr_len;

        return tot;
    }

    uint32_t SerializedSize() const override {
        uint32_t tot = 0;

        for (int i = 0; i < 6; ++i) {
            if (ctx_models[i].vptr_len) {
                tot += 3*sizeof(int) + sizeof(uint32_t) + ctx_models[i].vptr_len;
            }
        }
        for (int i = 0; i < 2; ++i) {
            if (ctx_partitions[i].vptr_len) {
                tot += 3*sizeof(int) + sizeof(uint32_t) + ctx_partitions[i].vptr_len;
            }
        }

        return tot;
    }

    uint16_t GetController() const override {
        uint16_t raw = 0;
        djinn_ctx_ctrl_t* ctrl = (djinn_ctx_ctrl_t*)&raw;
        ctrl->di2mc1 = ctx_models[0].vptr_len > 0;
        ctrl->di2mc2 = ctx_models[1].vptr_len > 0;
        ctrl->di2mm1 = ctx_models[2].vptr_len > 0;
        ctrl->di2mm2 = ctx_models[3].vptr_len > 0;
        ctrl->di2x1  = ctx_models[4].vptr_len > 0;
        ctrl->di2x2  = ctx_models[5].vptr_len > 0;
        ctrl->di2mc_bit1 = ctx_partitions[0].vptr_len > 0;
        ctrl->di2mc_bit2 = ctx_partitions[1].vptr_len > 0;

        return raw;
    }

    void reset() {
        for (int i = 0; i < 6; ++i) ctx_models[i].reset();
        for (int i = 0; i < 2; ++i) ctx_partitions[i].reset();
    }

    int Serialize(std::ostream& stream) const override {
        int ret = 0;
        for (int i = 0; i < 6; ++i) ret += ctx_models[i].Serialize(stream);
        for (int i = 0; i < 2; ++i) ret += ctx_partitions[i].Serialize(stream);

        return ret;
    }

    int Serialize(uint8_t* dst) const override {
        if (dst == nullptr) return -1;
        int ret = 0;
        for (int i = 0; i < 6; ++i) ret += ctx_models[i].Serialize(&dst[ret]);
        for (int i = 0; i < 2; ++i) ret += ctx_partitions[i].Serialize(&dst[ret]);
        return ret;
    }

    int Deserialize(std::istream& stream, uint16_t ctrl) override {
        if (stream.good() == false) return -1;
        djinn_ctx_ctrl_t* controller = (djinn_ctx_ctrl_t*)&ctrl;

        int ret = 0;
        ctrl >>= 1;
        for (int i = 0; i < 6; ++i) {
            if (ctrl & 1) {
                ret += ctx_models[i].Deserialize(stream);
            }
            ctrl >>= 1;
        }
        for (int i = 0; i < 2; ++i) {
            if (ctrl & 1) {
                ret += ctx_partitions[i].Deserialize(stream);
            }
            ctrl >>= 1;
        }
        return ret;
    }

public:
    djinn_data_container_t ctx_models[6];
    djinn_data_container_t ctx_partitions[2];
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
    
    void reset() {
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

    djinn_block_t() : type(BlockType::UNKNOWN), n_rcds(0), p_len(0), ctrl(0), data(nullptr){}
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