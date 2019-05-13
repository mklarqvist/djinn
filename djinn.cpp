#include "djinn.h"

namespace djinn {

int djinn_ctx_block_t::Deserialize(std::istream& stream) {
    if (stream.good() == false) return -1;
    stream.read((char*)&n_rcds, sizeof(size_t));
    stream.read((char*)&p_len, sizeof(size_t));
    stream.read((char*)&ctrl, sizeof(uint16_t));

    if (data == nullptr) data = new djinn_ctx_t();
    djinn_ctx_t* d = (djinn_ctx_t*)data;
    int ret = d->Deserialize(stream, ctrl);

    return sizeof(int) + 2*sizeof(size_t) + sizeof(uint16_t) + ret;
}

int djinn_wah_block_t::Deserialize(std::istream& stream) {
    if (stream.good() == false) return -1;
    stream.read((char*)&n_rcds, sizeof(size_t));
    stream.read((char*)&p_len, sizeof(size_t));
    stream.read((char*)&ctrl, sizeof(uint16_t));
    std::cerr << std::bitset<16>(ctrl) << std::endl;

    if (data == nullptr) data = new djinn_wah_t();
    djinn_wah_t* d = (djinn_wah_t*)data;
    int ret = d->Deserialize(stream, ctrl);

    return sizeof(int) + 2*sizeof(size_t) + sizeof(uint16_t) + ret;
}

int djinn_block_t::Deserialize(std::istream& stream) { 
    if (stream.good() == false) return -1;
    int out_type = 0;
    stream.read((char*)&out_type, sizeof(int));
    type = BlockType(out_type);

    if (out_type == 0) {
        std::cerr << "illeal crap" << std::endl;
        exit(1);
        return -2;
    }
    else if (out_type == 1) {
        // Illegal
        // static_cast<djinn_wah_block_t*>(this)->Deserialize(stream);

        if (stream.good() == false) return -1;
        stream.read((char*)&n_rcds, sizeof(size_t));
        stream.read((char*)&p_len,  sizeof(size_t));
        stream.read((char*)&ctrl,   sizeof(uint16_t));

        if (data == nullptr) data = new djinn_wah_t();
        djinn_wah_t* d = (djinn_wah_t*)data;
        int ret = d->Deserialize(stream, ctrl);

        return sizeof(int) + 2*sizeof(size_t) + sizeof(uint16_t) + ret;
    } else if (out_type == 2) {
        if (stream.good() == false) return -1;
        stream.read((char*)&n_rcds, sizeof(size_t));
        stream.read((char*)&p_len,  sizeof(size_t));
        stream.read((char*)&ctrl,   sizeof(uint16_t));

        if (data == nullptr) data = new djinn_wah_t();
        djinn_ctx_t* d = (djinn_ctx_t*)data;
        int ret = d->Deserialize(stream, ctrl);
        return sizeof(int) + 2*sizeof(size_t) + sizeof(uint16_t) + ret;
        
    } else {    
        std::cerr << "more illegal crap" << std::endl;
        return -2;
        exit(1);
    }

    return -1; 
}

}