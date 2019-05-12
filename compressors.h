#include <zstd.h>
#include <zstd_errors.h>
#include "lz4.h" // lz4
#include "lz4hc.h"

namespace djinn {

static
int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
}

static
int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int compressed_data_size = ZSTD_decompress(out, out_capacity, in, n_in);
    if (compressed_data_size < 0) {
        std::cerr << "A negative result from ZSTD_decompress indicates a failure trying to compress the data.  See exit code (echo $?) for value returned." << std::endl;
        exit(1);
    }
    return(compressed_data_size);
}

static
int Lz4Compress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int compressed_data_size = LZ4_compress_HC((const char*)in, (char*)out, n_in, out_capacity, c_level);

    if (compressed_data_size < 0) {
        std::cerr << "A negative result from LZ4_compress_default indicates a failure trying to compress the data.  See exit code (echo $?) for value returned." << std::endl;
        exit(1);
    }
        
    if (compressed_data_size == 0) {
        std::cerr << "A result of 0 means compression worked, but was stopped because the destination buffer couldn't hold all the information." << std::endl;
        exit(1);
    }


    return(compressed_data_size);
}

static
int Lz4Decompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int32_t decompressed_size = LZ4_decompress_safe((char*)in, (char*)out, n_in, out_capacity);
    // int decompressed_data_size = LZ4_decompress((const char*)in, (char*)out, n_in, out_capacity);
    if (decompressed_size < 0) {
        std::cerr << "A negative result from LZ4_decompress_safe indicates a failure trying to decompress the data.  See exit code (echo $?) for value returned." << std::endl;
        exit(1);
    }
    if (decompressed_size == 0) {
        std::cerr << "I'm not sure this function can ever return 0.  Documentation in lz4.h doesn't indicate so.";
        exit(1);
    }


    return(decompressed_size);
}

}