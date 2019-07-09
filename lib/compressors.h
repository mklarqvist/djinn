#if defined HAVE_ZSTD
#include <zstd.h>
#include <zstd_errors.h>
#endif

#if defined HAVE_LZ4
#include <lz4.h> // lz4
#include <lz4hc.h> // lz4hc
#endif

namespace djinn {

typedef int (*GeneralCompressor)(const uint8_t*, uint32_t, uint8_t*, uint32_t, const int32_t);
typedef int (*GeneralDecompressor)(const uint8_t*, uint32_t, uint8_t*, uint32_t);

static
int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
#if defined HAVE_ZSTD
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
#else
    std::cerr << "Program was not compiled with Zstd support!" << std::endl;
    return -1;
#endif
}

static
int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
#if defined HAVE_ZSTD
    int compressed_data_size = ZSTD_decompress(out, out_capacity, in, n_in);
    
    if (compressed_data_size < 0) {
        std::cerr << "A negative result from ZSTD_decompress indicates a failure trying to compress the data.  See exit code (echo $?) for value returned." << std::endl;
        exit(1);
    }

    return(compressed_data_size);
#else
    std::cerr << "Program was not compiled with Zstd support!" << std::endl;
    return -1;
#endif
}

static
int Lz4Compress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
#if defined HAVE_LZ4
    int compressed_data_size = LZ4_compress_HC((const char*)in, (char*)out, n_in, out_capacity, c_level);

    if (compressed_data_size < 0) {
        std::cerr << "A negative result from LZ4_compress_default indicates a failure trying to compress the data.  See exit code (echo $?) for value returned." << std::endl;
        exit(compressed_data_size);
    }
        
    if (compressed_data_size == 0) {
        std::cerr << "A result of 0 means compression worked, but was stopped because the destination buffer couldn't hold all the information." << std::endl;
        exit(compressed_data_size);
    }

    return(compressed_data_size);
#else
    std::cerr << "Program was not compiled with LZ4 support!" << std::endl;
    return -1;
#endif
}

static
int Lz4Decompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
#if defined HAVE_LZ4
    int32_t decompressed_size = LZ4_decompress_safe((const char*)in, (char*)out, n_in, out_capacity);
    
    if (decompressed_size < 0) {
        std::cerr << "A negative result from LZ4_decompress_safe indicates a failure trying to decompress the data.  See exit code (echo $?) for value returned." << std::endl;
        exit(decompressed_size);
    }

    if (decompressed_size == 0) {
        std::cerr << "I'm not sure this function can ever return 0. Documentation in lz4.h doesn't indicate so.";
        exit(decompressed_size);
    }

    return(decompressed_size);
#else
    std::cerr << "Program was not compiled with LZ4 support!" << std::endl;
    return -1;
#endif
}

}