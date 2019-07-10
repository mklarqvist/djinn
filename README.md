[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

![screenshot](DJINN.png)

Djinn is an open-source single-header C++ library for efficiently storing and analyzing large collections of genotypes and/or haplotypes. 

Implemented algorithms:

- [x] PBWT-preprocessor with RLE-bitmap hybrid (EWAH) compression.
- [x] PBWT-preprocessor with higher-order context modelling.
- [ ] Genotype PBWT (gtPBWT) with RLE-bitmap hybrid compression (from [Tachyon](https://github.com/mklarqvist/tachyon)). 

## Performance evaluation

For reference, our proposed compression algorithms were tested and compared against each other and against the incumbent standard Vcf and Bcf formats using a disk-based benchmark. The host architecture used is a 10 nm Cannon Lake Core i3-8121U with gcc (GCC) 7.3.1 20180303 (Red Hat 7.3.1-5).

| Method                   | File size (MB) | Comp. Time (s) | Ratio (VCF) | Ratio (BCF) | Decompress (MB/s) | Inflate (MB/s) | Output VCF (s) |
|--------------------------|----------------|----------------|------------------|------------------|-------------------|----------------|----------------|
| djinn PBWT+EWAH+ZSTD-19  | 21.36          | 119.9        | 814.2          | 10.7           | 24119           | 1140        | 31.0         |
| djinn PBWT+EWAH+LZ4-HC-9 | 26.48          | 93.9         | 656.7          | 8.6            | 25917           | 1143        | 28.2         |
| djinn PBWT+CTX           | 16.19          | 83.0         | 1074.1          | 14.2             | 3757           | 907         | 34.1         |
| djinn EWAH+LZ4-HC-9      | 57.05          | 84.9         | 304.8          | 4.0             | 15241           | 2827         | 19.0          |
| djinn EWAH+ZSTD-19       | 41.84          | 207.4        | 415.6          | 5.4            | 12853           | 2667        | 19.9         |
| djinn CTX                | 86.7           | 73.6         | 200.5          | 2.6            | 1118           | 858        | 33.2         |
| BCF                      | 229.9          | NA             | 80.5            | 1                | NA                | NA             | 148.3        |
| VCF.gz                   | 281.34         | NA             | 65.8            | 0.82             | NA                | NA             | 457.6         |

Djinn can offer strong compression with slower decompression speeds using statistical models (CTX). Reversely, the EWAH-based algorithms provide stronger decompression and query speeds with lower compression rates. The optimal algorithm dependends on its application context: whether query speeds or compression matters the most.

### Building

For full support, building requires [zstd](https://github.com/facebook/zstd) or [lz4](https://github.com/lz4/lz4). If you intend on using the optional support class `VcfReader.h` then [htslib](https://github.com/samtools/htslib) is required. Build with `make` in the root directory.


## Developers' Guide

Djinn is designed as a programming library providing simple C++ APIs to build/load an archieve and to query it. The `examples/` directory contains a variety of examples that demonstrates typical uses of the C++ APIs. The library requires only a single header: `djinn.h`. This file contains additional detailed API documentation. 
Djinn aims to keep the APIs for `djinn_*` classes in `djinn.h` stable.

### Quick guide

A simple quick start guide to encoding:
```C++
#include <djinn.h>

// Initiate new model of interest.
djinn::djinn_model* djn_ctx = new djinn::djinn_ctx_model();

// Encode Bcf-encoded data
int ret = djn_ctx->EncodeBcf(data_in, data_len, ploidy, n_allele);

// Finish encoding and serialize to stdcout
djn_ctx->FinishEncoding();
int serial_size = djn_ctx->Serialize(std::cout);
delete djn_ctx
```

Changing the desired algorithm is easy: replace `djinn::djinn_ctx_model()` with for example,
 `djinn::djinn_ewah_model(djinn::CompressionStrategy::LZ4, 9)` or `djinn::djinn_ewah_model(djinn::CompressionStrategy::ZSTD, 21)`.

Decoding data is equally simple

```C++
#include <djinn.h>

// Initiate new model of interest.
djinn::djinn_model* djn_decode = new djinn::djinn_ctx_model();

// Deserialize model from cin
int decode_ctx_ret = djn_decode->Deserialize(std::cin);
// Primary return structure in Djinn
djinn::djinn_variant_t* variant = nullptr;

// Start decoding and retrieve variants iteratively
djn_ctx->StartDecoding();
for (int i = 0; i < djn_decode->n_variants; ++i) {
    int objs = djn_decode->DecodeNext(variant);
    assert(objs > 0);
}
delete variant;
delete djn_decode
```



## Getting help

We are actively developing Djinn and are always interested in improving its quality. If you run into an issue, please report the problem on our [Issue tracker](https://github.com/mklarqvist/djinn/issues). Be sure to add enough detail to your report that we can reproduce the problem and address it.

## Limitations

* Djinn is a framework for storing and querying sequence variant data only. This excludes all other meta data such as positional information (contig and position) and reference/alternative allele encodings.
* Djinn is currently limited to at most 16 alternative alleles.
  
### Note

This is a collaborative effort between Marcus D. R. Klarqvist ([@klarqvist](https://github.com/mklarqvist/)) and James Bonfield ([@jkbonfield](https://github.com/jkbonfield)).