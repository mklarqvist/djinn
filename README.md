<p align="center"><img src="https://raw.githubusercontent.com/mklarqvist/djinn/auto/DJINN.png" alt="Djinn" width="350px"></p>

Djinn is an open-source single-header C++ library for efficiently storing and analyzing large collections of genotypes and/or haplotypes.

### Branch status:

[![Build Status](https://travis-ci.com/mklarqvist/djinn.svg?branch=master)](https://travis-ci.com/mklarqvist/djinn)
[![Release](https://img.shields.io/badge/Release-beta_0.0.1-blue.svg)](https://github.com/mklarqvist/djinn/releases)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

## Performance evaluation

For reference, our proposed compression algorithms were tested and compared against each other and against the incumbent standard Vcf and Bcf formats using a disk-based benchmark. The host architecture used is a 10 nm Cannon Lake Core i3-8121U with gcc (GCC) 7.3.1 20180303 (Red Hat 7.3.1-5).

| Method                   | File size (MB) | Comp. Time (s) | Ratio (VCF) | Ratio (BCF) | Decompress* (MB/s) | Inflate** (MB/s) | Output VCF (s) |
|--------------------------|----------------|----------------|------------------|------------------|-------------------|----------------|----------------|
| djinn PBWT+EWAH+ZSTD-19  | 21.36          | 119.9        | 814.2          | 10.7           | 24119           | 1140        | 31.0         |
| djinn PBWT+EWAH+LZ4-HC-9 | 26.48          | 93.9         | 656.7          | 8.6            | **25917**           | 1143        | 28.2         |
| djinn PBWT+CTX           | **16.19**          | **83.0**         | **1074.1**          | **14.2**             | 3757           | 907         | 34.1         |
| djinn EWAH+LZ4-HC-9      | 57.05          | 84.9         | 304.8          | 4.0             | 15241           | **2827**         | **19.0**          |
| djinn EWAH+ZSTD-19       | 41.84          | 207.4        | 415.6          | 5.4            | 12853           | 2667        | 19.9         |
| djinn CTX                | 86.7           | 73.6         | 200.5          | 2.6            | 1118           | 858        | 33.2         |
| BCF                      | 229.9          | NA             | 80.5            | 1                | NA                | NA             | 148.3        |
| VCF.gz                   | 281.34         | NA             | 65.8            | 0.82             | NA                | NA             | 457.6         |

\*Decompress: decompressing into the succinct EWAH data structure. \*\*Inflate: decompress and inflate into byte arrays of length N (as in uBcf).

Djinn can offer strong compression with slower decompression speeds using statistical models (CTX). Reversely, the EWAH-based algorithms provide stronger decompression and query speeds with lower compression rates. The optimal algorithm dependends on its application context: whether query speeds or compression matters the most.

### Building

Building Djinn requires no additional packages. For full support, building requires [zstd](https://github.com/facebook/zstd) or [lz4](https://github.com/lz4/lz4). If you intend on using the optional support class `vcf_reader.h` for consuming `htslib`-based files then [htslib](https://github.com/samtools/htslib) is required. Build with autotools:

```bash
./autogen.sh   # Generate the header template
./configure    # Generate Makefile
make           # Build libraries and place in .libs directory
make install   # Install headers and libraries
```

## Developers' Guide

Djinn is designed as a programming library providing simple C++ APIs to build/load an archieve and to query it. The `examples/` directory contains a variety of examples that demonstrates typical uses of the C++ APIs. The library requires only a single header: `djinn.h`. This file contains additional detailed API documentation. 
Djinn aims to keep the APIs for `djinn_*` classes in `djinn.h` stable.

Implemented algorithms:

- [x] PBWT-preprocessor with RLE-bitmap hybrid (EWAH) compression.
- [x] PBWT-preprocessor with higher-order context modelling.
- [ ] Genotype PBWT (gtPBWT) with RLE-bitmap hybrid compression (from [Tachyon](https://github.com/mklarqvist/tachyon)).

### Quick guide

A simple quick start guide to encoding:

```C++
#include <djinn/djinn.h>

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
#include <djinn/djinn.h>

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

## Algorithmic overview

The extended word-aligned hybrid (EWAH) data model was first described by [@lemire](https://github.com/lemire/):

- Lemire D, Kaser O, Aouiche K. Sorting improves word-aligned bitmap indexes.Data & Knowledge Engineering2010;69(1):3–28, doi:[10.1016/j.datak.2009.08.006](https://arxiv.org/abs/0901.3751).

The positional Burrows-Wheeler transform (PBWT) was described by [@richarddurbin](https://github.com/richarddurbin/):

- Richard Durbin, Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (PBWT), Bioinformatics, Volume 30, Issue 9, 1 May 2014, Pages 1266–1272, [https://doi.org/10.1093/bioinformatics/btu014](https://doi.org/10.1093/bioinformatics/btu014)

The idea of using statistical models in compression has been around for decades but the ideas and original code used in this project was first described by [@jkbonfield](https://github.com/jkbonfield):

- Bonfield JK, Mahoney MV (2013) Compression of FASTQ and SAM Format Sequencing Data. PLoS ONE 8(3): e59190. [https://doi.org/10.1371/journal.pone.0059190](https://doi.org/10.1371/journal.pone.0059190)

## Limitations

- Djinn is a framework for storing and querying sequence variant data only. This excludes all other meta data such as positional information (contig and position) and reference/alternative allele encodings.
- Djinn is currently limited to at most 16 alternative alleles.
- There is no default support for random access in a data block. If desired, this must be implemented by the user. For EWAH-enoded data this is trivial: store the virtual byte offset to the beginning of each variant site in an array.
  
## Note

This is a collaborative effort between Marcus D. R. Klarqvist ([@klarqvist](https://github.com/mklarqvist/)) and James Bonfield ([@jkbonfield](https://github.com/jkbonfield)).

## License

Djinn is licensed under [MIT](LICENSE)
