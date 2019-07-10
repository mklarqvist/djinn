[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

![screenshot](DJINN.png)

Djinn is an open-source single-header C++ library for efficiently storing and analyzing large collections of genotypes and/or haplotypes. 

Unlike previous efforts, the Djinn algortihms support multiple alternative alleles (up to 16 by default), mixed phasing, and are embarrassingly parallel because of its blocked memory layout. Note: this library is for genotypes/haplotypes only---there are no included functions to store contig, position, and ref/alt information. The library is currently limited to diploid and haploid
samples.

Implemented algorithms:

- [x] PBWT-preprocessor with RLE-bitmap hybrid (EWAH) compression.
- [x] PBWT-preprocessor with higher-order context modelling.
- [ ] Genotype PBWT (gtPBWT) with RLE-bitmap hybrid compression (from [Tachyon](https://github.com/mklarqvist/tachyon)). 

### Building

For full support, building requires [zstd](https://github.com/facebook/zstd) or [lz4](https://github.com/lz4/lz4). [htslib](https://github.com/samtools/htslib) is required if you intend on using the optional support class `VcfReader.h`. Build with `make` in the root directory.

### Usage



### Note

This is a collaborative effort between Marcus D. R. Klarqvist ([@klarqvist](https://github.com/mklarqvist/)) and James Bonfield ([@jkbonfield](https://github.com/jkbonfield)).

## Performance evaluation

1kgp3-chr20-hg38 (n = 2,548, m = 1,706,442)

| Approach                                | Output (MB) | Compression ratio (VCF) | Import time | Output VCF (time) | Output BCF (time) |
|-----------------------------------------|-------------|-------------------------|-------------|-------------------|-------------------|
| hPBWT + EWAH + LZ4-HC-9                 | 25.36       | 687.81                  | 1m47.015s   | 23.905s           | 22.174s           |
| hPBWT + EWAH + ZSTD-19                  | 22.12       | 788.56                  | 1m52.906s   | 23.838s           | 22.571s           |
| hPBWT + arithmetic                      | 18.49       | 943.37                  | 1m57.831s   | 26.179s           | 25.472s           |
| hpPBWT + arithmetic + no-random-access† | 14.71       | 1185.79                 | 1m42.815s   |                   |                   |
| EWAH + LZ4-HC-9                         | 60.36       | 288.98                  | 1m6.636s    | 10.080s           | 10.176s           |
| EWAH + ZSTD-19                          | 48.48       | 359.8                   | 1m7.568s    | 10.381s           | 10.152s           |
| arithmetic only                         | 80.57       | 216.49                  | 1m57.925s   | 11.852s           | 11.021s           |
| PBWT                                    | 29.81       | 585.14                  | 2m4.855s    | 2m51.207s         | 1m20.485s         |
| BGT*                                    | 64.47       | 270.56                  | 2m20.834s   | CORRUPTED         | CORRUPTED         |
| GTShark                                 | 12.6        | 1384.36                 | 2m23.363s   | 7m55.707s         | 2m19.970s         |
| GQT                                     | 233.7       | 74.64                   | 8m13.750s   | N/A               | N/A               |
| tsinfer‡                                | 184.89      | 94.34                   | 298m51.235s | N/A               | N/A               |
| BCF                                     | 196.88      | 88.6                    | -           | 2m22.953s         | 22.482s (uBCF)    |
| uBCF                                    | 8683.26     | 2.01                    | -           | 2m4.911s          | 2m17.367s         |
| VCF.gz**                                | 240.23      | 72.61                   | -           | 1m10.408s (zcat)  | 7m26.180s         |
| VCF**                                   | 17442.97    | 1                       | -           | -                 | 6m59.067s         |

\* Listing size of the genotype component (`pbf` file); ^ listing size of the genotype component only (`_gt` file).

HRC-chr11 (n = 32,470, m = 1,936,990). 

| Approach                 | Output (MB) | Compression ratio (VCF) | Import time | Output VCF (time)       | Output BCF (time) |
|--------------------------|-------------|-------------------------|-------------|-------------------------|-------------------|
| hPBWT + RLE-H + LZ4-HC-9 | 182.29      | 1380.38                 | 25m56.414s  | 6m22.417s               | 6m14.151s         |
| hPBWT + RLE-H + ZSTD-19  | 159         | 1582.58                 | 26m37.005s  | 6m22.953s               | 6m18.962s         |
| hPBWT + arithmetic       | 123.19      | 2042.61                 | 41m16.079s  | 10m41.234s              | 10m32.878s        |
| WAH + LZ4-HC-9           | 1087.25     | 231.44                  | 16m26.147s  | 2m24.646s               | 2m21.152s         |
| WAH + ZSTD-19            | 656.04      | 383.56                  | 31m27.405s  | 2m27.597s               | 2m25.958s         |
| arithmetic only          | 1923.25     | 130.84                  | 41m3.145s   | -                       | -                 |
| BGT*                     | 365.65      | 688.17                  | 33m10.116s  | 25m11.121s              | 12m5.086s         |
| GTShark^                 | 89.78       | 2802.74                 | 33m38.250s  | -                       | 36m10.505s        |
| GQT                      | 4877.69     | 51.59                   | 195m16.803s | N/A                     | N/A               |
| BCF                      | 3476.07     | 72.39                   | -           | 36m15.195s              | 5m16.029s (uBCF)  |
| uBCF                     | 125549.35   | 2                       | -           | 31m56.717s              | 41m22.754s        |
| VCF.gz**                 | 3940.51     | 63.86                   | -           | 104m19.148s             | -                 |
| VCF**                    | 251629.59   | 1                       | -           | 95m48.053s (14m28.331s) | -                 |

More 1kgp3-hg38

| Chr | hPBWT + arithmetic | BCF    | Ratio |
|-----|--------------------------|--------|-------|
| 1   | 58.25                    | 773.03 | 13.27 |
| 2   | 62.62                    | 834.94 | 13.33 |
| 3   | 52.05                    | 710.43 | 13.65 |
| 4   | 48.98                    | 711.16 | 14.52 |
| 5   | 45.57                    | 640.92 | 14.06 |
| 6   | 43.83                    | 640.84 | 14.62 |
| 7   | 41.84                    | 580.63 | 13.88 |
| 8   | 39.74                    | 554.35 | 13.95 |
| 9   | 32.43                    | 430.84 | 13.29 |
| 10  | 36.2                     | 501.42 | 13.85 |
| 11  | 35.23                    | 491.32 | 13.95 |
| 12  | 34.49                    | 475.37 | 13.78 |
| 13  | 25.47                    | 358.16 | 14.06 |
| 14  | 23.63                    | 323.14 | 13.67 |
| 15  | 22.48                    | 290.02 | 12.9  |
| 16  | 24.96                    | 314.68 | 12.61 |
| 17  | 21.79                    | 272.88 | 12.52 |
| 18  | 20.92                    | 281.59 | 13.46 |
| 19  | 17.99                    | 227.26 | 12.63 |
| 20  | 18.5                     | 229.9  | 12.43 |
| 21  | 11.07                    | 139.3  | 12.58 |
| 22  | 11.3                     | 138.28 | 12.24 |

### Stripped vcf files

To make comparisons as fair as possible to other tools we created GT only BCF files:

```bash
bcftools view /media/mdrk/NVMe/1kgp3/bcf/hg38/1kgp3_chr20_hg38.bcf -h > /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.vcf
bcftools view /media/mdrk/NVMe/1kgp3/bcf/hg38/1kgp3_chr20_hg38.bcf -m2 -M2 | bcftools query -f "%CHROM\t%POS\t.\t%REF\t%ALT\t.\t.\t.\tGT[\t%GT]\n" >> /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.vcf
bcftools view /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.vcf -O b -o /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.bcf --threads 28
time bgzip --threads 28 -c /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.vcf > /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.vcf.gz
time bcftools view /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.bcf -O u -o /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.ubcf --threads 28
```

Extract genotypes only:

```bash
time bcftools query -f "[\t%GT]\n" /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare.bcf | bgzip --threads 28 -c > /media/mdrk/NVMe/1kgp3/hg38_1kgp3_chr20_bare_gt.vcf.gz
```

For diploid and biallelic sites the memory cost for genotypes in uncompressed bcf is NM bytes and M(4N-1) bytes in raw vcf.

## Algorithm overview

## Getting help

We are actively developing Djinn and are always interested in improving its quality. If you run into an issue, please report the problem on our [Issue tracker](https://github.com/mklarqvist/djinn/issues). Be sure to add enough detail to your report that we can reproduce the problem and address it.

## Developers' Guide

Djinn is designed as a programming library providing simple C++ APIs to build/load an archieve and to query it. The `examples/` directory contains a variety of examples that demonstrates typical uses of the C++ APIs. The library requires only a single header: `djinn.h`. This file contains additional detailed API documentation. 
Djinn aims to keep the APIs for `djinn_*` classes in `djinn.h` stable.

## Limitations

* Djinn is a framework for storing and querying sequence variant data only. This excludes all other meta data such as positional information (contig and position) and reference/alternative allele encodings.