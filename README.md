[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

![screenshot](DJINN.png)

Djinn is an open-source C++ library for efficiently storing and analyzing whole-genome genotypes and/or haplotypes. Unlike previous efforts, the Djinn algortihms support multiple alternative alleles (up to 16 by default), mixed phasing, and are embarrassingly parallel because of its blocked memory layout. Note: this library is for genotypes/haplotypes only---there are no included functions to store contig, position, and ref/alt information. The library is currently limited to diploid and haploid
samples.

Implemented algorithms:

- [x] Haplotype PBWT (hPWBT) with RLE-bitmap hybrid (RLE-H) compression.
- [x] Haplotype PBWT with higher-order context modelling.
- [ ] Genotype PBWT (gtPBWT) with RLE-bitmap hybrid compression (from [Tachyon](https://github.com/mklarqvist/tachyon)). 
- [x] RLE-bitmap hybrid compression for use with the gtOcc data structure.
- [ ] Block-wise genotype PBWT with word-aligned hybrid compression.

### Building

Building requires the packages: [htslib](https://github.com/samtools/htslib), [zstd](https://github.com/facebook/zstd), [lz4](https://github.com/lz4/lz4), and `openssl` if debugging. Build with `make` and build with debug flags by passing invoking `make debug`. Compile with `make debug_size` to get compression performance dumped to the console.

### Note

This is a collaborative effort between Marcus D. R. Klarqvist ([@klarqvist](https://github.com/mklarqvist/)) and James Bonfield ([@jkbonfield](https://github.com/jkbonfield)).

## Performance evaluation

1kgp3-chr20-hg38 (n = 2548, m = 1706442)

| Approach                 | Output (MB) | Compression ratio (VCF) | Import time | Output VCF (time) | Output BCF (time) |
|--------------------------|-------------|-------------------------|-------------|-------------------|-------------------|
| hPBWT + RLE-H + LZ4-HC-9 | 32.4        | 538.36                  | 2m12.883s   | 23.905s           | 22.174s           |
| hPBWT + RLE-H + ZSTD-19  | 27.79       | 627.67                  | 2m16.400s   | 23.838s           | 22.571s           |
| hPBWT + arithmetic       | 18.49       | 943.37                  | 1m57.831s   | 26.179s           | 25.472s           |
| WAH + LZ4-HC-9           | 60.36       | 288.98                  | 1m6.636s    | 10.080s           | 10.176s           |
| WAH + ZSTD-19            | 48.48       | 359.8                   | 1m7.568s    | 10.381s           | 10.152s           |
| arithmetic only          | 80.57       | 216.49                  | 1m57.925s   | 11.852s           | 11.021s           |
| BGT*                     | 64.47       | 270.56                  | 2m20.834s   | CORRUPTED         | CORRUPTED         |
| GTShark^                  | 12.6        | 1384.36                 | 2m23.363s   | 7m55.707s         | 2m19.970s         |
| GQT                      | 233.7       | 74.64                   | 8m13.750s   | N/A               | N/A               |
| BCF                      | 196.88      | 88.6                    | -           | 2m22.953s         | 22.482s (uBCF)    |
| uBCF                     | 8683.26     | 2.01                    | -           | 2m4.911s          | 2m17.367s         |
| VCF.gz**                 | 240.23      | 72.61                   | -           | 1m10.408s (zcat)  | 7m26.180s         |
| VCF**                    | 17442.97    | 1                       | -           | -                 | 6m59.067s         |

\* Listing size of the genotype component (`pbf` file); ^ listing size of the genotype component only (`_gt` file).

HRC-chr11 (n = 32470, m = 1936990)

| Approach                 | Output (MB) | Compression ratio (VCF) | Import time | Output VCF (time)       | Output BCF (time) |
|--------------------------|-------------|-------------------------|-------------|-------------------------|-------------------|
| hPBWT + RLE-H + LZ4-HC-9 | 182.29      | 1380.38                 | 25m56.414s  | 6m22.417s               | 6m14.151s         |
| hPBWT + RLE-H + ZSTD-19  | 159         | 1582.58                 | 26m37.005s  | 6m22.953s               | 6m18.962s         |
| hPBWT + arithmetic       | 123.19      | 2042.61                 | 41m16.079s  | 10m41.234s              | 10m32.878s        |
| WAH + LZ4-HC-9           | 1087.25     | 231.44                  | 16m26.147s  | 2m24.646s               | 2m21.152s         |
| WAH + ZSTD-19            | 656.04      | 383.56                  | 31m27.405s  | 2m27.597s               | 2m25.958s         |
| arithmetic only          | 1923.25     | 130.84                  | 41m3.145s   | -                       | -                 |
| BGT*                     | 365.65      | 688.17                  | 33m10.116s  | 25m11.121s              | 12m5.086s         |
| GTShark^                  | 89.78       | 2802.74                 | 33m38.250s  | -                       | 36m10.505s        |
| GQT                      | 4877.69     | 51.59                   | 195m16.803s | N/A                     | N/A               |
| BCF                      | 3476.07     | 72.39                   | -           | 36m15.195s              | 5m16.029s (uBCF)  |
| uBCF                     | 125549.35   | 2                       | -           | 31m56.717s              | 41m22.754s        |
| VCF.gz**                 | 3940.51     | 63.86                   | -           | 104m19.148s             | -                 |
| VCF**                    | 251629.59   | 1                       | -           | 95m48.053s (14m28.331s) | -                 |

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