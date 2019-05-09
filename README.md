[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

![screenshot](DJINN.png)

Building requires the packages: [htslib](https://github.com/samtools/htslib), [zstd](https://github.com/facebook/zstd), [lz4](https://github.com/lz4/lz4), and `openssl` if debugging.

### Note

This is a collaborative effort between Marcus D. R. Klarqvist ([@klarqvist](https://github.com/mklarqvist/)) and James Bonfield ([@jkbonfield](https://github.com/jkbonfield)).

## Performance evaluation

1KGP3-chr20-hg38 (n=2548)

| Approach              | Output (MB) | Compression ratio (VCF) | Import time | Output VCF (time) | Output BCF (time) |
|-----------------------|-------------|-------------------------|-------------|-------------------|-------------------|
| hPBWT + WAH + LZ4-HC9 | 31.62       | 551.64                  | 2m12.883s   | 26.104s           | 30.441s           |
| hPBWT + WAH + ZSTD-19 | 26.91       | 648.2                   | 2m16.400s   | 27.784s           | 30.560s           |
| hPBWT + context       | 18.45       | 945.42                  | 2m32.443s   | 1                 | 1                 |
| BGT*                  | 64.47       | 270.56                  | 2m20.834s   | CORRUPTED         | CORRUPTED         |
| GQT                   | 233.7       | 74.64                   | 8m13.750s   | N/A               | N/A               |
| BCF                   | 196.88      | 88.6                    | -           | 2m22.953s         | -                 |
| uBCF                  | 8683.26     | 2.01                    | -           | 2m4.911s          | 2m17.367s         |
| VCF.gz**              | 240.23      | 72.61                   | -           | 1m10.408s (zcat)  | 7m26.180s         |
| VCF**                 | 17442.97    | 1                       | -           | -                 | 6m59.067s         |

HRC-chr11 (n = 32470, m = 1936990)

| Approach              | Output (MB) | Compression ratio (VCF) | Import time | Output VCF (time) | Output BCF (time) |
|-----------------------|-------------|-------------------------|-------------|-------------------|-------------------|
| hPBWT + WAH + LZ4-HC9 | 181.26      | 1388.22                 | 32m2.310s   | 7m16.925s         | 7m27.964s         |
| hPBWT + WAH + ZSTD-19 | 156         | 1613.01                 | 33m0.165s   |                   |                   |
| hPBWT + context       | 123.19      | 2042.61                 |             |                   |                   |
| BGT*                  | 365.65      | 688.17                  | 33m10.116s  | 25m11.121s        | 12m5.086s         |
| GQT                   |             |                  |             | N/A               | N/A               |
| BCF                   | 3476.07     | 72.39                   |             |                   |                   |
| uBCF                  |             |                  |             |                   |                   |
| VCF.gz**              |             |                  |             |                   |                   |
| VCF**                 | 251629.59   | 1                       |             |                   |                   |

### Stripped vcf files

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