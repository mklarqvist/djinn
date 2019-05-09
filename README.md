[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Docs](https://img.shields.io/badge/Docs-Available-brightgreen.svg)](https://mklarqvist.github.io/tomahawk/)

![screenshot](DJINN.png)

Building requires the packages: htslib, zstd, lz4, and openssl if debugging.

## Performance evaluation

| Approach              | Output (MB) | Compression ratio (VCF) | Import time | Output VCF (time) | Output BCF (time) |
|-----------------------|-------------|-------------------------|-------------|-------------------|-------------------|
| hPBWT + WAH + LZ4-HC9 | 31.62       | 551.64                  | 2m12.883s   | 26.104s           | 30.441s           |
| hPBWT + WAH + ZSTD-19 | 26.91       | 648.2                   | 2m16.400s   | 0m27.784s         | 0m30.560s         |
| hPBWT + context       | 18.45       | 945.42                  | 2m32.443s   | 1                 | 1                 |
| BGT*                  | 64.47       | 270.56                  | 2m20.834s   | CORRUPTED         | CORRUPTED         |
| GQT                   | 233.7       | 74.64                   | 8m13.750s   | N/A               | N/A               |
| BCF                   | 196.88      | 88.6                    | -           | 2m22.953s         | -                 |
| uBCF                  | 8683.26     | 2.01                    | -           | 2m4.911s          | 2m17.367s         |
| VCF.gz**              | 240.23      | 72.61                   | -           | 1m10.408s (zcat)  | 7m26.180s         |
| VCF**                 | 17442.97    | 1                       | -           | 0m2.859s (cat)    | 6m59.067s         |

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

On average, genotypic memory cost for bcf is NM bytes and M(4N-1) bytes for vcf.