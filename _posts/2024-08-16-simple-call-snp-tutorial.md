---
layout: post
title:  "简单的 Call SNP 教程"
categories: 教程
---

```bash
bcftools mpileup \
    -O u -a AD,DP \
    -r {wildcards.chr} \
    -f {input.ref_genome} \
    {input.bam_list} \
    2> {log} \
| bcftools call \
    -v -m -O z -a GQ \
    -o {output.chr_vcf} \
    2>> {log}
```

```bash
bcftools concat \
    {input.chr_vcf_list} \
    2> {log} \
| bcftools norm \
    -m-both \
    -f {input.ref_genome} \
    2>> {log} \
| bgzip \
    -c \
    > {output.norm_vcf} \
    2>> {log}
tabix {output.norm_vcf} 2>> {log}
```

```bash
bcftools view \
    -e "QUAL<=100 | INFO/DP > 500 | INFO/DP < 10" \
    {input.vcf_norm} \
    > {output.vcf_filter} \
    2> {log}
```
