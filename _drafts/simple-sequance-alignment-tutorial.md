---
layout: post
title:  "简单的序列比对教程"
---

## 文件格式

在序列比对流程中，会遇到多种文件格式：

- 公司在测完序后会将数据返给我们，返的数据一般都是 gzip 压缩的 FASTQ 格式
- FASTQ 格式的数据经过序列比对软件处理后会变成 SAM/BAM 格式
- 而序列比对所需要的参考基因组数据通常是 FASTA 格式

了解这些文件格式，有助于更好地理解序列比对的过程中，数据都发生了什么变化。

### FASTQ

公司返的数据一般都是 FASTQ

## 质检和过滤

这里使用[简单的公共数据下载教程](/_drafts/simple-download-public-data-tutorial.md)中下载的 FASTQ 数据作为例子。

不同的平台，不同的公司可能产生不同的数据，所以对数据进行初步的质检和过滤是必要的。质检和过滤 FASTQ 的工具有很多，例如 [fastp](https://doi.org/10.1093/bioinformatics/bty560)，[trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) 等。这里以 fastp 为例。

使用下面的命令安装 fastp。

```bash
sudo apt install fastp
```

使用 fastp 过滤数据的例子如下。

```bash
fastp \
    --thread 4 \
    --in1 SRR15724031_1.fastq.gz \
    --in2 SRR15724031_2.fastq.gz \
    --out1 DJ-3_R1.fastq.gz \
    --out2 DJ-3_R2.fastq.gz \
    --html DJ-3.html \
    --json DJ-3.json \
    --report_title DJ-3
```

再过滤完成后，生成的 `DJ-1_R1.fastq.gz` 和 `DJ-3_R2.fastq.gz` 就是过滤后的 gzip 压缩的 FASTQ 数据。可在浏览器中打开 `DJ-3.html` 以查看质检报告。

## 序列比对

### 基因组序列比对

### 转录组序列比对

可用于比对转录组测序数据的软件有很多，例如 [STAR](https://doi.org/10.1093/bioinformatics/bts635)、[HISAT2](https://www.nature.com/articles/s41587-019-0201-4) 等。这里以 STAR 为例。

STAR 的编译好的二进制可执行文件可在 STAR 的 Github 仓库上找到：[点这里](https://github.com/alexdobin/STAR/releases)。将其下载到 WSL 中，解压后就可使用。

```bash
# 下载
wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip

# 解压
unzip STAR_2.7.11b.zip
```

解压出来的二进制文件（例如 `./STAR_2.7.11b/Linux_x86_64_static/STAR`）可放到 PATH 中，也可直接使用。

使用 STAR 进行转录组序列比对前需要建个索引，而建索引需要 FASTA 格式的参考基因组序列和 GTF 格式的参考基因组注释数据。

以水稻为例，这些数据可在 RAP-DB 上找到。将其下载下来并解压备用。

```bash
wget "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz"
wget "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_transcript_exon_2024-07-12.gtf.gz"
gzip -d IRGSP-1.0_genome.fasta.gz
gzip -d IRGSP-1.0_representative_transcript_exon_2024-07-12.gtf.gz
```

得到的 `IRGSP-1.0_genome.fasta` 和 `IRGSP-1.0_representative_transcript_exon_2024-07-12.gtf` 便可用来建索引了。

```bash
./STAR_2.7.11b/Linux_x86_64_static/STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir os_star_index \
    --genomeFastaFiles IRGSP-1.0_genome.fasta \
    --sjdbGTFfile IRGSP-1.0_representative_transcript_exon_2024-07-12.gtf \
    --genomeSAindexNbases 13
```

建完索引后，便可进行序列比对。

```bash
./STAR_2.7.11b/Linux_x86_64_static/STAR \
    --runThreadN 16 \
    --genomeDir os_star_index \
    --readFilesCommand zcat \
    --readFilesIn \
        DJ-3_R1.fastq.gz \
        DJ-3_R2.fastq.gz \
    --outFileNamePrefix DJ-3_ \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 6 \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --quantMode GeneCounts
```

获得的 `DJ-3_Aligned.sortedByCoord.out.bam` 就是比对结果的 BAM 文件，而 `DJ-3_ReadsPerGene.out.tab` 中则有每个基因上有多少个 reads 比对上的 counts 数据。其内部有 4 列：

1. 基因 ID
2. Map 到基因上的 no strand 的 reads 数（Ambiguous 的 reads 不会被计数）
3. Map 到基因上的 1st strand 的 reads 数
4. Map 到基因上的 2nd strand 的 reads 数

根据 STAR 文档中的介绍，这个文件的内容与默认设置的 `HTSeq-count` 的结果一致。根据[这个提问](https://www.biostars.org/p/218995/)，一般的 non-strand specific 的数据的差异表达分析应该用第二列 no strand 的 reads 数进行。

当完成了所有样本的比对，可将所有的 `ReadsPerGene.out.tab` 文件的内容合并为一个 counts 矩阵，用于进行后续的差异表达分析（参考：[简单的差异表达分析教程](/_drafts/simple-de-analysis-tutorial.md)。

