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

### FASTA

序列比对中用到的参考基因组通常是 FASTA 格式的。FASTA 格式非常简单，下面是一个 FASTA 文件的例子。

```
>chr02
ACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAACCCTAAACCCCTAACCCT
AAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAC
CCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCctaaaccctaaaccctaaacccta
aaccctaaaccctaaaaccctaaaccctaaaccctaaaccctaaaccctaaaccctaaac
cctaaaccctaaaccctaaccctaaaccctaaaccctaaaccctaaaccctaaaccctaa
>chr05
CCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTA
AACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAACCCTAAACCC
TAAAACCCTAAACCCTAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAA
CCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCT
AAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAACCCTAAACCCTAAACCCTAAACC
```

其中 `>` 代表一段序列的开始，其后跟着的 `chr02` 是这段序列的名字，而之后跟着的行就是序列的具体内容，换行符会被忽略，直到下一个 `>`。

### FASTQ

未经比对的原始测序数据通常是 FASTQ 格式的。下面是一个例子（`SRR15724032_2.fastq.gz` 文件的开头）。

```
@SRR15724032.1 1/2
CGCCCTCCTTAACCGTGTCGACAAGTCCTCCAGTAGATGAGCAGATGGGAACCACTCCATATCTCATCCCTTGCAATTGGATGAGACCACATGGCTCAAACCTACTAGGAACAATTATAAAATCAGCAGATCGGAAGAGCGTCGTGTAGG
+
FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SRR15724032.2 2/2
TTTCAGTTGACAGTTCTTCTAGCAGATACAGCTTCTTTGGCAAAAGGCCCAATCTGATGCAAAAGCTCAAGAAGTGGGGAAGGGGCAAGGATGACGGAAGCAGCTTAGCTTCACCGACACAGTCCTTCACTAGTGACTCCCCAAAGAGCG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF
```

一条 read 由 4 行组成。第一行以 `@` 开头，后面紧跟这条 read 的名字；第二行是这条 read 的序列；第三行通常以 `+` 作为分隔符，第四行则是代表测序质量的字符串，与第二行等长。

现在很多测序数据都是双端测序的，即测序时是从一个 DNA 片段的两端向中间测，所以一个测序数据由两个 FASTQ 文件组成（例如 `SRR15724032_1.fastq.gz` 和 `SRR15724032_2.fastq.gz`）。两个 FASTQ 里 reads 的顺序是匹配的，即文件 1 里的第 1 个 read 与 文件 2 的第 1 个 read 配对，文件 1 里的第 2 个 read 与 文件 2 的第 2 个 read 配对……

现在的序列比对软件对双端测序数据进行特殊处理，使得比对结果能够利用双端数据的特性变得更准。

### SAM/BAM

SAM 格式的全称是 Sequence Alignment/Map，它在 FASTQ 的基础上添加了每个 reads 的比对信息。而 BAM 就是 Binary 的 SAM。BAM 可保留 SAM 的所有信息，并且体积更小，所以一般都会将比对的结果转化为 BAM 保存。

## 质检和过滤

这里使用[简单的公共数据下载教程](/_drafts/simple-download-public-data-tutorial.md)中下载的 FASTQ 数据作为例子。

不同的平台，不同的公司可能产生不同的数据，所以对数据进行初步的质检和过滤是必要的。质检和过滤 FASTQ 的工具有很多，例如 [fastp](https://doi.org/10.1093/bioinformatics/bty560)，[trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) 等。这里以 fastp 为例。

fastp 可以用 apt 安装。

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

我们可以将 RUN ID 和样本名称写进一个关联数组中，然后用 for 循环进行批量处理。

```bash
declare -A srr_to_sample

srr_to_sample=(
    ["SRR15724033"]="DJ-1"
    ["SRR15724032"]="DJ-2"
    ["SRR15724031"]="DJ-3"
    ["SRR15724036"]="a-1"
    ["SRR15724035"]="a-2"
    ["SRR15724034"]="a-3"
    ["SRR15724028"]="A1-Input"
    ["SRR15724030"]="A1-1"
    ["SRR15724029"]="A1-2"
)

for srr in ${!srr_to_sample[@]}
do
    sample=${srr_to_sample[$srr]}

    if [[ -f ${sample}.html ]] ; then
        continue  # 跳过已经处理完的
    fi

    fastp \
        --thread 4 \
        --in1 ${srr}_1.fastq.gz \
        --in2 ${srr}_2.fastq.gz \
        --out1 ${sample}_R1.fastq.gz \
        --out2 ${sample}_R2.fastq.gz \
        --html ${sample}.html \
        --json ${sample}.json \
        --report_title ${sample}
done
```

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
    --quantMode GeneCounts \
    --outTmpDir /tmp/star
```

获得的 `DJ-3_Aligned.sortedByCoord.out.bam` 就是比对结果的 BAM 文件。

通过添加 `--quantMode GeneCounts` 选项，可以获得 `DJ-3_ReadsPerGene.out.tab` 文件，这个文件中有代表基因表达量的 counts 数据。其内部有 4 列：

1. 基因 ID
2. Map 到基因上的 no strand 的 reads 数（Ambiguous 的 reads 不会被计数）
3. Map 到基因上的 1st strand 的 reads 数
4. Map 到基因上的 2nd strand 的 reads 数

根据 STAR 文档中的介绍，这个文件的内容与默认设置的 `HTSeq-count` 的结果一致。根据[这个提问](https://www.biostars.org/p/218995/)，一般的 non-strand specific 的数据的差异表达分析应该用第二列 no strand 的 reads 数进行。

同样的，也可以写一个循环进行批量处理。

```bash
declare -a rna_samples

rna_samples=( "DJ-1" "DJ-2" "DJ-3" "a-1" "a-2" "a-3" )

for sample in ${rna_samples[@]}
do
    if [[ -f ${sample}_ReadsPerGene.out.tab ]] ; then
        continue
    fi
    
    ./STAR_2.7.11b/Linux_x86_64_static/STAR \
        --runThreadN 16 \
        --genomeDir os_star_index \
        --readFilesCommand zcat \
        --readFilesIn \
            ${sample}_R1.fastq.gz \
            ${sample}_R2.fastq.gz \
        --outFileNamePrefix ${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN 6 \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --quantMode GeneCounts \
        --outTmpDir /tmp/star
done
```

当完成了所有样本的比对，可将所有的 `ReadsPerGene.out.tab` 文件的内容合并为一个 counts 矩阵，用于进行后续的差异表达分析（参考：[简单的差异表达分析教程](/_drafts/simple-de-analysis-tutorial.md)。

### 基因组序列比对

可用于比对基因组测序数据的软件有很多，例如 [bwa](https://github.com/lh3/bwa)、[bowtie2](https://github.com/BenLangmead/bowtie2) 等。这里以 bwa 为例。

bwa 可以用 apt 安装，由于 bwa 只提供序列比对功能，所以还需要安装 samtools 来对比对结果进行排序。

```bash
sudo apt install bwa samtools
```

在用 bwa 进行序列比对之前，也需要先建一个索引。

```bash
bwa index IRGSP-1.0_genome.fasta
```

建完索引后，便可进行序列比对。

```bash
sample="A1-Input"
bwa_log=${sample}.bwa.log
out_bam=${sample}.bam

bwa mem -M -t 16 -R "@RG\\tID:${sample}\\tSM:${sample}" \
    IRGSP-1.0_genome.fasta ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz 2>> ${bwa_log} \
| samtools fixmate -u -m - - 2>> ${bwa_log} \
| samtools sort -u -@2 2>> ${bwa_log} \
| samtools markdup -u -@16 - - 2>> ${bwa_log} \
| samtools view -b -h -@16 -o ${out_bam} 2>> ${bwa_log}

samtools index ${out_bam} 2>> ${bwa_log}
```

同样的，也可以写一个循环进行批量处理。

```bash
declare -a samples

samples=( "A1-Input" "A1-1" "A1-2" )

for sample in ${samples[@]}
do
    bwa_log=${sample}.bwa.log
    out_bam=${sample}.bam

    if [[ -f $out_bam ]]; then
        continue
    fi

    bwa mem -M -t 16 -R "@RG\\tID:${sample}\\tSM:${sample}" \
        IRGSP-1.0_genome.fasta ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz 2>> ${bwa_log} \
    | samtools fixmate -u -m - - 2>> ${bwa_log} \
    | samtools sort -u -@2 2>> ${bwa_log} \
    | samtools markdup -u -@16 - - 2>> ${bwa_log} \
    | samtools view -b -h -@16 -o ${out_bam} 2>> ${bwa_log}

    samtools index ${out_bam} 2>> ${bwa_log}
done
```

