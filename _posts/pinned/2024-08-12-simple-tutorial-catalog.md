---
layout: post
title:  "简单的测序数据分析教程的目录"
categories: 教程
pinned: true
---

在三代测序逐渐普及，二代测序费用白菜化的现在，越来越多的研究工作中用到了 RNA-seq、ChIP-seq、BSA-seq 等技术。为了方便新手入门，这里总结了一些常用的数据分析教程。

首先，要进行数据分析，需要了解基本的运行软件的环境。虽然现在可用的软件包越来越丰富，但是大部分测序分析流程还是无法完全避免在 Linux 命令行下进行操作。微软从 Windows 10 开始，就为 Windows 内置了 Windows Subsystem for Linux（WSL）的功能，使在 Windows 下也能方便的运行 Linux 的软件。关于 WSL 的安装和使用可以参考[简单的 WSL 教程]({% post_url 2024-08-12-simple-wsl-tutorial %})。

然后，要进行测序数据的分析，需要有数据，关于数据的获取可以参考[简单的公共数据下载教程]({% post_url 2024-08-12-simple-public-data-download-tutorial %})。

在获得了原始测序数据后，几乎所有对测序数据进行分析的流程的起点都是序列比对，关于序列比对可以参考[简单的序列比对教程]({% post_url 2024-08-12-simple-sequance-alignment-tutorial %})

对于转录组数据分析，在序列比对完之后，会：

- 用 DESeq2 等软件进行差异表达分析（参考[简单的差异表达分析教程]({% post_url 2024-08-12-simple-de-analysis-tutorial %})）
- 用 clusterProfiler 等软件按进行差异表达基因的功能富集（参考[简单的富集分析教程]({% post_url 2024-08-12-simple-enrichment-analysis-tutorial %})）

对于 ChIP-seq 数据分析，在序列比对完之后，会：

- 用 MACS 等软件 call peak（参考[简单的 Call Peak 教程]({% post_url 2024-08-16-simple-call-peak-tutorial %})）
- 用 deeptools 等软件进行 peak 的可视化（参考[简单的 deeptools 教程]({% post_url 2024-08-16-simple-deeptools-tutorial %})）
- 用 MEME-ChIP 等软件分析结合的 motif（参考[简单的 MEME-ChIP 教程]({% post_url 2024-08-16-simple-meme-chip-tutorial %})）

对于 BSA-seq 数据分析，在序列比对完之后，会：

- 用 bcftools 等软件 call SNP（参考[简单的 Call SNP 教程]({% post_url 2024-08-16-simple-call-snp-tutorial %})）
- 用 VEP 等软件注释 SNP（参考[简单的 VEP 教程]({% post_url 2024-08-16-simple-vep-tutorial %})）
- 用 PyBSASeq 等软件计算连锁统计量，并绘图（参考[简单的 PyBSASeq 教程]({% post_url 2024-08-16-simple-pybsaseq-tutorial %})）
