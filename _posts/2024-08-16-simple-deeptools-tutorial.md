---
layout: post
title:  "简单的 deeptools 教程"
---

[deeptools](https://deeptools.readthedocs.io/en/stable/) 是一个常用的对序列比对结果的 BAM 进行进一步处理和可视化的软件。

在进行序列比对和 call peak 后，

这里我们以[简单的序列比对教程]({% post_url 2024-08-12-simple-sequance-alignment-tutorial %})中基因组序列比对的结果，以及[简单的 Call Peak 教程]({% post_url 2024-08-16-simple-call-peak-tutorial %})中 call 的 peak 为例，介绍一下 deeptools 的使用方法。

## 0. 安装 deeptools

安装 deeptools 的方法可参考 [deeptools 文档的 Installation 部分](https://deeptools.readthedocs.io/en/stable/content/installation.html)。除此之外，deeptools 也可用 `apt` 进行安装。

```bash
sudo apt install python3-deeptools
```

## 1. 将 BAM 转换为 BigWig 格式

使用 deeptools 进行可视化之前，需要将 BAM 文件转换为 BigWig 文件，deeptools 中负责这项操作的是 `bamCoverage` 命令。

一个抄自 [HBC Traning](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html)的使用例子如下。

```bash
bamCoverage \
    --bam A1-Input.bam \
    --outFileName A1-Input.bw \
    --binSize 20 \
    --normalizeUsing BPM \
    --smoothLength 60 \
    --extendReads \
    --centerReads \
    --numberOfProcessors 6
```

其中：

- `--normalizeUsing BPM` 使 reads 数被标准化，类似于 RNA-seq 标准化中的 TPM
- `--smoothLength 60` 使每一个 bin 的值是由其周围的一个范围取平均值得到的，这使获得的数值更平滑
- `--extendReads` 和 `--centerReads` 将双端的 reads 延伸为两端对应的 fragment，然后将位置移到 fragment 的中间，这使富集的区域更集中

获得的 BigWig 文件（`A1-Input.bw`）可以直接拖到 IGV 中进行可视化。

## 2. 

