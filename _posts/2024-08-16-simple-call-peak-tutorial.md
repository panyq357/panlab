---
layout: post
title:  "简单的 Call Peak 教程"
---



## 使用 MACS Call Peak

[MACS](https://doi.org/10.1186/gb-2008-9-9-r137) 是最常见的 Call Peak 的软件, MACS 可用 apt 安装.

```bash
sudo apt install macs
```

安装完成后便可用 MACS call peak 了.

```bash
macs2 callpeak \
    --treatment A1-1.bam \
    --control A1-Input.bam \
    --format BAMPE \
    --gsize 370000000 \
    --name A1-1 \
    --outdir A1-1_peak \
    --qvalue 0.01
```

各个参数的意义如下：

- `--treatment`：要 call peak 的 BAM 文件
- `--control`：对应的对照或 input 的 BAM 文件，如果没有对照，也可以不要这项
- `--format`：输入文件的类型, 可选 `BAM`（单端测序）、`BAMPE`（双端测序）等
- `--gsize`：拟南芥大概 `120000000`, 水稻大概 `370000000`
- `--name`：输出文件的前缀
- `--outdir`：指定输出的目录
- `--qvalue`：peak 的 q-value 上限
- `--tempdir`：指定运行过程中临时文件存放目录

通过键入上面的命令，macs2 将在指定的输出目录下生成五个前缀相同，后缀不同的文件。

- `{prefix}_peaks.xls`：一个包含 peak 信息及分析时的各种参数的表格，其中的坐标与 0-base 的 BED 文件不同，是 1-base 的
- `{prefix}_summits.bed`：包含 peak 最高点位置信息的 BED 格式文件
- `{prefix}_peaks.narrowPeak`：包含 peak 位置和统计信息的 BED 6+4 格式文件
