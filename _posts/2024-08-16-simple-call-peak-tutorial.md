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
    --treatment ChIP.bam \
    --control Control.bam \
    --format BAMPE \
    --gsize 370000000 \
    --name my_sample \
    --outdir my_outdir \
    --qvalue 0.01
```

各个参数的意义如下：

- `--treatment`：要 call peak 的 BAM 文件
- `--control`：对应的 Input 的 BAM 文件
- `--format`：输入文件的类型, 可选 `BAM`（单端测序）、`BAMPE`（双端测序）等
- `--gsize`：拟南芥大概 `120000000`, 水稻大概 `370000000`
- `--name`：输出文件的前缀
- `--outdir`：指定输出的目录
- `--qvalue`：peak 的 q-value 上限
- `--tempdir`：指定运行过程中临时文件存放目录

