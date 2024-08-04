---
layout: post
title:  "简单的公共数据下载教程"
---

现在许多文章都会将自己的原始数据提交到公共数据库中，本篇教程一这篇文章为例简单地介绍以下如何下载这些数据。

> Hua Wei, Hang Xu, Chen Su, Xiling Wang, Lei Wang, Rice CIRCADIAN CLOCK ASSOCIATED 1 transcriptionally regulates ABA signaling to confer multiple abiotic stress tolerance, Plant Physiology, Volume 190, Issue 2, October 2022, Pages 1057–1073, https://doi.org/10.1093/plphys/kiac196

在文章的 Materials and methods 中的的 Accession numbers 部分，可以找到测序数据在 NCBI 数据库中的编号：PRJNA754696。

这是一个 BioProject 编号，文章提交的 RNA-seq 和 DAP-seq 的数据都被包括在这一个 BioProject 中。

我们可以在 [EBI](https://www.ebi.ac.uk/ena/browser/home) 上搜索这个编号，找到对应的 BioProject（就是[这个](https://www.ebi.ac.uk/ena/browser/view/PRJNA754696)）。

然后找到 `Generated FASTQ files: FTP` 这一栏目，点击上面 `Download All` 按钮，便可下载下来一个 bash 脚本。

在安装好了 WSL 和 Windows Terminal 后（参考：[简单的 WSL 教程](/_drafts/simple-wsl-tutorial.md)），将这个脚本放到要装数据的文件夹里后，右键空白处，选择“在终端中打开”，便可在此处启动 WSL，然后键入如下命令：

```bash
bash <下载下来的脚本名>
```

便可用 wget 批量地下载 FASTQ 格式的数据了。

下载下来的文件名都是 SRR 开头的编号，要知道那个数据对应的是那个样本，可以在 EBI 的页面上点击 Show Column Selection，然后勾上 sample_alias，便可显示样本的名字了。（不同的数据提交者可能在不同的字段里写样品信息，可多勾选几个其他的字段看看）
