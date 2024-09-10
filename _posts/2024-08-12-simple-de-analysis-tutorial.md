---
layout: post
title:  "简单的差异表达分析教程"
categories: 教程
---

## 获得 counts 矩阵

在[简单的序列比对教程]({% post_url 2024-08-12-simple-sequance-alignment-tutorial %})中，我们用 STAR 对转录组数据进行了序列比对，并通过添加 `--quantMode GeneCounts` 选项，得到了每个基因的 counts。

接着，我们便可以将数据读进 R，将数据合并成为一个 counts 矩阵。

```r
# 获得所有 STAR 生成的文件的文件名
tabs <- list.files(".", pattern="*_ReadsPerGene.out.tab")

# 从文件名中提取样本名字
names(tabs) <- sub("(.*)_ReadsPerGene.out.tab", "\\1", tabs)

# 使用 lapply 处理每一个文件，最终获得一个列表，列表中的每个元素是一个样本的 counts
counts_list <- tabs |> 
    lapply(function(x) {
        df <- read.table(x, header=F)
        df <- df[-(1:4),]
        counts <- setNames(df[[2]], df[[1]])
        return(counts)
    })

# 生成一个空 data.frame，将行名设为所有样本中测到的基因的 ID
counts <- data.frame(row.names = lapply(counts_list, names) |> unlist() |> unique() |> sort())

# 合并所有样本的 counts
for (name in names(counts_list)) {
    counts[[name]] <- counts_list[[name]][row.names(counts)]
}

# 输出结果
write.csv(counts, "counts.csv")
```

这样便可以获得一个 counts 矩阵 csv 文件，其形状如下所示。

```csv
"","a-1","a-2","a-3","DJ-1","DJ-2","DJ-3"
"Os01g0100100",673,703,604,598,834,739
"Os01g0100200",0,0,0,0,0,1
"Os01g0100300",0,0,0,0,0,0
"Os01g0100400",62,42,53,79,84,82
"Os01g0100466",0,0,0,0,0,0
"Os01g0100500",1668,1405,1259,1675,1673,1447
"Os01g0100600",472,386,377,381,381,357
"Os01g0100650",0,0,0,0,0,0
"Os01g0100700",474,550,459,556,631,730
```

每一行是一个基因，每一列是一个样本，矩阵中的数字就是这个基因在这个样本中测到的次数。

## 差异表达分析

在获得了 counts 矩阵后，便可进行差异表达分析。差异表达分析的工具有很多，例如 [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)、[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) 等，这里以 DESeq2 为例。

```r
library(DESeq2)

# 读入 counts 矩阵
cts <- read.csv("counts.csv", check.names=F, row.names=1)

# 生成一个装有样本分组信息的 data.frame
coldata <- data.frame(
    row.names = names(cts),  # 行名要和矩阵的列名一样
    Group = factor(sub("([^-]+)-.*", "\\1", names(cts)), levels=c("DJ", "a"))  # 提取样本 ID 中横杠 "-" 前的部分作为样本的分组名
)

# 对 counts 矩阵进行初步过滤，去掉在所有样本中被测到的总次数不大于 10 的基因
cts <- cts[rowSums(cts) >= 10,]

# 生成一个 DESeq 对象，并告知用 coldata 中的 Group 作为差异表达分组的依据
dds <- DESeqDataSetFromMatrix(cts, coldata, ~ Group)

# 运行 DESeq
dds <- DESeq(dds)

# 表达量不高的基因会容易有较高的 Log2FoldChange（例如从 1 > 2 相比于从 11 > 12 翻的倍数要更多）
# 使用 lfcShrink 可对表达量低的基因的 Log2FoldChange 进行缩减，使得获得的数值更有意义
# 处理后的 Log2FoldChange 也更适合下游的 GSEA 分析
# 其中 colnames(coef(dds))[2] 的实际值是 Group_a_vs_DJ
res <- lfcShrink(dds, colnames(coef(dds))[2], type="apeglm")

# 将结果转换成 data.frame
res_df <- data.frame(res, check.names=F)

# 将 DESeq 标准化后的 counts 矩阵贴在表格的后面
nc <- counts(dds, normalized=TRUE)
res_df <- cbind(res_df, nc[row.names(res_df),])

# 按 padj 列进行排序
res_df <- res_df[order(res_df$padj, decreasing=FALSE),]

# 输出结果
saveRDS(res, "deseq2_res.rds")
write.csv(res_df, "deseq2_res.csv")
```

这样获得的 `deseq2_res.csv` 是一个表格，内有每一个基因的检验结果以及标准化后的 counts 矩阵。

## 贴注释

`deseq2_res.csv` 文件中只有基因的基因号，这样看起来不方便。我们可以从 RAP-DB 下载水稻基因的注释数据，然后把基因注释也贴在表格的后面，这样看起来就方便了。

首先使用 wget 下载注释数据。

```bash
wget "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_annotation_2024-07-12.tsv.gz"
```

然后用 R 语言脚本将基因的注释添加到差异表达分析结果的后面。

```r
# 读入注释数据
rap_anno <- readr::read_tsv("IRGSP-1.0_representative_annotation_2024-07-12.tsv.gz")

# 读入之前的结果
deseq2_res <- read.csv("deseq2_res.csv", row.names=1, check.names=F)

# 按基因的 ID，将注释数据添加到表的后面
deseq2_res[c("Name", "Description")] <- rap_anno[
    match(row.names(deseq2_res), rap_anno$Locus_ID),
    c("CGSNL Gene Name", "Description")
]

# 输出结果
write.csv(deseq2_res, "deseq2_res_anno.csv")
```

## 后续

在获得了差异表达分析结果后，我们可以：

- 用差异表达分析结果中的 Log2FoldChange 和 padj，画一个火山图（参考[简单的 ggplot2 教程]({% post_url 2024-08-16-simple-ggplot2-tutorial %})）
- 提取有差异的基因的 ID，进行富集分析（参考[简单的富集分析教程]({% post_url 2024-08-12-simple-enrichment-analysis-tutorial %})）
