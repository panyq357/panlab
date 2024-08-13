---
layout: post
title:  "简单的差异表达分析教程"
categories: 教程
---

## 获得 counts 矩阵

在[简单的序列比对教程](/_drafts/simple-sequance-alignment-tutorial.md)中，我们用 STAR 对转录组数据进行了序列比对，并通过添加 `--quantMode GeneCounts` 选项，得到了每个基因的 counts。

接着，我们便可以将数据读进 R，将数据合并成为一个 counts 矩阵。

```r
tabs <- list.files(".", pattern="*_ReadsPerGene.out.tab")

names(tabs) <- sub("(.*)_ReadsPerGene.out.tab", "\\1", tabs)

counts_list <- tabs |> 
    lapply(function(x) {
        df <- read.table(x, header=F)
        df <- df[-(1:4),]
        counts <- setNames(df[[2]], df[[1]])
        return(counts)
    })

counts <- data.frame(row.names = lapply(counts_list, names) |> unlist() |> unique() |> sort())

for (name in names(counts_list)) {
    counts[[name]] <- counts_list[[name]][row.names(counts)]
}

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

## 差异表达分析

差异表达分析的工具有很多，例如 [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)、[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) 等，这里以 DESeq2 为例。

```r
library(DESeq2)

cts <- read.csv("counts.csv", check.names=F, row.names=1)

coldata <- data.frame(
    row.names = names(cts),
    Group = factor(sub("([^-]+)-.*", "\\1", names(cts)), levels=c("DJ", "a"))
)

dds <- DESeqDataSetFromMatrix(cts, coldata, ~ Group)

# Pre-filtering.
dds <- dds[rowSums(counts(dds)) >= 10,]

# DESeq.
dds <- DESeq(dds)

# Use lfcShrink to get a better Log2FoldChange for downstream GSEA.
# colnames(coef(dds))[2] looks like: "Group_mutent_vs_wildtype".
res <- lfcShrink(dds, colnames(coef(dds))[2], type="apeglm")

res_df <- data.frame(res, check.names=F)

# Combine normalized counts.
nc <- counts(dds, normalized=TRUE)
res_df <- cbind(res_df, nc[row.names(res_df),])

# Sort by padj
res_df <- res_df[order(res_df$padj, decreasing=FALSE),]

## Output DESeq results.
saveRDS(res, "deseq2_res.rds")
write.csv(res_df, "deseq2_res.csv")
```

获得的 `deseq2_res.csv` 是一个表格，内有每一个基因的检验结果以及 normalized counts。

## 差异基因注释

### 贴注释

我们可以从 RAP-DB 下载水稻基因的注释。

```bash
wget "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_annotation_2024-07-12.tsv.gz"
```

然后用 R 语言脚本将基因的注释添加到差异表达分析结果的后面。

```r
rap_anno <- readr::read_tsv("IRGSP-1.0_representative_annotation_2024-07-12.tsv.gz")

deseq2_res <- read.csv("deseq2_res.csv", row.names=1, check.names=F)

deseq2_res[c("Name", "Description")] <- rap_anno[
    match(row.names(deseq2_res), rap_anno$Locus_ID),
    c("CGSNL Gene Name", "Description")
]

write.csv(deseq2_res, "deseq2_res_anno.csv")
```
