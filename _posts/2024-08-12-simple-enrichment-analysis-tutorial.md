---
layout: post
title:  "简单的富集分析教程"
categories: 教程
---

## 富集分析的原理

假设有两个罐子，一个罐子里有 8 个白球和 12 个黑球。

从这个罐子里随机无放回地抽 10 个球。如果是 4 白 6 黑，我们会觉得很正常。但如果是 7 白 3 黑，我们会觉得这是个极小概率事件，进而怀疑这个罐子是不是有问题？或者抽样方法有问题？

我们可以用 R 语言计算抽 10 个球有 7 个或更多白球的概率。

```r
1 - phyper(q=7, m=8, n=12, k=10)
# [1] 0.0003572279
```

现在把球换成基因，假设我们通过[差异表达分析]({% post_url 2024-08-12-simple-de-analysis-tutorial %})获得了 1000 个差异表达的基因，其中有 100 个基因调控株高，而已知水稻中有 30000 个基因，其中有 300 个基因调控株高。我们想知道差异表达的基因中，调控株高的是不是明显比随机抽样要更多，我们也可以计算随机抽样的情况下，抽 1000 个基因中有 100 个是调控株高的基因的概率。

```r
1 - phyper(q=100, m=300, n=30000-300, k=1000)
# [1] 0
```

于是我们知道，随机从 30000 个基因里抽 1000 个基因，其中有 100 个或更多基因是调控株高的基因的概率为 0。进而我们就可以说，这 1000 个差异表达的基因中调控株高的基因显著地要多。

现在我们不只想知道调控株高的基因是否显著地多，还想知道其他功能的基因是否显著地多。于是我们便可以这样做：

1. 统计这 1000 个基因中每种功能都有几个基因（株高有几个，分蘖有几个……）
2. 统计全基因组范围的所有的 30000 个基因中，每种功能都有几个基因
3. 对每种功能，都用上面的方法计算一遍概率
4. 将结果汇总在一个表里，按概率从低到高排，便知道这 1000 个基因中，调控哪些功能的基因的数量显著地要多了

这便是富集分析所要做的事了。

## 基因的注释数据

以水稻为例，在 [Oryzabase](https://shigen.nig.ac.jp/rice/oryzabase/) 中搜索一个基因，我们可以搜到这个基因的 Gene Ontology，Trait Ontology 和 Plant Ontology 注释。

这些 Ontology 注释就是描述这个基因的功能的，一种标准化的标签。

所有基因的所有注释数据可在 [Oryzabase 的下载区](https://shigen.nig.ac.jp/rice/oryzabase/download/gene)下载。

而我也整理了一份可直接用于富集分析的数据放在了一个 [Github 仓库](https://github.com/panyq357/onto-scripts)里（[点此下载](https://github.com/panyq357/onto-scripts/releases/download/2024-08-08/oryzabase-ontologies.xlsx)）。

## 用 clusterProfiler 进行富集分析

[clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html) 是一个专门用来进行富集分析、GSEA 等分析的 R 包，可用下面的命令在 R 语言中安装。

```r
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
```

然后，就可以用 clusterProfiler 提供的函数和整理好的注释数据，对未知功能的基因集进行注释了。

```r
# 读入注释数据
# Excel 文件里有不同的 sheet，对应不同的 Ontology 和不同的基因 ID
# 比如将 sheet 换成 "MSU_TO"，便可用 MSU 的 ID 及进行 Trait Ontology 的富集分析了
onto = readxl::read_excel("oryzabase-ontologies.xlsx", sheet="RAP_GO")

# 这里给要富集分析的基因集的基因号
gene = c("Os01g0118100", "Os01g0549700", "Os02g0710800", "Os03g0108600", "Os03g0158200", "Os03g0746500")

universe = NULL  # NULL 意为将所有基因作为背景，更常见的做法是将 RNA-seq 所有能测到的基因作为背景

# 富集分析
enrich_res = clusterProfiler::enricher(
    gene=gene,
    universe=universe,
    TERM2GENE=onto[c("OntoID", "GeneID")],
    TERM2NAME=onto[c("OntoID", "Description")]
)

# 导出结果
write.csv(as.data.frame(enrich_res), "enrich_res.csv")
```

输出的 `enrich_res.csv` 文件里就写了每种 Ontology 的富集结果了。

clusterProfiler 也提供了一些方便的绘图函数，比如绘制气泡图的 `dotplot()`。

```r
svg("demo_dotplot.svg")
clusterProfiler::dotplot(enrich_res)
dev.off()
```
