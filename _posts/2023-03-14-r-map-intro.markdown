---
layout: post
title:  "用 R 画地图"
date:   2023-03-14 21:30:00 +0800
categories: 教程
---

本文将简单的介绍一下如何在 R 中绘制世界行政区划，以及在地图上画点。

## 所需的 R 包

```r
library(tidyverse)
library(sf)
```

以上 R 包均可用 `install.packages()` 安装。

## 获取数据

首先，去地理所的数据中心下载[全球国家行政边界数据](https://www.resdc.cn/data.aspx?DATAID=205)（需要注册一下）。

解压后有以下几个文件。

```
世界国家
├── 世界国家.dbf
├── 世界国家.jpg
├── 世界国家.shp
└── 世界国家.shx
```

这一组文件就是是我们将要使用的 [shapefile](https://zh.wikipedia.org/zh-cn/Shapefile) 文件。

## 将数据读进 R

[sf](https://r-spatial.github.io/sf/) 包可将 shapefile 读进 R。

一个 shapefile 可能有很多个图层，可以用 `st_layers()` 函数列出 shapefile 的所有图层。

```r
st_layers("世界国家/世界国家.shp")
# Driver: ESRI Shapefile
# Available layers:
#   layer_name geometry_type features fields crs_name
# 1   世界国家       Polygon      247     10     <NA>
```

可以看到该 shapefile 只有一个图层。所以直接读入这个图层即可。

```r
countries <- st_read("世界国家/世界国家.shp", layer="世界国家")
```

## 绘制地图

我们可以用 ggplot 将世界地图画出来。

```r
library(ggplot2)
pdf("map1.pdf", width=12, height=7)
ggplot() +
    geom_sf(data=countries)
dev.off()
```

![map1](/assets/images/r-map-intro/map1.png)

可以看到图上的 x 轴和 y 轴对应着经度和纬度，而用于表示行政区划的多边形已经画到了对应的位置。

## 添加坐标点

如果要在地图上画点的话，只需再用 `geom_point()` 添加一个图层即可。

```r
points <- read_csv(
   "X,   Y,  Where
    116, 39, Beijing
    121, 31, Shanghai
    106, 29, Chongqing"
)

pdf("map2.pdf", width=12, height=7)
ggplot() +
    geom_sf(data=countries) +
    geom_point(data=points, aes(x=X, y=Y), color="blue", size=5)
dev.off()
```

![map2](/assets/images/r-map-intro/map2.png)

## 改变区块中的颜色

刚才读入的 `countries` 实际就是一个格式稍微特殊点的 data.frame。

```r
class(countries)
# [1] "sf"         "data.frame"
```

所以我们既可以向这个 data.frame 中赋值，也可以用其中的其他列来画图。

比如其中有一列 `POP`，应该对应的是国家的人口。

```r
pdf("map3.pdf", width=12, height=7)
ggplot() +
    geom_sf(data=countries, aes(fill=POP)) +
    scale_fill_viridis_b() +
    geom_point(data=points, aes(x=X, y=Y), color="blue", size=5)
dev.off()
```

![map2](/assets/images/r-map-intro/map3.png)
