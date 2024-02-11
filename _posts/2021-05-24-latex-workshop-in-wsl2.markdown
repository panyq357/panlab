---
layout: post
title:  "在 WSL 2 中使用 LaTeX Workshop"
date:   2021-05-25 20:29:00 +0800
categories: 笔记
---

本文记录了在 WSL 2 中安装 TeX Live 以及配置 VS Code 的 LaTeX Workshop 插件的过程。

## 1. 安装 TeX Live

首先下载 TeX Live 的安装镜像文件，这里我选择从[中科大的镜像](http://mirrors.ustc.edu.cn/CTAN/systems/texlive/Images/)下载。

```bash
wget 'http://mirrors.ustc.edu.cn/CTAN/systems/texlive/Images/texlive.iso'
```

下载完成后，挂载这个镜像文件。

```bash
sudo mkdir /mnt/texlive
sudo mount texlive.iso /mnt/texlive/
```

然后执行其中的 ``install-tl`` 脚本。

```bash
cd /mnt/texlive/
sudo ./install-tl
```

```txt
 <O> options:
   [ ] use letter size instead of A4 by default
   [X] allow execution of restricted list of programs via \write18
   [X] create all format files
   [X] install macro/font doc tree
   [X] install macro/font source tree
   [ ] create symlinks to standard directories
   [X] after install, set CTAN as source for package updates

 <V> set up for portable installation

Actions:
 <I> start installation to hard disk
 <P> save installation profile to 'texlive.profile' and exit
 <Q> quit

Enter command:
```

启动脚本后，终端中会呈现出一个安装界面，这里需要先键入 ``O``，修改一个选项。

```txt
===============================================================================
Options customization:

 <P> use letter size instead of A4 by default: [ ]
 <E> execution of restricted list of programs: [X]
 <F> create all format files:                  [X]
 <D> install font/macro doc tree:              [X]
 <S> install font/macro source tree:           [X]
 <L> create symlinks in standard directories:  [ ]
            binaries to:
            manpages to:
                info to:
 <Y> after install, set CTAN as source for package updates: [X]

Actions: (disk space required: 6963 MB)
 <R> return to main menu
 <Q> quit

Enter command:
```

键入 ``L``，然后**连敲三次回车**（使用默认的路径），这样就能让安装脚本自动创建软连接，省去安装后手动配置环境变量。

然后再键入 ``R`` 返回主安装界面，再键入 ``I`` 选项即可开始安装，安装时间大概需要十几分钟。

## 2. 让 LaTeX Workshop 能在 WSL2 中工作

其实只要安装完了 TeX Live 并配置好了环境变量，LaTeX Workshop 应该是开箱即用的，
但由于我用的是 WSL 2，所以根据 LaTeX Workshop 的 [FAQ](https://github.com/James-Yu/LaTeX-Workshop/wiki/FAQ#using-latex-workshop-with-wsl)，
还需要在 VS Code 的 ``settings.json`` 里加上一条。

```json
"latex-workshop.latex.watch.usePolling ": true
```

然后 LaTeX Workshop 就能正常工作了。

## 3. 让 LaTeX Workshop 能编译中文文档

出于某种原因，``ctex`` 宏包和 ``pdflatex`` 并不兼容，但 LaTeX Workshop 默认是用 ``pdflatex`` 的，要想用 ``ctex`` 宏包书写中文的内容，需要修改下 LaTeX Workshop 的编译设置。

将下列内容放入 VS Code 的 ``settings.json`` 中，即可用一个简单的，可生成中文 PDF 的编译方法覆盖掉原来那些复杂的方法。

```json
"latex-workshop.latex.recipes":[
    {
        "name": "simple", // 一个简单的编译方法。
        "tools": [
            "latexmk_xelatex"  // 这个编译方法会去找下面叫 latexmk_xelatex 的 tool。
        ]
    }
],
"latex-workshop.latex.tools":[
    {
        "name": "latexmk_xelatex",  // 被上面的 recipe 调用。
        "command": "latexmk",  // 用 latexmk 进行编译。
        "args": [
            "-xelatex",  // 改用 xelatex 生成 PDF。
        ],
    }
]
```

上面这些内容的意思是，当按下 LaTeX Workshop 的编译按钮时，它就会在终端里替你输入这一条命令

```bash
latexmk -xelatex <源文件>
```

这样就能生成中文的 PDF 了。现在试着编译这一段吧。

```latex
\documentclass{article}
\usepackage{ctex}
\begin{document}
    你好，世界！
\end{document}
```

## 参考

1. 知乎上的 TexLive 2020 安装指南: <https://zhuanlan.zhihu.com/p/136931926>
2. TeX Live - Quick install: <https://www.tug.org/texlive/quickinstall.html>
3. LaTeX Workshop 的 FAQ: <https://github.com/James-Yu/LaTeX-Workshop/wiki/FAQ#using-latex-workshop-with-wsl>