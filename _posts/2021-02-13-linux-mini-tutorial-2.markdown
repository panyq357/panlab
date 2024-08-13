---
layout: post
title:  "Linux 入门小教程（续）"
categories: 教程
---

在上一个小教程中，我们学习了目录的“增”、“删”、“改”、“查”和文件的“增”、“删”、“改”、“查”。

而最后我们通过 ``ls -al`` 命令看到了类似这样的一些东西。

````bash
$ ls -al
total 8
drwxr-xr-x 1 panyq panyq 4096 Feb  8 12:53 .
drwxr-xr-x 1 root  root  4096 Feb  8 12:46 ..
-rw------- 1 panyq panyq  128 Feb  8 12:53 .bash_history
-rw-r--r-- 1 panyq panyq  220 Feb  8 12:46 .bash_logout
-rw-r--r-- 1 panyq panyq 3771 Feb  8 12:46 .bashrc
drwxr-xr-x 1 panyq panyq 4096 Feb  8 12:46 .landscape
drwxr-xr-x 1 panyq panyq 4096 Feb  8 12:51 .local
-rw-r--r-- 1 panyq panyq    0 Feb  8 12:46 .motd_shown
-rw-r--r-- 1 panyq panyq  807 Feb  8 12:46 .profile
````

在输出结果中，除了第一行 ``total 8`` 以外，其余每一行都展现了当前目录下的一个文件或目录的详细信息。

以第 6 行为例，它可以被这样解读。

| 权限信息 | 连接数 | 所有者 | 所属组 | 大小（字节数） | 最后修改时间 | 名字 |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| ``-rw-r--r--`` | ``1`` | ``panyq`` | ``panyq`` | ``3771`` | ``Feb  8 12:46`` | ``.bashrc`` |

接下来简要介绍一下 Linux 的权限。

> 第一行的 ``total 8`` 和第二列的连接数有点复杂，这里就不讲了。😜

# 文件的权限

让我们 ``touch`` 一个空文件，并用 ``ls -l`` 查看它的详细信息。

````bash
$ touch abc
$ ls -l abc  # 只列出 abc 文件的详细信息。
-rw-r--r-- 1 panyq panyq 0 Feb 11 10:23 abc
````

让我们把输出中的权限信息 ``-rw-r--r--`` 放大一下。

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:--:|
| ``-`` | ``r`` | ``w`` | ``-`` | ``r`` | ``-`` | ``-`` | ``r`` | ``-`` | ``-`` |

第一个字符标识了该文件的类型（例如，普通文件是 ``-``，目录是 ``d``）。  
第 2、3、4 个字符标识了所有者的权限，依次为读（**r**ead）、写（**w**rite）和执行（e**x**ecute）。  
第 5、6、7 个字符标识了所属组的权限。（顺序同上）  
第 8、9、10 个字符标识了其他人的权限。

所以，把上面这串字符翻译成普通话就是：“这是一个普通文件，如果你是它的所有者，那么你可以读或写这个文件，但不能执行它。如果你不是它的所有者，那么这个文件对你来说是只读的。”

众所周知，Linux 是一个多用户的操作系统，所以为了保证每个用户的文件不会轻易地被别的用户修改，新建文件的默认权限就被设成了“只有所有者能读写，其他人只能读”了。

那么什么是执行权限呢？

在 Window 中，我们常常能碰到 ``.exe`` 结尾的可执行文件，它们通常是二进制的，可以被机器当作程序运行的文件。在 Linux 中也有这样的二进制地文件，但是为了安全起见，让机器运行这样的二进制文件也是需要权限的。

Linux 自带了许多可执行文件，比如位于 ``/usr/bin/`` 目录下的 ``ls`` 文件。

````bash
$ ls /usr/bin/ls -l
-rwxr-xr-x 1 root root 142144 Sep  5  2019 /usr/bin/ls
$
````

是的没错，这个 ``ls`` 文件就是 ``ls`` 命令的可执行文件，当我们键入 ``ls`` 命令时，系统就会去找这个 ``ls`` 文件，执行它，并打印它的输出结果。

通过查看它的权限信息，我们了解到：对于所有的用户，这个 ``ls`` 文件都是可执行的。

> 目录的权限与普通文件的权限有所不同，鸟哥曾把目录比喻为抽屉，把“读”比喻为“查看抽屉外面贴的标签”，把“写”比喻为“处置抽屉里的东西”，把“执行”比喻为“打开抽屉”。详情 -> [点这儿](https://linux.vbird.org/linux_basic/centos7/0210filepermission.php#filepermission_dir)

关于权限的内容先告一段落，现在让我们学学如何装软件。

## 软件安装

如同 Windows 10 有 Microsoft Store 一样，大部分 Linux 发行版都自带了像应用商店般的软件包管理程序。

我们所安装的发行版为 Ubuntu，它的“应用商店”叫 ``apt``（**A**dvanced **P**ackaging **T**ools）。

``apt`` 命令的用法十分简单，想要安装软件，只需键入 ``apt install <软件名字>`` 即可。

让我们来试试安装一个叫 ``tree`` 的软件。

````bash
$ apt install tree
E: Could not open lock file /var/lib/dpkg/lock-frontend - open (13: Permission denied)
E: Unable to acquire the dpkg frontend lock (/var/lib/dpkg/lock-frontend), are you root?
$
````

报错了，我们可以从错误信息中的 ``Permission denied`` 得知我们没有相应的权限。

就如同在 Windows 中安装软件常常需要管理员身份一样，在调用 ``apt`` 安装软件时，也需要管理员权限。

在错误信息的最后，系统提示我们 ``are you root?``，``root`` 即为 Linux 中管理员账户的名字。

管理员也被称为 Superuser，在 Linux 中，我们可以用 ``sudo`` 命令（**s**uper**u**ser **do**）临时获得管理员权限。

让我们再试一次，这次在命令的最前面加上 ``sudo``。（过程中会要求输入密码，就是那个安装 WSL 时设置的密码啦。😉）

````bash
$ sudo apt install tree
[sudo] password for panyq:  # 在此输入密码
Reading package lists... Done
Building dependency tree
Reading state information... Done
E: Unable to locate package tree
$
````

嗯，虽然这回我们有管理员权限了，但是又报错了。原来，我们想要安装的软件都放在互联网上的软件仓库里，要想用 ``apt`` 安装软件，需要先让 ``apt`` 知道远程的软件仓库里有哪些软件，它才能帮我们安装，所以我们需要更新一下软件源。

更新软件源的命令是 ``apt update``。（也需要管理员权限）

但如果就这样更新软件源，``apt`` 就会去连接 Ubuntu 官方的服务器，这个服务器在国外，从上面下载软件的速度会很慢。

所幸，国内有许多的软件源镜像站，比如清华大学的 [TUNA](https://mirrors.tuna.tsinghua.edu.cn/)。

要想让 ``apt`` 去连接 [TUNA](https://mirrors.tuna.tsinghua.edu.cn/) 而不是国外的服务器，我们需要进行一些设置。（以下内容来自 [TUNA 的教程](https://mirrors.tuna.tsinghua.edu.cn/help/ubuntu/)）

首先，我们需要以管理员身份，编辑软件源配置文件。

````bash
$ sudo mv /etc/apt/sources.list /etc/apt/sources.list.bak  # 备份下原始的配置文件。
$ sudo nano /etc/apt/sources.list  # 用 nano 编辑配置文件
````

然后把配置信息复制进去。（复制后用鼠标右键在黑框框内点一下即可粘贴）

````bash
# 默认注释了源码镜像以提高 apt update 速度，如有需要可自行取消注释
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal-updates main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal-updates main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal-backports main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal-backports main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal-security main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal-security main restricted universe multiverse

# 预发布软件源，不建议启用
# deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal-proposed main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ focal-proposed main restricted universe multiverse
````

保存退出后再用 ``sudo apt update`` 命令更新软件源。

````bash
$ sudo apt update
Get:1 https://mirrors.tuna.tsinghua.edu.cn/ubuntu focal InRelease [265 kB]
# 此处省略很多行……
Get:46 https://mirrors.tuna.tsinghua.edu.cn/ubuntu focal-security/multiverse amd64 c-n-f Metadata [284 B]
Fetched 19.9 MB in 5s (3863 kB/s)
Reading package lists... Done
Building dependency tree
Reading state information... Done
162 packages can be upgraded. Run 'apt list --upgradable' to see them.
$
````

这样，软件源就更新好了，现在用下面的命令安装 ``tree``。

````bash
sudo apt install tree
````

回车后安装会自动进行，当命令提示符再次出现时，安装就完成了。

那么 ``tree`` 是个做什么用的软件呢？请看下面这个例子。

````bash
$ mkdir outer
$ mkdir outer/inner
$ touch outer/inner/a_file
$ tree outer
outer
└── inner
    └── a_file

1 directory, 1 file
$
````

对，``tree`` 就是一个把文件和目录的结构以树形图的形式打印出来的小软件。

> 要卸载软件也很简单，只需 ``sudo apt remove <要卸载的软件名>`` 即可。

但是，不是所有软件都能在 ``apt`` 的软件源里找到的，有些软件需要我们手动从网上下载下来。所以接下来以 [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 和 [PLINK](http://www.cog-genomics.org/plink2/) 为例，介绍一下如何下载和解压网上的软件压缩包。

## 下载压缩包

从网上（**W**eb）获取（**Get**）压缩包，最简单的方法是使用 ``wget`` 命令，格式如下。

````bash
wget <下载链接>
````

让我们用它来试着下载 BLAST+ 和 PLINK。

![下载 BLAST+ 和 PLINK](/assets/images/linux-mini-tutorial/下载BLAST+和PLINK.gif)

下载完成后，用 ``ls`` 即可发现当前目录下多了两个文件。

````bash
$ ls
ncbi-blast-2.11.0+-x64-linux.tar.gz  plink_linux_x86_64_20201019.zip
````

一个以 ``.tar.gz`` 结尾，另一个以 ``.zip`` 结尾。这是两个不同格式的安装包，要用不同的命令进行解压。

首先是 ``ncbi-blast-2.11.0+-x64-linux.tar.gz``，要用 ``tar`` 命令进行解压。

````bash
$ tar -xf ncbi-blast-2.11.0+-x64-linux.tar.gz
$ # 没有进度条。提示符再次出现时，解压就完成了。
````

这里，``-x`` 表示我们要从一个压缩包中提取（e**x**tract）文件，``-f`` 选项则告诉 ``tar`` 要解压的文件的路径是 ``ncbi-blast-2.11.0+-x64-linux.tar.gz``。

> 注意，``-f`` 选项必须紧跟压缩包的路径，即不能写成 ``tar -fx ncbi-blast-2.11.0+-x64-linux.tar.gz‍``（之间隔了个 ``-x`` 选项）。

解压完成后，当前目录下就会多出一个 ``ncbi-blast-2.11.0+`` 目录，让我们用刚装的 ``tree`` 看看里头是啥。

````bash
$ tree ncbi-blast-2.11.0+
ncbi-blast-2.11.0+
├── ChangeLog
├── LICENSE
├── README
├── bin
│   ├── blast_formatter
│   ├── blastdb_aliastool
│   ├── blastdbcheck
│   ├── blastdbcmd
│   ├── blastn
│   ├── blastp
│   ├── blastx
│   ├── cleanup-blastdb-volumes.py
│   ├── convert2blastmask
│   ├── deltablast
│   ├── dustmasker
│   ├── get_species_taxids.sh
│   ├── legacy_blast.pl
│   ├── makeblastdb
│   ├── makembindex
│   ├── makeprofiledb
│   ├── psiblast
│   ├── rpsblast
│   ├── rpstblastn
│   ├── segmasker
│   ├── tblastn
│   ├── tblastx
│   ├── update_blastdb.pl
│   └── windowmasker
├── doc
│   └── README.txt
└── ncbi_package_info

2 directories, 29 files
$
````

嗯，可执行的二进制（**Bin**ary）文件都在其中的 ``bin`` 目录下了，我们我们只需输入这些软件的路径，即可调用这些软件。

````bash
$ ncbi-blast-2.11.0+/bin/blastn -h  # 用 -h 选项查看使用说明。
USAGE
  blastn [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
# 此处省略很多行。
    [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
    [-num_threads int_value] [-remote] [-version]

DESCRIPTION
   Nucleotide-Nucleotide BLAST 2.11.0+

Use '-help' to print detailed descriptions of command line arguments
$
````

BLAST+ 装好了，现在来看 PLINK。PLINK 的安装包是一个在 Windows 中常见的 ``.zip`` 格式的压缩包，在 Linux 中，并没有软件可以直接解压它，我们得安装 ``unzip`` 软件。

````bash
sudo apt install unzip
````

然后再用 ``unzip`` 解压压缩包。

````bash
panyq@CFDJH:~$ unzip plink_linux_x86_64_20201019.zip -d plink
Archive:  plink_linux_x86_64_20201019.zip
  inflating: plink/plink
  inflating: plink/LICENSE
  inflating: plink/toy.ped
  inflating: plink/toy.map
  inflating: plink/prettify
panyq@CFDJH:~$
````

这里，``-d`` 参数的作用是让 ``unzip`` 解压到指定的目录下（不然就都解压到当前目录下了😉）。

> 注：``tar`` 也可以指定解压的路径，比如可以写成这样：``tar -xf ncbi-blast-2.11.0+-x64-linux.tar.gz -C blast``（``blast`` 目录要事先存在）。

再用 ``tree`` 看看解压出的内容。

````bash
$ tree plink
plink
├── LICENSE
├── plink
├── prettify
├── toy.map
└── toy.ped

0 directories, 5 files
$
````

这个在 ``plink`` 目录下的 ``plink`` 文件就是 PLINK 的二进制可执行文件了，让我们来看看它的使用帮助。

````bash
$ plink/plink -h
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3

  plink <input flag(s)...> [command flag(s)...] [other flag(s)...]
  plink --help [flag name(s)...]

Commands include --make-bed, --recode, --flip-scan, --merge-list,
--write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,
--hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,
--blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,
--rel-cutoff, --cluster, --pca, --neighbour, --ibs-test, --regress-distance,
--model, --bd, --gxe, --logistic, --dosage, --lasso, --test-missing,
--make-perm-pheno, --tdt, --qfam, --annotate, --clump, --gene-report,
--meta-analysis, --epistasis, --fast-epistasis, and --score.

"plink --help | more" describes all functions (warning: long).
$
````

下载、解压软件压缩包的方法就介绍这么多。可以看出，比起用 ``apt`` 一键安装，手动下载解压要麻烦不少，所以接下来介绍一个安装生信常用软件的利器——conda。

## conda 和 bioconda

就如同安卓手机可以安装第三方应用商店一样，在 Linux 中也可以安装别的软件包管理程序，[conda](https://docs.conda.io/en/latest/) 就是这样的一个东西。

conda 最大的特点是可以创建虚拟环境。比如，我们想尝试两个不同的数据分析 protocol，但一个要求 python 的版本为 3.8，另一个要求 python 的版本为 3.2，我们就可以创建两个虚拟环境，将不同版本的 python 隔离开来，方便管理。

也正因如此，conda 受到了许多生信行业人员的青睐。为了能方便地安装各种生信软件，bioconda 应运而生。bioconda 是 conda 的一个 “channel”，囊括了许多的生信软件，这使得我们能通过 conda 一键安装各种常用的生信软件。

### 安装 conda

[Miniconda](https://docs.conda.io/en/latest/miniconda.html) 是一个包含着 conda 的安装包，它只包含 conda 以及 python 等 conda 的依赖软件，下面介绍下它的安装方法。

首先，在 [Miniconda 的官网](https://docs.conda.io/en/latest/miniconda.html)找到适用于 Linux 的安装包，复制它的下载链接，然后下载它。

````bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
````

安装包是一个 bash 脚本文件，它可以被这样打开。

````bash
$ bash Miniconda3-latest-Linux-x86_64.sh

Welcome to Miniconda3 py38_4.9.2

In order to continue the installation process, please review the license
agreement.
Please, press ENTER to continue
>>>
````

然后就进入安装界面了，根据提示按下回车，安装程序会调用一个叫 ``more`` 的文档查看器软件来显示用户协议，按住空格即可快速跳过。然后会询问是否接受用户协议。

````bash
Do you accept the license terms? [yes|no]
[no] >>>
````

输入 ``yes`` 后回车。然后会询问安装路径。

````bash
Miniconda3 will now be installed into this location:
/home/panyq/miniconda3

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/panyq/miniconda3] >>>
````

这里直接回车，按默认的来。之后还会提示是否初始化，``yes`` 即可😉。

安装完成后重启 WSL，命令提示符就变成了像这样的东西。

````bash
(base) panyq@CFDJH:~$
````

最左侧的 ``(base)`` 告诉我们，一个叫 ``base`` 的 conda 虚拟环境已被激活了，让我们再 ``ls`` 一下。

````bash
(base) panyq@CFDJH:~$ ls
Miniconda3-latest-Linux-x86_64.sh  miniconda3
````

这个 miniconda3 目录就是刚才安装的内容了。

### 配置 conda 的软件源。

就如同配置 ``apt`` 的软件源一样，我们也可以通过修改配置文件让 ``conda`` 到 TUNA 镜像站下载软件。

在 ``~`` 目录下新建一个叫 ``.condarc`` 的文件，然后把这些内容复制进去（内容修改自 [TUNA 的教程](https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/)）。

````yaml
channels:
  - bioconda
  - conda-forge
  - defaults
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
````

保存退出即可。

### 创建环境和安装软件

经过刚才的配置，我们就能用 ``conda`` 从 bioconda channel 中一键安装软件了，为了方便管理新安装的软件，让我们用下面的命令创建一个虚拟环境。

````bash
conda create --name temp
````

``--name`` 参数用来指定虚拟环境的名字。创建完成后，让我们根据提示，激活这个环境。

````bash
(base) panyq@CFDJH:~$ conda activate temp
(temp) panyq@CFDJH:~$
````

提示符最左侧发生了变化，说明 ``temp`` 虚拟环境被激活了。

> 这时，执行 ``conda deactivate`` 命令即可返回 base 环境，再来一次就会退出 base 环境了。

接下来让我们在 ``temp`` 虚拟环境中安装 BLAST+ 和 PLINK，用下面的命令即可。

```bash
conda install blast plink
```

回车后，conda 会开始计算需要安装哪些依赖，等它计算完成后，只需一个 ``y`` 确认，安装就自动进行了。

安装完成后，我们可以用 ``conda list`` 命令列出当前环境下安装的所有软件。

```bash
(temp) panyq@CFDJH:~$ conda list
# packages in environment at /home/panyq/miniconda3/envs/temp:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       1_gnu    conda-forge
blast                     2.9.0                h20b68b9_1    bioconda
boost                     1.68.0          py37h8619c78_1001    conda-forge
boost-cpp                 1.68.0            h11c811c_1000    conda-forge
bzip2                     1.0.8                h7f98852_4    conda-forge
ca-certificates           2020.12.5            ha878542_0    conda-forge
certifi                   2020.12.5        py37h89c1867_1    conda-forge
gmp                       6.2.1                h58526e2_0    conda-forge
gnutls                    3.6.13               h79a8f9a_0    conda-forge
icu                       58.2              hf484d3e_1000    conda-forge
ld_impl_linux-64          2.35.1               hea4e1c9_2    conda-forge
libblas                   3.9.0                8_openblas    conda-forge
libcblas                  3.9.0                8_openblas    conda-forge
libffi                    3.3                  h58526e2_2    conda-forge
libgcc-ng                 9.3.0               h2828fa1_18    conda-forge
libgfortran-ng            9.3.0               hff62375_18    conda-forge
libgfortran5              9.3.0               hff62375_18    conda-forge
libgomp                   9.3.0               h2828fa1_18    conda-forge
liblapack                 3.9.0                8_openblas    conda-forge
libopenblas               0.3.12          pthreads_h4812303_1    conda-forge
libstdcxx-ng              9.3.0               h6de172a_18    conda-forge
ncurses                   6.2                  h58526e2_4    conda-forge
nettle                    3.4.1             h1bed415_1002    conda-forge
numpy                     1.20.1           py37haa41c4c_0    conda-forge
openssl                   1.1.1i               h7f98852_0    conda-forge
pcre                      8.44                 he1b5a44_0    conda-forge
perl                      5.26.2            h36c2ea0_1008    conda-forge
perl-archive-tar          2.32                    pl526_0    bioconda
perl-carp                 1.38                    pl526_3    bioconda
perl-common-sense         3.74                    pl526_2    bioconda
perl-compress-raw-bzip2   2.087           pl526he1b5a44_0    bioconda
perl-compress-raw-zlib    2.087           pl526hc9558a2_0    bioconda
perl-exporter             5.72                    pl526_1    bioconda
perl-exporter-tiny        1.002001                pl526_0    bioconda
perl-extutils-makemaker   7.36                    pl526_1    bioconda
perl-io-compress          2.087           pl526he1b5a44_0    bioconda
perl-io-zlib              1.10                    pl526_2    bioconda
perl-json                 4.02                    pl526_0    bioconda
perl-json-xs              2.34            pl526h6bb024c_3    bioconda
perl-list-moreutils       0.428                   pl526_1    bioconda
perl-list-moreutils-xs    0.428                   pl526_0    bioconda
perl-pathtools            3.75            pl526h14c3975_1    bioconda
perl-scalar-list-utils    1.52            pl526h516909a_0    bioconda
perl-types-serialiser     1.0                     pl526_2    bioconda
perl-xsloader             0.24                    pl526_0    bioconda
pip                       21.0.1             pyhd8ed1ab_0    conda-forge
plink                     1.90b6.18            h516909a_0    bioconda
python                    3.7.9           hffdb5ce_0_cpython    conda-forge
python_abi                3.7                     1_cp37m    conda-forge
readline                  8.0                  he28a2e2_2    conda-forge
setuptools                49.6.0           py37h89c1867_3    conda-forge
sqlite                    3.34.0               h74cdb3f_0    conda-forge
tk                        8.6.10               h21135ba_1    conda-forge
wheel                     0.36.2             pyhd3deb0d_0    conda-forge
xz                        5.2.5                h516909a_1    conda-forge
zlib                      1.2.11            h516909a_1010    conda-forge
(temp) panyq@CFDJH:~$
```

这里头就有 ``blast`` 和 ``plink``，以及所有的依赖软件，是不是很方便？😉

### 卸载软件和删除环境

这里就直接贴命令啦😜。

````bash
conda remove plink  # 卸载当前环境中的 PLINK。
conda remove --name temp blast  # 卸载 temp 环境中的 BLAST+。
conda remove --name temp --all  # 删除 temp 环境，以及其中的所以软件。
````

其他常用指令可通过 [Cheat sheet](https://docs.conda.io/projects/conda/en/latest/_downloads/843d9e0198f2a193a3484886fa28163c/conda-cheatsheet.pdf) 快速查询。

这就是有关软件安装的全部内容啦。
