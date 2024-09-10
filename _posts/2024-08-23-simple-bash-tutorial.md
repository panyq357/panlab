---
layout: post
title:  "简单的 Bash 教程"
categories: 教程
---

这篇教程介绍了一些常用的 Bash 终端使用方法。

## 0. 内核、外壳和终端

首先需要解释清楚一些基本的概念。

Linux 一词既可指 Linux 操作系统，也可指该操作系统的**内核**（kernel）。内核是一组软件，负责计算机硬件（如 CPU、内存、硬盘等）的调度。

内核的结构很复杂，普通的计算机用户难以直接操作内核，所以在 Linux 系统的内核之上，还有一层外壳，即 **shell**，它是介于内核与各种应用程序的中间层。

常见的 shell 有 Bash、Fish、Zsh 等。在 shell 中可以运行命令和脚本，对于操作系统中一些自启动的基本的程序，它们会启动一些没有用户界面的 shell，来运行这些程序。

而我们人类用户想要运行命令时，我们通常是通过**终端**这个用户界面，将命令输入到 shell 中进行运行的。

在通过[简单的 WSL 教程]({% post_url 2024-08-12-simple-wsl-tutorial %})，安装并启动 WSL 后，就会进入一个终端界面，这是一个 Bash 的 shell 的终端。

在这个终端中，最底下一行应该是这样的。

```bash
<用户名>@<机器名>:~$ 
```

这是一个 prompt，`<用户名>@<机器名>` 告诉了你：你是谁，以及这是哪台机器；`:` 与 `$` 符号中间的 `~` 告诉了你：当前所在的目录是家目录；`$` 后面闪烁的光标则告诉了你：可以在此输入**命令**。

## 1. 命令

普通的命令由单词和空格组成。以下是一些例子（`#` 之后的内容为注释，它们不会被 shell 解释执行）。

```bash
# 显示当前目录的绝对路径
pwd

# 详细显示显示当前目录中的所有内容，并将文件大小显示为人类可读的样式
ls -alh

# 新建一个名为 docs 的目录
mkdir docs

# 在 docs 目录下生成 3 个文件：file1、file2 和 file3
touch file{1..3}

# 复制 file1 文件为 file4
cp file1 file4

# 将文件名以 file 开头的文件移动至 docs 目录内
mv file* docs

# 将 docs 目录及其内容压缩为 archive.tar.gz
tar -cf archive.tar.gz docs

# 强制删除 archive.tar.gz 文件、docs 目录以及 docs 目录里的所有内容
rm -rf archive.tar.gz docs
```

在上述例子中，像 `pwd`、`ls`、`tar` 等就是命令的名字，`-alh`、`-rf` 等以横杠 `-` 开头的是命令的选项，之后的 `file1`、`archive.tar.gz` 等则是命令的参数。命令、选项和参数之间通常至少要有 1 个空格分隔。

以 `ls -alh` 命令为例，这是一个 `ls` 命令，以及 3 个选项：`-a`（显示 `.` 开头的隐藏文件和目录）、`-l`（显示文件详细信息）和 `-h`（将文件大小显示为人类可读的带单位的格式）。3 个选项可以分开写（例如 `ls -a -l -h`），也可以缩写在一起（即 `ls -alh`）。

而在 `tar -cf archive.tar.gz docs` 中，用了 2 个选项：`-c`（进行压缩操作）和 `-f`（指定压缩包文件名）。`-f` 选项后必须立即跟上对应的参数（即压缩包的名字 `archive.tar.gz`）。这条命令也可以写成 `tar -c docs -f archive.tar.gz`（注意 `-f` 选项和压缩包名字 `archive.tar.gz` 的相对位置）。

可使用 `man` 命令查看某个命令的手册（manual），例如 `man tar`。

## 2. 路径

在 Linux 中，一切皆为文件。Linux 的目录是一个以根目录 `/` 为起始的树状结构。

```bash
# 列出根目录的内容
ls /
```

Linux 系统中有着数不清的目录，而我们正处于其中一个目录之中，这个目录也被称作**工作目录**。我们可以用 `pwd` 打印工作目录的路径，也可以用 `cd` 切换工作目录。

```bash
# 切换工作目录到根目录下的 mnt 目录中
cd /mnt
```

如果你使用的是 WSL，现在再键入 `ls`，你会发现 `/mnt` 有着 `c`、`d` 等目录，这些目录就是 Windows 中的磁盘在 WSL 中对应的目录。

对于每个 Linux 中的用户，都有一个“家目录”，除此之外，对于“刚才那个目录”、“上一层目录”、“当前工作目录”等常用的目录，也有对应的简写。

```bash
# 切换工作目录至家目录
cd ~

# 切换工作目录至刚才那个目录
cd -

# 列出工作目录上一层的内容
ls ..

# 列出当前工作目录的内容
ls .
```

Linux 中存在有两种路径：绝对路径和相对路径。绝对路径是相对于系统根目录的路径，而相对路径是相对于当前工作目录的路径。

绝对路径一般以 `/` 开头（例如 `/home/user/file1`），而相对路径则一般直接开始写路径或以 `.` 开头（例如 `docs/file1` 或 `./docs/file1`）。

> 注意：Linux 的路径中目录的分隔符是 `/`，而 Windows 的路径中的分隔符是 `\`，不要搞混了。

## 3. 命令中的特殊字符

下面来看一个比较长的命令例子。

```bash
cat /etc/apt/sources.list \
| gzip \
> sources.list.gz
```

上面这个例子中，第 1、2 行末尾的 `\` 字符用于**转义**其后的换行符，使 Bash 忽略换行符。Bash 实际读到的命令就像下面这样。

> 注意：用 `\` 转义换行符时，其后必须直接接换行，不能再写任何东西。

```bash
cat /etc/apt/sources.list | gzip > sources.list.gz 
```

`gzip` 命令前的 `|` 是**管道符**，它的意思是将前面的命令的**标准输出**（STDIN），作为后一个命令的**标准输入**（STDOUT）。

`gzip` 命令后的 `>` 为**输出重定向**符号，它的意思是将命令的标准输出写入指定名字的文件中（这里是 `sources.list.gz`）。

相比于将一个命令的输出写入一个文件，再用另一个命令读取这个文件，用管道符将命令的输入和输出直接串起来能很大程度上提高程序的运行效率（前一个命令产出的数据直接在内存中被交给了下一个命令）。

## 4. 安装软件

我们所安装的 Debian 自带有一个名为 `apt` 的软件包管理器。

`apt` 会将软件安装到系统的目录下，所以一般需要搭配 `sudo` 命令获取管理员权限。在使用 `apt` 安装软件前，通常需要用如下命令更新软件包索引。

```bash
sudo apt update
```

由于使用了 `sudo` 进行提权，在键入命令后，会提示输入管理员密码（即安装 WSL 时设置的密码）。

更新完索引后，便可安装软件。以下是一个用 `apt` 安装 `git` 的例子。

```bash
sudo apt install git
```

之后会询问是否确认安装，键入 `y` 即可。

可使用 `sudo apt remove` 命令来删除软件。

绝大部分常用的软件都可用 `apt install` 来安装。然而对于一些生信分析中用到的软件，可能并没有打包到 `apt` 上，这时便需要仔细阅读软件作者提供的安装说明进行安装。常见的安装方法有直接下载可执行二进制文件、下载源代码编译，使用 conda 安装等。

## 5. Bash 编程

除了可以直接输入命令执行，我们也可以在 Bash 里运行一些 Bash 脚本语句，方便我们一次性处理多个文件。

网上有很多的教程可供深入学习（例如：[阮一峰的 Bash 脚本教程](https://wangdoc.com/bash/)），这里只介绍一些基本的语法。

### 变量与字符串

首先，定义一个变量 `x`，它的值是一个字符串 `"world"`（注意：等号的左右不能有空格）。

```bash
x="world"
```

我们可以用 `$` 符号提取这个变量的值。

```bash
echo hello, $x!
# hello, world!
```

Bash 中有两种字符串：用双引号引起来的，和用单引号引起来的（注意：是半角的英文引号）。

双引号字符串内可用 `$` 符号将变量替换为对应的值，而单引号则不会进行这些操作。

```bash
echo "hello, $x!"
# hello, world!
echo 'hello, $x!'
# hello, $x!
```

有时，为了避免变量名后面的字母引发歧义，我们会在用 `$` 替换变量值时用花括号 `{}` 将其括起来。

```bash
xy="me"
echo "hello, $xy!"
# hello, me!
echo "hello, ${x}y!"
# hello, worldy!"
```

`$` 符号除了可以替换变量的值，也可以用来替换命令的输出结果。

```bash
touch file1 file2 file3
echo "Files in this directory: $(ls -m)"
# Files in this directory: file1, file2, file3
```

要删除一个变量，可使用 `unset` 命令。

```bash
unset xy
```

### 数组与循环

Bash 中可以用 `()` 定义数组。普通的索引数组的定义方法如下所示。

```bash
index_array=(a b c)
```

对于索引数组，元素的索引是从 0 开始的整数。我们可以用 `${}` 加上 `[]` 和元素的索引，从数组中提取元素。

```bash
echo "1st: ${index_array[0]}, 2nd: ${index_array[1]}"
# 1st: a, 2nd: b
```

将 `@` 符号放入 `[]` 中，可提取数组的所有元素。

```bash
echo "All items: ${index_array[@]}"
# All items: a b c
```

配合使用 for-in 循环，可依次遍历数组的元素。

```bash
for item in ${index_array[@]}
do
    echo "Iterating item: ${item}."
done
# Iterating item: a.
# Iterating item: b.
# Iterating item: c.
```

Bash 中还有一种可自定义索引名称的数组：关联数组。定义关联数组必须使用 `declare -A` 命令声明。

```bash
declare -A assoc_array
assoc_array=(["A"]=file1 ["B"]=file2 ["C"]=file3)
```

关联数组提取所有元素的方法与索引数组一样，除此之外，还可以在关联数组变量名前添加一个 `!`，以提取所有的索引名。

```bash
echo "All values: ${assoc_array[@]}."
# All values: file3 file2 file1.
echo "All keys: ${!assoc_array[@]}."
# All keys: C B A.
```

搭配 for-in 循环，便可同时遍历数组的元素及其索引了。

```bash
for index in ${!assoc_array[@]}
do
    echo "Key: ${index}, value: ${assoc_array[$index]}"
done
# Key: C, value: file3
# Key: B, value: file2
# Key: A, value: file1
```

### 模式扩展

很多时候，我们要处理的文件的名字中存在重复的部分（例如 `file1`、`file2` 和 `file3`）。

Bash 中有一些模式扩展的方法，以下是一些例子。

```bash
ls  # 列出当前目录下的文件
# file1  file2  file3

touch file{4..6}  # 新建以 file 开头，以 4、5 和 6 结尾的文件

ls
# file1  file2  file3  file4  file5  file6

gzip file*

ls
# file1.gz  file2.gz  file3.gz  file4.gz  file5.gz  file6.gz
```
