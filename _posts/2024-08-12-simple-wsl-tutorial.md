---
layout: post
title:  "简单的 WSL 教程"
categories: 教程
---

## 安装 WSL

首先需要开启 Windows 的相关功能。按快捷键 `win` + `Q`，搜索 `启用或关闭 Windows 功能`。打开后，找到 `适用于 Linux 的 Windows 子系统` 和 `虚拟机平台` 这两个选项，勾选上后确定退出，然后等待系统完成更新后重启。

接着是安装一个 Linux 发行版。按快捷键 `win` + `Q`，搜索 `PowerShell` 并打开，然后依次输入以下命令：

```powershell
wsl --update
wsl --install -d Debian
```

在安装完成后，会自动打开 WSL，根据提示输入 WSL 的用户名和密码后即可进入 WSL 的终端。

若退出后想要重新进入，在开始菜单中找到 Debian 重新进入即可。

## Linux 入门

网上的 Linux 入门教程很多，所以这里只做基本的介绍。

### 基础文件操作

刚一进入 WSL，所处的目录被称作“家目录”。在这个目录下，让我们先练习以下常见的文件操作：

0. 键入 `pwd` 可显示当前目录的绝对路径
1. 键入 `ls` 可列出当前目录下的内容（目前是空的）
2. 键入 `mkdir Documents` 可新建一个名为 `Documents` 的目录
3. 键入 `touch file1` 可新建一个名为 `file1` 的空文件
4. 键入 `cp file1 file2` 可将 `file1` 复制为 `file2`
5. 再键入 `ls`，当前目录下就有 `Documents`、`file1` 和 `file2` 这三个东西了
6. 键入 `mv file1 file2 Documents` 将 `file1` 和 `file2` 都移动进 `Documents` 目录
7. 再键入 `ls`，当前目录下只剩 `Documents` 目录了
8. 使用 `cd Documents` 切换当前工作目录至 `Documents` 目录下
9. 再键入 `ls`，便可列出 `Documents` 下的 `file1` 和 `file2` 这两个文件
10. 键入 `cd ..` 可退回上一层目录
11. 键入 `tar -cf archive.tar.gz Documents` 可将 `Documents` 压缩为一个名为 `archive.tar.gz` 的压缩包
12. 键入 `rm -r Documents` 可删除 `Documents` 目录及其内部的所有文件
13. 再键入 `ls`，当前目录下只剩 `archive.tar.gz` 压缩包了
14. 键入 `tar -xf archive.tar.gz`，解压 `archive.tar.gz`
15. 再键入 `ls`，`Documents` 目录已经被解压出来了

总结一下，在以上的练习中，我们了解了以下命令的功能：

- `pwd` - 显示当前目录的绝对路径
- `ls` - 列出目录中的内容
- `mkdir` - 新建目录
- `touch` - 新建文件
- `cp` - 复制文件
- `mv` - 移动文件
- `cd` - 切换工作目录
- `tar` - 压缩和解压缩
- `rm` - 删除文件

你可以在任意一个大语言模型 AI 处获得这些命令的详细使用方法及示例，所以在此不做更多的介绍了。

### 安装软件

我们所安装的 Debian 自带有一个名为 `apt` 的软件包管理器。

`apt` 会将软件安装到系统的目录下，所以一般需要搭配 `sudo` 命令获取管理员权限。在使用 `apt` 安装软件前，通常需要用如下命令更新软件包索引。

```bash
sudo apt update
```

更新完索引后，便可安装软件。以下是一个用 `apt` 安装 `git` 的例子。

```bash
sudo apt install git
```

由于使用了 `sudo` 进行提权，在键入命令后，会提示输入管理员密码（即安装 WSL 时设置的密码）。

之后会询问是否确认安装，键入 `y` 即可。

若要卸载软件，使用 `sudo apt remove` 命令可用 `apt` 删除软件。

绝大部分常用的软件都可用 `apt install` 来安装。然而对于一些生信分析中用到的软件，可能并没有打包到 `apt` 上，这时便需要仔细阅读软件作者提供的安装说明进行安装。常见的安装方法有直接下载可执行二进制文件、下载源代码编译，使用 conda 安装等。

### 路径

在 Linux 中输入命令时，常涉及到文件的路径，准确无误地输入路径是成功运行命令所必需的。

Linux 中存在有两种路径：绝对路径和相对路径。绝对路径是相对于系统根目录的路径，而相对路径是相对于当前工作目录的路径。

绝对路径一般以 `/` 开头（例如 `/home/user/Documents/file1`），而相对路径则一般直接开始写路径（例如 `Documents/file1`）。

注意：Linux 的路径中目录的分隔符是 `/`，而 Windows 的路径中的分隔符是 `\`，不要搞混了。

除此之外，还有一些特殊路径的符号，例如 `./` 是当前工作目录，`../` 是上级目录，而在 Shell 中，`~/` 可被扩展成家目录。

以下是一些例子：

- `/home/user/` - 家目录的绝对路径
- `./` - 当前工作目录的相对路径
- `file1`（或 `./file1`） - 当前目录下 `file1` 文件的相对路径
- `~/Documents/file1` 家目录下 `file1` 文件的绝对路径，等于 `/home/user/Documents/file1`

### Bash 基本语法

刚才用来输入命令的终端是一个 Bash 终端。除了可以直接输入命令执行，我们也可以在 Bash 里运行一些 Bash 脚本语句，方便我们一次性处理多个文件。

网上有很多的教程可供深入学习（例如：[阮一峰的 Bash 脚本教程](https://wangdoc.com/bash/)），这里只介绍一些基本的语法。

#### 变量与字符串

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

#### 数组与循环

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

#### 模式扩展

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

## WSL 与 Windows 的协作

### 右键+在终端中打开

大部分时候，我们想要做的是在 WSL 中运行软件，来处理 Windows 文件系统下的数据。

在 WSL 中，Windows 的各个磁盘的映射位于 `/mnt/` 目录下，例如 C 盘就是 `/mnt/c/`，D 盘就是 `/mnt/d/`。但如果每次都要先打开 WSL，再一路 `cd` 到数据存放的目录，还是太麻烦了。

通过安装 Windows Terminal，可以让我们在 Windows 的文件资源管理器中，通过 `右键` + `在终端中打开`，在当前目录下打开 WSL，并切换好工作目录。

Windows 11 默认是安装有 Windows Terminal 软件的，而 Windows 10 也可以在 Microsoft Store 搜索安装。

在安装完 Windows Terminal 后，文件资源管理器的右键菜单里便会出现 `在终端中打开` 的选项，但这时默认打开的是 Windows 的 Powershell。

可在 Windows Terminal 里按 `Ctrl` + `,`，在 `设置 > 启动 > 默认配置文件` 处将默认启动的终端修改成 Debian，这样默认打开的就是 WSL 了。

### 在 WSL 中启动 Windows 上的软件

在 WSL 中，也可使用一些命令来启动安装在 Windows 上的软件。

例如在 WSL 中键入 `explorer.exe .`（注意要有个 `.`）便可用文件资源管理器中打开当前工作目录。

如果 Windows 下安装有 [Visual Studio Code](https://code.visualstudio.com/)，在 WSL 中键入 `code .`，便可在当前目录启动 VSCode。而 VSCode 里也有名为微软开发的 WSL 的插件，安装后可进一步方便在 WSL 中的开发工作。
