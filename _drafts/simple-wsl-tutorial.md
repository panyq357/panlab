---
layout: post
title:  "简单的WSL教程"
---

## 安装 WSL

首先需要开启 Windows 的相关功能。按快捷键 `win` + `Q`，搜索 `PowerShell`，找到 `适用于 Linux 的 Windows 子系统` 和 `虚拟机平台` 这两个选项，勾选上后确定退出，然后等待系统完成更新后重启。

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

以下是一个用 `apt` 安装 `git` 的例子。

```bash
sudo apt update
sudo apt install git
```

键入命令后，会提示输入密码（即安装 WSL 时设置的密码），以及询问是否确认安装（键入 `y` 即可）。

命令开头的 `sudo` 意为提升权限至超级管理员，这对“用 `apt` 安装软件”这种修改系统内部文件的操作来说是必须的。

使用 `sudo apt remove` 命令可用 `apt` 删除软件。

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

## WSL 与 Windows 的协作

大部分时候，我们想要做的是在 WSL 中运行软件，来处理 Windows 文件系统下的数据。

在 WSL 中，Windows 的各个磁盘的映射位于 `/mnt/` 目录下，例如 C 盘就是 `/mnt/c/`，D 盘就是 `/mnt/d/`。但如果每次都要先打开 WSL，再一路 `cd` 到数据存放的目录，还是太麻烦了。

通过安装 Windows Terminal 可以解决这个问题。Windows 11 默认是安装有 Windows Terminal 软件的，而 Windows 10 也可以在 Microsoft Store 搜索安装。

在安装完 Windows Terminal 后，文件资源管理器的右键菜单里便会出现“在终端中打开”的选项，但这时默认打开的是 Windows 的 Powershell。可在 Windows Terminal 里按 `Ctrl` + `,`，在 `设置 > 启动 > 默认配置文件` 处将默认启动的终端修改成 Debian，这样默认打开的就是 WSL 了，而且 WSL 的工作目录也会自动切换到当前目录下了。

在 WSL 中，也可使用一些命令来启动安装在 Windows 上的软件。

例如在 WSL 中键入 `explorer.exe .`（注意要有个 `.`）便可在文件资源管理器中打开当前目录。

如果 Windows 下安装有 [Visual Studio Code](https://code.visualstudio.com/)，在 WSL 中键入 `code .`，便可在当前目录启动 VSCode。而 VSCode 里也有名为微软开发的 WSL 的插件，安装后可进一步方便在 WSL 中的开发工作。
