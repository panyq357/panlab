---
layout: post
title:  "简单的 WSL 教程"
categories: 教程
---

## 安装 WSL

首先需要开启 Windows 的相关功能。按快捷键 `win` + `Q`，搜索 `启用或关闭 Windows 功能`。打开后，找到 `适用于 Linux 的 Windows 子系统` 和 `虚拟机平台` 这两个选项，勾选上后确定退出，然后等待系统完成更新后重启。

接着是安装一个 Linux 发行版，这里我们安装的发行版是 Debian。按快捷键 `win` + `Q`，搜索 `PowerShell` 并打开，然后依次输入以下命令：

```powershell
wsl --update
wsl --install -d Debian
```

在安装完成后，会自动打开 WSL，根据提示输入 WSL 的用户名和密码后即可进入 Debian 默认的 Bash 终端。

如果你从未接触过 Linux 的 Bash 终端，请阅读[简单的 Bash 教程]({% post_url 2024-08-23-simple-bash-tutorial %})。

若退出后想要重新进入，在开始菜单中找到 Debian 重新进入即可。

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
