---
layout: post
title:  "在 WSL 中安装 RStudio Server"
categories: 笔记
---

本文记载了以下内容：
- 在 WSL-Ubuntu-20.04 中安装 R 和 RStudio Server；
- 通过各种设置以达到“优雅地”使用 RStudio Server 的目的。

## 安装 R 和 RStudio Server

首先按照添加 TUNA 的 CRAN 软件源。

```bash
sudo vim /etc/apt/sources.list.d/r-tuna.list
```

在其中添加以下内容。

```
deb https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/linux/ubuntu focal-cran40/
```

保存退出后在终端中执行以下内容。

```bash
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt-get update
```

这样 TUNA 的 CRAN 软件源就添加完了。

然后安装 R，以及其他一堆 RStudio Server 依赖的软件。

```bash
sudo apt install -y r-base r-base-core r-recommended r-base-dev gdebi-core build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
```

最后，安装 RStudio Server。

```bash
wget https://rstudio.org/download/latest/stable/server/bionic/rstudio-server-latest-amd64.deb
sudo gdebi rstudio-server-latest-amd64.deb
```

安装完成后通过以下命令即可启动 RStudio Server。

```bash
sudo rstudio-server start
```

然后用浏览器打开 <http://127.0.0.1:8787/>。

![浏览器中的 RStudio](/assets/images/install-rstudio-server-in-wsl/rstudio-server-installed.png)

当当~！一个在浏览器里运行的 RStudio 就装好了。🎉

用户名和密码就是 WSL 的用户名和密码。

## 关闭 RStudio Server 的自动登出

在刚才的登录界面中有这样一行小字。

```
You will automatically be signed out after 60 minutes of inactivity.
```

翻译：“超过 60 分钟没进行操作，就会自动登出。”

RStudio Server 原本是一个部署在服务器上的软件，所以出于安全考虑，默认开启这个自动登出的功能。

但对于我们来说，每次都要输密码太麻烦了，让我们关掉它。

编辑 RStudio Server 的配置文件 ``/etc/rstudio/rserver.conf``

```bash
sudo vim /etc/rstudio/rserver.conf
```

在其中添加以下语句。

```txt
auth-timeout-minutes=0
auth-stay-signed-in-days=30
```

保存退出后重启 RStudio Server。

```bash
sudo rstudio-server restart
```

再重新通过浏览器登录 RStudio Server 时，小字就消失了。

## 创建一个 chrome 应用

现在的 RStudio Server 看上去就是个运行在浏览器里的 RStudio，是一个网页，并不像一个应用。

我们可以用 chrome 来让它看上去更像一个独立的 RStudio 应用。

首先用 chrome 打开 <http://localhost:8787/> 然后点击右上角的菜单，选择“更多工具” -> “创建快捷方式”，勾选“在窗口中打开”，确定后，开始菜单和桌面上就会出现 RStudio Server 的快捷方式了。

少了标签栏和地址栏，看上去就更像个应用了。

## 在 WSL 中一键启动 RStudio Server

现在要使用 RStudio Server，需要以下步骤

1. 启动 WSL，运行 ``sudo rstudio-server start``
2. 点击桌面上的 RStudio Server 快捷方式，启动网页应用。

要想在 WSL 中用一个命令直接启动 RStudio Server，我们可以在 ``~/.bashrc`` 文件末尾添加个函数。

```bash
vim ~/.bashrc
```

```bash
# RStudio Server
rstudio(){
  sudo rstudio-server start
  '/mnt/c/Program Files/Google/Chrome/Application/chrome_proxy.exe' \
    --profile-directory=Default \
    --app-id=cfmligackgegaglpgmjfndpcappofaie  # 这里换成自己的
}
```

注意：其中的 ``--app-id`` 参数因人而异，是通过右击 RStudio Server 的快捷方式，点击“属性”，然后从“目标”栏的末尾复制的。

![app-id](/assets/images/install-rstudio-server-in-wsl/app-id.png)

保存退出后重启个终端，再键入 ``rstudio``，输入完 ``sudo`` 要求的密码后，RStudio Server 的网页应用就直接弹出来了。

## 我连 ``sudo`` 的密码都不想输……

新增一个 ``sudo`` 的配置文件。

> ``<用户名>`` 换成自己的用户名。

```bash
sudo vim /etc/sudoers.d/<用户名>
```

在里面写上这一句。

```bash
<用户名> ALL=(ALL) NOPASSWD:ALL
```

保存退出后，再使用 ``sudo`` 命令就不需要密码了。

## 参考

- TUNA 的 CRAN 镜像使用帮助：<https://mirrors.tuna.tsinghua.edu.cn/help/CRAN/>
- RStudio Support 中的安装教程：<https://support.rstudio.com/hc/en-us/articles/360049776974-Using-RStudio-Server-in-Windows-WSL2>
- 如何关闭 RStudio Server 的自动登出：<https://github.com/rstudio/rstudio/issues/5449#issuecomment-637586731>
- How to disable WSL password：<https://superuser.com/questions/1470035/how-to-disable-wsl-password)>

