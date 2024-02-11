---
layout: post
title:  "用 Captura 录屏"
date:   2021-03-22 23:00:00 +0800
categories: 教程
---

[Captura](https://mathewsachin.github.io/Captura/) 是一个免费开源的录屏软件。它功能强大，但需要搭配 [FFmpeg](https://www.ffmpeg.org/) （一个几乎万能的音视频流处理软件）使用。

本文记录了 Captura 安装过程和使用时的注意事项。

## 下载软件

首先到 [Captura 的官网](https://mathewsachin.github.io/Captura/)下载 Portable 版的 Captura 软件本体。

![download-captura](/assets/images/screen-recording/download-captura.png)

然后，从 [FFmpeg 的官网](https://www.ffmpeg.org/)开始，找到并下载 [gyan.dev](https://www.gyan.dev/ffmpeg/builds/)  编译的 FFmpeg。

![download-ffmpeg](/assets/images/screen-recording/download-ffmpeg.png)

> 懒得到网页上找下载链接？右键另存为吧：[Captura-Portable.zip](https://github.com/MathewSachin/Captura/releases/download/v8.0.0/Captura-Portable.zip)，[ffmpeg-4.3.2-2021-02-27-full_build.7z](https://github.com/GyanD/codexffmpeg/releases/download/4.3.2-2021-02-27/ffmpeg-4.3.2-2021-02-27-full_build.7z)

## 设置一下

都下载下来后，解压两个压缩包。

![extracted](/assets/images/screen-recording/decompressed.png)

双击从 Captura-Portable.zip 解压出来的 captura.exe 即可启动该软件。但启动之后并不能直接使用。我们得把 Captura 和 FFmpeg 这两个软件连接起来。

单击右上角的齿轮，然后修改 FFmpeg 的位置为刚解压出来的 ``bin`` 文件夹。

![set-ffmpeg-location](/assets/images/screen-recording/set-ffmpeg-location.png)

现在，回到 Captura 主界面，点击左上角的红圆点就可以录屏了。

## 使用 Captura 时的注意事项

Captura 功能强大，但也有一些 BUG 需要注意。

首先，录屏前需要点击主界面中下部的麦克风按钮设置音频，如果不设置，录出来的是没有声音的。

![set-audio](/assets/images/screen-recording/set-audio.png)

如果电脑插着耳机，那么在音频设置中要选择相应的耳机的音频，而且在录制过程中，**耳机不能拔掉**，否则即使再连上，之后的录屏也是没有声音的。

### 参考

- Captura 官网: <https://mathewsachin.github.io/Captura/>
- FFmpeg 官网: <https://www.ffmpeg.org/>
- gyan.dev 编译的 FFmpeg: <https://www.gyan.dev/ffmpeg/builds/>