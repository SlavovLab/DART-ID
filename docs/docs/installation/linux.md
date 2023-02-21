---
layout: default
title: Linux Installation
parent: Installation
nav_order: 2
---

# Linux Installation

## Getting Python 3.7

DART-ID is tested on Ubuntu 14.04 (Trusty64) and CentOS 7, but should run on most distros.

The most important step is getting Python = 3.7.6, as DART-ID and some of its dependencies will fail to build on earlier versions.

With Ubuntu 14.04, Python3 is not installed by default, and to install you can follow instructions
[here](https://askubuntu.com/questions/865554/how-do-i-install-python-3-6-using-apt-get).
Ubuntu 16.10 and above should have Python3 installed by default.

## Getting `pip` or `conda`

DART-ID does not require the Anaconda/Miniconda environment, but it is useful to have for other reasons. To get Anaconda/Miniconda, follow the instructions [here](https://conda.io/miniconda.html)

To get `pip`, follow the instructions [here](https://pip.pypa.io/en/stable/installing/)

## Installing DART-ID and dependencies

```bash
python3 -m pip install dart-id
```

If you are working on a managed machine or computing cluster, you may have to install to your user directory instead of the system python directory. Pass the `--user` flag to the install command:

```bash
python3 -m pip install dart-id --user
```
