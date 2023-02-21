---
layout: default
title: Windows Installation
parent: Installation
nav_order: 1
---

# Windows Installation

DART-ID has been tested on Windows 7/8/10. However, we recommend running this software on Linux/OSX instead, if possible.

## Requirements

- Git source control management -- [download here](https://git-scm.com/downloads)
- Python = 3.7.6
  - We strongly recommend using the Miniconda distribution: [download here](https://conda.io/miniconda.html)
  - The Anaconda distribution works as well but its 700+ packages are not strictly necessary -- [download here](https://www.anaconda.com/download/#windows)
  - Vanilla distribution of Python, if you want to start from scratch. [Download here](https://www.python.org/download/releases/3.0/)
- `pip`
  - `pip` should be packaged with Anaconda/Miniconda.
  - If using a vanilla installation, `pip` can be [installed here](https://pip.pypa.io/en/stable/installing/)

## Installation

Run:

```bash
pip install dart-id
```

After installation, you may need to add your user scripts directory to your PATH environment variable,
in order to run it from the command line. You can skip this step by navigating to the folder of the `dart_id.exe` executable, but this is tedious and can be avoided.

Adding the pip user scripts directory to the system PATH is different between Windows versions, so that needs to be looked up. Once the user scripts folder is loaded onto the PATH, restart your shell and try to run `dart_id -h` to confirm that the PATH is loaded correctly and that the installation worked.

## Running DART-ID

DART-ID can be run in the same way, with the same command structure, as it is in Linux/OSX. Make sure to account for path delimiter changes ("path\\to\\folder" instead of "path/to/folder" to separate folders) in both the command line and the configuration file.
