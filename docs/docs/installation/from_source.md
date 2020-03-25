---
layout: default
title: Install from Source
parent: Installation
nav_order: 4
---

# Install from Source (advanced)

Please skip this section if you have already installed with ```pip```. This is intended for more advanced users intending to modify the source code of DART-ID.

We recommend using a [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) to run DART-ID from source, so as to avoid any package versioning issues. You can create the same environment as us by using either the provided ```environment.yml``` spec or ```spec-file.txt``` (Mac OSX only). 

```bash
conda create --file environment.yml -n dart python=3.7 && conda activate dart
```

Then download the source code and run the provided shell script:

```bash
git clone https://github.com/SlavovLab/DART-ID && cd DART-ID
./dart_id.sh
```
