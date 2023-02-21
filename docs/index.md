---
layout: default
title: Home
nav_order: 1
description: 'Documentation for the DART-ID tool, from the paper: "DART-ID Increases Single-Cell Proteome Coverage".'
permalink: /
---

# DART-ID

{: .fs-9 }

DART-ID is a computational tool that leverages reproducible retention times to increase confident peptide identifications in LC-MS/MS proteomics experiments.
{: .fs-6 .fw-300 }

This method is particularly useful for methods relying on MS2-based quantification, i.e., isobaric tags such as TMT.

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/SlavovLab/DART-ID){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2} [PLOS Comp Biol Article](https://doi.org/10.1371/journal.pcbi.1007082){: .btn .fs-5 .mb-4 .mb-md-0 } [bioRxiv Preprint](https://www.biorxiv.org/content/10.1101/399121v3){: .btn .fs-5 .mb-4 .mb-md-0 }

![GitHub release](https://img.shields.io/github/release/SlavovLab/DART-ID.svg)
![PyPI](https://img.shields.io/pypi/v/dart-id.svg)
![PyPI - Downloads](https://img.shields.io/pypi/dm/dart-id.svg)
![GitHub](https://img.shields.io/github/license/SlavovLab/DART-ID.svg)

---

## Getting started

### Dependencies

DART-ID requires Python = 3.7.6 (64-bit - [miniconda distribution recommended](https://docs.conda.io/en/latest/miniconda.html)), and has been tested on Windows 8 / OSX Mojave 10.14 / Centos 7 / Ubuntu 14.04.

### Installation

DART-ID is available on [PyPI](https://pypi.org/project/dart-id/) and can be installed with [pip](https://pip.pypa.io/en/stable/installing/).

```bash
pip install dart-id
```

### Usage

DART-ID requires a YAML-formatted configuration file to run. An example annotated config file can be found in [example/config_annotated.yaml](https://github.com/SlavovLab/DART-ID/blob/master/example/config_annotated.yaml). You can specify input files and the output folder on the command line, if that's what you prefer.

View the command-line arguments anytime by running: `dart_id -h`.

```
usage: dart_id [-h] [-i INPUT [INPUT ...]] [-o OUTPUT] [-v] [--version] -c
                 CONFIG_FILE

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Input file(s) from search engine output (e.g.,
                        MaxQuant evidence.txt). Not required if input files
                        are specified in the config file
  -o OUTPUT, --output OUTPUT
                        Path to output folder
  -v, --verbose
  --version             Display the program's version
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        Path to config file (required). See
                        example/config_example.yaml
```

### Example runs

An example configuration file can be downloaded from GitHub: [https://github.com/SlavovLab/DART-ID/blob/master/config_files/example_sqc_67_95_varied.yaml](https://raw.githubusercontent.com/SlavovLab/DART-ID/master/config_files/example_sqc_67_95_varied.yaml).

The first few lines of the above configuration file specify the path to the input file:

```yaml
## Input
## ==========================

input:
  - /path/to/SQC_67_95_Varied/evidence.txt
```

You can download the `evidence.txt` file from MassIVE: [ftp://massive.ucsd.edu/MSV000083149/other/MaxQuant/SQC_67_95_Varied/evidence.txt](ftp://massive.ucsd.edu/MSV000083149/other/MaxQuant/SQC_67_95_Varied/evidence.txt).

Then edit the path to the file downloaded, and run the following command:

```bash
dart_id -c config_files/example_sqc_67_95_varied.yaml -o ~/DART_ID/SQC_67_95_varied_20181206
```

The `-o` parameter points to the output folder for DART-ID. You can also specify this path in the config file.

An example analysis of the data and configuration file specified above is available publicly at [ftp://massive.ucsd.edu/MSV000083149/other/Alignments/SQC_varied_20180711_4/](ftp://massive.ucsd.edu/MSV000083149/other/Alignments/SQC_varied_20180711_4/).

---

## About the project

DART-ID is a project developed in the [Slavov Laboratory](https://web.northeastern.edu/slavovlab/) at [Northeastern University](https://www.northeastern.edu/) [Bioengineering](http://www.bioe.neu.edu/), and was authored by [Albert Tian Chen](https://atchen.me), [Alexander Franks](http://afranks.com/) (of [UCSB Statistics and Applied Probability](https://www.pstat.ucsb.edu/)), and [Nikolai Slavov](https://web.northeastern.edu/slavovlab/).

The article for DART-ID is freely available on PLOS Computational Biology: [https://doi.org/10.1371/journal.pcbi.1007082](https://doi.org/10.1371/journal.pcbi.1007082).

The preprint is also available on bioRxiv: [https://www.biorxiv.org/content/10.1101/399121v3](https://www.biorxiv.org/content/10.1101/399121v3).

Contact the authors by email: [nslavov\{at\}northeastern.edu](mailto:nslavov@northeastern.edu).

### License

DART-ID is distributed by an [MIT license](https://github.com/SlavovLab/DART-ID/blob/master/LICENSE.txt).

### Contributing

Please feel free to contribute to this project by opening an issue or pull request in the [GitHub repository](https://github.com/SlavovLab/DART-ID).

### Data

All data used for the manuscript is available on [UCSD's MassIVE Repository](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ed5a1ab37dc34985bbedbf3d9a945535)

### Figures/Analysis

Scripts for the figures in the DART-ID manuscript are available in a separate GitHub repository, [https://github.com/SlavovLab/DART-ID_2018](https://github.com/SlavovLab/DART-ID_2018)
