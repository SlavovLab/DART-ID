DART-ID
=============

Intro
-----

The DART-ID code goal is to make the process as simple as possible. Run ```dart_id```, point it at the inputs, and expect an additional PEP_new column at the end. There are still some parameters that you will need to tweak manually as the default settings have not been fully generalized yet.

# Installation

```
pip install git+https://github.com/SlavovLab/DART-ID --user --no-cache-dir
```

Uninstall the package with:

```
pip uninstall dart_id
```

The installation requires Python >= 3.6, and has been tested on Windows 8 / OSX Mojave 10.14 / Ubuntu 14.04. Ubuntu/Linux instructions can be found in linux_install.md

# Known bugs:

- I/O issues when outputting to a Google File Stream folder. Output to a local folder instead, and then drag it into Google File Stream later

Usage
----------

View parameters anytime yourself by typing "dart_id -h".

Input and output files are optional in the command line, and can be specified in the config file instead, if that's what you prefer.

An example config file, annotated and unannotated, can be found at ```example/config.yaml```

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

Example runs
============

```
dart_id -i /path/to/search_engine_output/evidence.txt -o /path/to/output/folder -v -c /path/to/config/file.yaml
```

If, for example, the input files are already listed out in the config file, then you can skip the ```-i``` option.

```
dart_id -c /path/to/config/file.yaml -o /path/to/output/folder
```
