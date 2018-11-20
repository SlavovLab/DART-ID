# Linux/Ubuntu Install Instructions

## Getting Python 3.6

DART-ID is tested on Ubuntu 14.04 (Trusty64), but should run on most distros. 

The most important step is getting Python >=3.6, as DART-ID and some of its dependencies will fail to build on earlier versions.

With Ubuntu 14.04/16.04, Python3.6 is not installed by default, and to install you can follow instructions 
[here](https://askubuntu.com/questions/865554/how-do-i-install-python-3-6-using-apt-get). 
Ubuntu 16.10 and above should have Python3.6 installed by default.

## Getting `pip` or `conda`

DART-ID does not require the Anaconda/Miniconda environment, but it is useful to have for other reasons.
To get Anaconda/Miniconda, follow the instructions [here](https://conda.io/miniconda.html)
To get `pip`, follow the instructions [here](https://pip.pypa.io/en/stable/installing/)

## Installing DART-ID and dependencies

```python3 -m pip install git+https://github.com/SlavovLab/DART-ID --user```

-----------

### Misc. Notes

If using a Python installation from `apt` on Ubuntu, make sure to also install `python3-dev` 
as the `pystan` dependency will need to access some of those tools.
