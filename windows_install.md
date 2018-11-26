Windows Installation
====================

DART-ID has been tested on Windows 7/8/10, but the installation process is longer, more complex, and prone to errors. 
We strongly recommend running this software on Linux/OSX instead.

## Requirements

* Git source control management -- [download here](https://git-scm.com/downloads)
* Python >= 3.6
    * Strongly recommend the Anaconda distribution of Python -- [download here](https://www.anaconda.com/download/#windows)
    * Miniconda (a version packaged with less stuff) also works too: [download here](https://conda.io/miniconda.html)
    * Vanilla distribution of Python, if you want to start from scratch. [Download here](https://www.python.org/download/releases/3.0/)
* ```pip```
    * ```pip``` should be packaged with Anaconda/Miniconda.
    * If using a vanilla installation, ```pip``` can be [installed here](https://pip.pypa.io/en/stable/installing/)
* ```pystan```
    * On Linux/OSX, pystan in installed with the ```dart_id``` installation command, but in Windows we need to take a few extra steps
    * Follow [the ```pystan``` installation instructions](https://pystan.readthedocs.io/en/latest/windows.html#windows) to get the library to work after installing ```dart_id```

## Installation

Run:

```
pip install git+https://github.com/SlavovLab/DART-ID --user --no-cache-dir
``` 

After installation, you may need to add your user scripts directory to your PATH environment variable, 
in order to run it from the command line. You can skip this step by navigating to the folder of the ```dart_id.exe```
executable, but this is tedious and can be avoided.

Adding the pip user scripts directory to the system PATH is different between Windows versions, so that needs to be looked up.
Once the user scripts folder is loaded onto the PATH, restart your shell and try to run ```dart_id -h``` to confirm that
the PATH is loaded correctly and that the installation worked.

## Running DART-ID

DART-ID can be run in the same way, with the same command structure, as it is in Linux/OSX. 
Make sure to account for path delimiter changes ("\\" instead of "/" to separate folders) in both the command line and
the configuration file.

