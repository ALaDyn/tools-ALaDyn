# LAMPy

LAMPy (**L**ight and fast **A**ladyn data **M**anipulation in **Py**thon) is a python package that provides modules to access and plot ALaDyn output files.

The idea behind its development is to provide a selfconsistent and portable environment, allowing any user to elaborate ALaDyn simulations results.
Other possibility are also available (C++ data converters, gnuplot analysis, IDL routines), but they are mostly single file analyzers. With LAMPy, we would like to offer a fast and modern way to perform a complete *inline* data analysis, taking advantage of the simplicity of the Python language which requires a minimum effort even for a beginner user. 

## Package

LAMPy is still under development, which means problems may occur for some untested output configuration. The recommended Python version for LAMPy is Python 3.7. Tests on previous versions are ongoing.

## Installation

First, you should get a Python distribution, such as [Anaconda](https://www.anaconda.com/distribution/#download-section). You could also use a system installation of Python, however, due to the portability of Anaconda and it's ability to manage different environments, we recommend it as a preferred option, in particular for beginner users.

Create a Python 3.7 environment, if it's not the default one', by typing

`conda create -n py3 python=3.7`

To install LAMPy, you can build it from the source code, or you can download the latest stable version from Pypi, installing it via pip

`pip install lampy`

## Build LAMPy from source

To build LAMPy from source, you must run the following commands

- `python setup.py build_clib` that compiles the C libraries needed to read the data files using a system C compiler,
- `python setup.py develop` that builds the cython and python sources and installs the package in the `$PYTHONENV`.

