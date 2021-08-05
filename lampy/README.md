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

## How to use

First import the package into the python interpreter (we suggest using `ipython` for an inline analysis)
```
import lampy as lpy
```
which will provide a starting help message:
```
Light and fast ALaDyn's data Manipulation in PYthon (LAMPy)

This package is developed by the ALaDyn Collaboration

You can read, plot and elaborate ALaDyn's output files.
To import the simulation folder, type:

>>> sim=lampy.Simulation(path)

where 'path' is the relative or absolute path of interest.
To read more, please use the command:

>>> help(lampy.Simulation)
```

To load the data into a simulation object and start the manipulation, use the command:
```
s=lpy.Simulation(path)
```
Be sure that the folder specified by `path` contains the output folders of the simulation (typically called `0000`, `0001`, etc...) and the `input_**.nml` or `input_**.json` file that contains the simulation infomation that `lampy` uses to decode the data.

Once the simulation is loaded, three main classes will be available to manipulate the data:

- `Field` provides methods to read, plot and extract data related to the electromagnetic fields
- `Particles` provides methods to read, plot and extract data related to the simulation particles, if the corresponding outputs have been produced
- `Diagnostics` provides methods to read the diagnostic files produced in the folder `diagnostics`

For more details on any of those classes and for the complete list of the methods they provide, please use the corresponding help command:
```
>>> help(s.Field)
>>> help(s.Particles)
>>> help(s.Diagnostics)
```
### Examples

---
**NOTE** : `ALaDyn` electromagnetic fields are all normalized to `e/mc^2`, which means that are in units of `\mu m^-1`.
To obtain the fields in `TV/m`, multiply for `m_e [MeV] = 0.511`, to obtain the fields in plasma units (*e.g* in units of the wavebreaking field), divide them by the plasma wavenumber.
---

We can plot the longitudinal electric field on axis using:
```
>>> s.Field.lineout('Ex', 0)
```
where `'Ex'` is the field chosen and `0` is the plotting time **in microns**.
The function supports different arguments (like off-axis plot, different plotting orientation, normalization, etc...) which you can read about using the help function
```
>>> help(s.Field.lineout)
```

It is also possible to produce a 2d color map (heatmap) of the fields,
using the function
```
>>> s.Field.map_2d('Ex', 0)
```
which again supports other arguments, available in the help function.

If the particle files are available, you can produce a phase
space histogram of a given phase space component via the command
```
>>> s.Particles.histogram('phase_space_electrons', 0)
```
where you have to specify an available phase space as a first argument. You can read more about it in the help function.