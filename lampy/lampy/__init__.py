from .Simulation import Simulation
from .__about__ import __version__

__doc__ = """
Light and fast ALaDyn's data Manipulation in PYthon (LAMPy)

This package is developed by the ALaDyn Collaboration

You can read, plot and elaborate ALaDyn's output files.
To import the simulation folder, type

>>> sim=lampy.Simulation(path)

where 'path' is the relative or absolute path of interest.
"""

print("LAMPy v"+__version__)
print(__doc__)
