from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import os

about = dict()
with open(os.path.join('lampy', '__about__.py')) as fp:
    exec(fp.read(), about)

__version__ = about['__version__']

with open('README.md') as f:
    long_description = f.read()

extensions = [Extension(
          name='read_field',
          sources=[os.path.join('lampy', 'fastread', 'read_field.pyx')]),
          Extension(
          name='read_phase_space',
          sources=[os.path.join('lampy', 'fastread', 'read_phase_space.pyx')]),
          Extension(
          name='read_tracking',
          sources=[os.path.join('lampy', 'fastread', 'read_tracking.pyx')])
]

lib_read_binary = ('lib_read_binary',
                   {'sources': [os.path.join('lampy', 'fastread', 'clibs',
                                             'lib_read_binary.c')]})
lib_read_phase_space = ('lib_read_phase_space',
                        {'sources':
                         [os.path.join('lampy', 'fastread', 'clibs',
                                       'lib_read_phase_space.c')]})
lib_read_tracking = ('lib_read_tracking',
                     {'sources':
                      [os.path.join('lampy', 'fastread', 'clibs',
                                    'lib_read_tracking.c')]})

setup(
    version=__version__,
    ext_modules=cythonize(extensions),
    ext_package='compiled_cython',
    libraries=[lib_read_binary, lib_read_phase_space, lib_read_tracking],
    include_dirs = [numpy.get_include()]
)
