from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import os

about = dict()
with open("lampy/__about__.py") as fp:
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

requirements = ['numpy', 'matplotlib', 'Cython', 'scipy']
setup(
    name='lampy',
    version=__version__,
    description="Python suite to access, manipulate and plot ALaDyn's datas",
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages('.'),
    author='ALaDyn Collaboration',
    author_email='dterzani@lbl.gov',
    maintainer='Davide Terzani',
    maintainer_email='dterzani@lbl.gov',
    license='GNU GPLv3',
    url='https://github.com/ALaDyn/ALaDyn',
    classifiers=[
        'Programming Language :: Python',
        'Natural Language :: English',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Physics',
        'Programming Language :: Python :: 3.7'],
    ext_modules=cythonize(extensions),
    install_requires=[requirements],
    ext_package=os.path.join('lampy', 'compiled_cython'),
    libraries=[lib_read_binary, lib_read_phase_space, lib_read_tracking],
    include_dirs=[numpy.get_include()]
)
