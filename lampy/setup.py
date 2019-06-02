from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [Extension(
          name="read_field",
          sources=["lampy/fastread/read_field.pyx","lampy/fastread/clibs/lib_read_binary.c"]),
          Extension(
          name="read_phase_space",
          sources=["lampy/fastread/read_phase_space.pyx","lampy/fastread/clibs/lib_read_phase_space.c"])
]

requirements = ['numpy', 'matplotlib', 'Cython', 'scipy']
setup(
    name='lampy',
    version='0.1.dev1',
    description="Python suite to access, manipulate and plot ALaDyn's datas",
    packages=find_packages('.'),
    author='ALaDyn Collaboration',
    maintainer='Davide Terzani',
    maintainer_email='davide.terzani@ino.it',
    license='GNU GPLv3',
    ext_modules=cythonize(extensions),
    ext_package='lampy/compiled_cython',
    url='https://github.com/ALaDyn/ALaDyn',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Physics',
        'Programming Language :: Python :: 3.7'],
    install_requires=[requirements]
)
