import sys
from setuptools import setup, find_packages
from setuptools import Extension
import os

SETUP_METADATA = \
               {
    "name": "ct-analizer",
    "description": "",
    "url": "https://github.com/",
    "author": "",
    "author_email": "",
    "version": "1.0.0",
    "license": "MIT",
    "packages": find_packages(),
    "entry_points": {'console_scripts': [
        'vamb = vamb.__main__:main'
        ]
    },
    "scripts": ['src/concatenate.py'],
    "ext_modules": [Extension("vamb._vambtools",
                               sources=["src/_vambtools.pyx"],
                               language="c")],
    "install_requires": ["numpy~=1.20", "torch~=1.8", "pysam~=0.14"],
    "setup_requires": ['Cython~=0.29', "setuptools~=58.0"],
    "python_requires": "~=3.7",
    "classifiers":[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    }

setup(**SETUP_METADATA)