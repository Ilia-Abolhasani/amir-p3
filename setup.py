from setuptools import setup, find_packages

SETUP_METADATA = {
    "name": "ct_analizer",
    "version": "1.0.0",
    "description": "",
    "url": "https://github.com/",
    "author": "",
    "author_email": "",    
    "license": "MIT",
    "packages": find_packages(where="src"),            
    "install_requires": ["numpy~=1.20"],
    "setup_requires": ["setuptools~=58.0"],
    "python_requires": ">=3.6",
    "classifiers":[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    }

setup(**SETUP_METADATA)