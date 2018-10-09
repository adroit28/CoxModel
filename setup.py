# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:45:20 2018

@author: A1298
"""

from distutils.core import setup

# Package meta-data.
NAME = 'CoxModelSelection'
DESCRIPTION = 'Choose a model to use using different selection criterion'
URL = 'https://github.com/me/myproject'
EMAIL = 'suraj_bansal@outlook.com'
AUTHOR = 'Suraj Bansal'
REQUIRES_PYTHON = '>=2.7.0'
VERSION = 0.1

# Packages required for this module to be executed
REQUIRED = [
    'lifelines', 'multiprocessing', 'pandas', 'numpy', 'copy', 'inspect', 'math'
]

# Import the README and use it as the long-description.
try:
    long_description = open('README.txt').read()
except:
    long_description = DESCRIPTION

setup(
    name = NAME,
    version = VERSION,
    packages = ['CoxModelSelection',],
    license = 'MIT',
    long_description = long_description,
    install_requires = REQUIRED,
    include_package_data = True
)