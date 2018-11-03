#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 12:16:41 2018

@author: jacobmarold
"""

import setuptools
from os import path

FILE_DIR = path.dirname(path.abspath(__file__))
README = path.join(path.dirname(FILE_DIR), 'README.md')

with open(README, "r") as fh:
    long_description = fh.read()

setuptools.setup(name='ising',
                 version='0.1',
                 description='For global fitting of systems represented by 1D Ising Models',
                 url='https://github.com/xtremejake/ising',
                 author='xtremejake',
                 author_email='xtremejake.usa@gmail.com',
                 license='Apache-2.0',
                 #packages=['ising'], # packages to include
                 #install_requires=[] # specifies dependencies of python packages in pip
                 #dependency_links=[] # for dependencies that require a weblink
                 #test_suite='nose.collector',
                 #tests_require=['nose'], # for unit testing
                 #include_package_data=True # allows things to be included which aren't .py files - put the cTPR and DPR data here
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 packages=setuptools.find_packages(),
                 zip_safe=True,
                 classifiers=[
                    "Programming Language :: Python :: 3",
                    "License :: OSI Approved :: Apache Software License",
                    "Operating System :: OS Independent",
                ],
                )   