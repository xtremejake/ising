#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 01:00:27 2017

@author: jacobmarold
"""
import os
import sys

BASE_DIR = os.path.dirname(__file__)
PACKAGE_DIR = os.path.dirname(BASE_DIR)
ASSETS_DIR = os.path.join(PACKAGE_DIR, 'assets')

if __name__ == "__main__":
    print ASSETS_DIR
