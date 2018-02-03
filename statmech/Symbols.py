#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 09:46:06 2017

@author: jacobmarold
"""

import sympy as sp

x = sp.Symbol('x')   
mi = sp.Symbol('mi')

# Interfacial energetic terms
W = sp.Symbol('W') # coupling energy for homopolymer model
Wab = sp.Symbol('Wab') # coupling energy between A and B
Wba = sp.Symbol('Wba') # coupling energy between B and A

# Intrinsic energetic terms


# Gibbs energy terms
dGAB = sp.Symbol('dGAB')
dGBA = sp.Symbol('dGBA')
