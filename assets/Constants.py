#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 07:35:23 2017

@author: jacobmarold
"""

gas_constants = {
        'J/mol': {'value': float(8.31459848), 
                  'units': 'J*K**-1*mol**-1'},
        
        'erg/mol': {'value': float(8.31459848*10**7), 
                    'units': 'erg*K**-1*mol**-1'},
        
        'amu/mol': {'value': float(8.31459848*10**-3), 
                    'units':'amu*(km/s)**2*K**-1'},
        }