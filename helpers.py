#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri June 3 17:18:26 2022

@author: rlm

helpers scripts for all other things
"""
#%
def namestr(obj):
    try:
        return([name for name in globals() if globals()[name] is obj][0])
    except IndexError:
        return([name for name in locals() if locals()[name] is obj][0])
#%%