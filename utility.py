# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 17:03:24 2022

@author: 王若宇
"""

import numpy as np

def iterative_print (ite_obj):
    '''print iteratives as string'''
    output=''
    for i in ite_obj:
        if type(i)==np.float64:
            output += (r' {:.16f}'.format(i))
        else:
            output += (r' {}'.format(i))
    return output