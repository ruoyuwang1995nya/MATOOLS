# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 18:02:37 2022

@author: 王若宇
"""

import numpy as np

# ==================== vector =============================
class vector_base():
    '''
    base class for vector
    dimension: 0: any dimension
    '''
    def __init__(self,vec,d=0):
        self._d=d
        self.vec=vec
        
    def __add__(self,other):
        if self.d==other.d:
            vec=self.vec+other.vec
            return vector_base(vec,d=self.d)
        else:
            print('The two vector must have the same dimension')
            raise RuntimeError
            
    def __sub__(self,other):
        if self.d==other.d:
            vec=self.vec-other.vec
            return vector_base(vec,d=self.d)
        else:
            print('The two vector must have the same dimension')
            raise RuntimeError
        
        
    @property
    def d(self):
        return self._d
    @d.setter
    def d(self,dimension):
        if isinstance(dimension,int):
            if dimension < 0:
                raise ValueError("Dimension cannot be negative")
            else:
                self._d=dimension
        else:
            raise TypeError("Dimension must be zero or a positive integer")
        
    @property
    def vec(self):
        return self._vec
    @vec.setter
    def vec(self,vec_a=np.array([1.,0.,0.])):
        if isinstance(vec_a,np.ndarray) or isinstance(vec_a,list) or isinstance(vec_a,tuple):  
            try:
                vec_a=np.array(vec_a).astype('float64')
            except ValueError:
                print('check input format')
                raise
            if self.d!=0:
                if vec_a.shape!=(self.d,):
                    raise ValueError("Wrong dimension!")
                else:
                    self._vec=vec_a
            else:
                self._vec=vec_a
        else:
            raise TypeError("Must be list, tuple or numpy array")
            
    
            
    def modulus(self):
        return sum(self.vec**2)**0.5
    

        
class lat_vec(vector_base):
    '''
    sub-class for lattvec
    '''
    def __init__(self,vec):
        super().__init__(vec,d=3)

    @property
    def vec(self):
        return super().vec
    @vec.setter
    def vec(self,new_vec):
        _value = new_vec
        super(type(self),lat_vec).vec.fset(self, _value)


def multiply(vec1:vector_base,vec2:vector_base):
    '''Multiply two vectors'''
    if vec1.d==vec2.d:    
        return sum([vec1.vec[i]*vec2.vec[i] for i in range(vec1.d)])
        
        
