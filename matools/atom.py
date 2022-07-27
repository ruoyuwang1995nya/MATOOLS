# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:11:09 2022

@author: 王若宇
"""
#class element():
elements={'H':(6.941,),'Ag':(107.9,),'In':(114.8,),'Te':(127.6,)}

class atom():
    '''
    atom properties, e.g., element, atomic mass, etc.
    '''
    def __init__(self,element):
        self._element=element
        #atom._mass=elements[self.element][0]
        
    @property
    def element(self):
        return self._element
    @element.setter
    def element(self,name:str):
        #if ele not in elements.keys():
        #    raise TypeError("Non-existent element!")
        #else:
            #self._element=ele
        self._element=name
        
    #@property
    #def mass(self):
    #    return self._mass
    #@mass.setter
    #def mass(self):
        
    
    