# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 19:02:49 2022

@author: I0001386

This script handles creation and all operations on ReferenceFrame model

"""


## ************************************************************************* ##
##                          Import statements                                ##
## ************************************************************************* ## 
from sympy.physics.mechanics import ReferenceFrame


## ************************************************************************* ##
##                           ReferenceFrame Handler                          ##
## ************************************************************************* ##
class ReferenceFrameModel(object):
    
    '''
    
    This class creates ReferenceFrame object from sympy mechanics
    
    '''
    
    def __init__(self, frame_name: str=''):
        self.__frame_name           = frame_name
        
        
    def createFrame(self):
        return ReferenceFrame(self.__frame_name)


