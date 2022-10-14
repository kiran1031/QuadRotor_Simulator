# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 20:24:57 2022

@author: I0001386

This script has predefined controllers for different manuevers
Manuevers available are
1. Start
2. Climb
3. Hover_At
4. Translate_(x,y,z)
5. Land

"""


## ************************************************************************* ##
##                          Import statements                                ##
## ************************************************************************* ## 
import warnings


## ************************************************************************* ##
##                     Class with all control methods                        ##
## ************************************************************************* ## 
class ControlHandler(object):
    
    '''
    
    1. start     --> thrust just balances the weight of the UAM
    2. climb     --> climbs with fixed climb thrust
    3. hover_at  --> hovers at a fixed altitude
    4. translate --> translates in x or y direction
    5. land      --> lands with fixed landing thrust
    
    '''

    def __init__(self, params:dict={}):
        self.__params                      = params
        self.manuevers                     = []
        self.curr_t                        = 0.
        
        
##---------------------------------------------------------------------------##
    def start(self, tf:float=0.0):
        if tf < self.curr_t or tf > self.__params['tf']:
            warnings.warn('Cannot execute start manuever as final time requested is invalid')
        self.manuevers.append(('start', self.curr_t, tf)) 
        self.curr_t                         = tf
        #self.curr_t                        = self.curr_t + tf
        
        
##---------------------------------------------------------------------------##
    def climb(self, tf:float=0.0):
        if tf < self.curr_t or tf > self.__params['tf']:
            warnings.warn('Cannot execute climbing manuever as final time requested is invalid')
        self.manuevers.append(('climb', self.curr_t, tf))
        self.curr_t                         = tf
        #self.curr_t                        = self.curr_t + tf
        
        
##---------------------------------------------------------------------------##
    def hover(self, tf:float=0.0, altitude:float=0.0):
        if tf < self.curr_t or tf > self.__params['tf']:
            warnings.warn('Cannot execute hovering manuever as final time requested is invalid')
        self.hover_alt                      = altitude
        self.manuevers.append(('hover', self.curr_t, tf))
        self.curr_t                         = tf
        #self.curr_t                        = self.curr_t + tf
        
        
##---------------------------------------------------------------------------##
    def climb_and_hover(self, tf:float=0.0, altitude:float=0.0):
        if tf < self.curr_t or tf > self.__params['tf']:
            warnings.warn('Cannot execute climb and hovering manuever as final time requested is invalid')    
        self.hover_alt                       = altitude
        self.manuevers.append(('climb_and_hover', self.curr_t, tf))
        self.curr_t                         = tf

        
##---------------------------------------------------------------------------##
    def descend(self):
        self.manuevers.append(('descend', self.curr_t, self.__params['tf']))
        #self.curr_t                        = self.curr_t + tf
        
        
##---------------------------------------------------------------------------##
    def update_control(self, t, X, params):
        control_identified                 = ''
        for elem in self.manuevers:
            if t >= elem[1] and t < elem[2]:
                control_identified         = elem[0]
                break
            
        if not control_identified:
            raise Exception('control failed')
            
        if control_identified == 'start':
            params['Ft/mg']                = 1.0
            params['L']                    = 0.0
            params['M']                    = 0.0
            params['N']                    = 0.0
            
        if control_identified == 'climb':
            params['Ft/mg']                = 1.05
            params['L']                    = 0.0
            params['M']                    = 0.0
            params['N']                    = 0.0
            
        if control_identified == 'hover':
            params['Ft/mg']                = 1.0
            params['L']                    = 0.0
            params['M']                    = 0.0
            params['N']                    = 0.0
            
        if control_identified == 'climb_and_hover':
            params['L']                    = 0.0
            params['M']                    = 0.0
            params['N']                    = 0.0
            
        if control_identified == 'descend':
            params['Ft/mg']                = 0.9
            params['L']                    = 0.0
            params['M']                    = 0.0
            params['N']                    = 0.0
        
        return params, control_identified
    
        
##---------------------------------------------------------------------------##
    def showManuevers(self):
        print('*******************Manuevers Executed*************************')
        for idx, elem in enumerate(self.manuevers):
            print(idx, '    -->    ', elem[0], ' executed from ', elem[1], ' till ', elem[2], ' seconds')

































