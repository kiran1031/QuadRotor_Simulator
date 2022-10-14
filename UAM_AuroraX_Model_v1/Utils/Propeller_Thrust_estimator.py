# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 16:19:49 2022

@author: I0001386

This script is dedicated to estimate thrust from propeller 

"""


## ************************************************************************* ##
##                          Import statements                                ##
## ************************************************************************* ## 
from numpy import pi, sqrt, arctan2, array, cross
from scipy.optimize import fsolve


## ************************************************************************* ##
##           Thrust using Blade element and momentum theory                  ##
## ************************************************************************* ## 
class PropellerThrust_BET_MT_Handler1(object):
    
    '''
    
    This class is a model for propeller using blade element and moment theory
    Reference: https://www.hindawi.com/journals/ijae/2018/9632942/
    
    name:  name of the propeller which identifies its location on the UAM
    pitch: pitch of the propeller
    dia:   diameter of the propeller
    Nb:    Number of propellers
    units: Whatever unit in which pitch and dia are mentioned
    
    '''
    
    def __init__(self, name:  str='',
                       pitch: float=0.0,
                       dia:   float=0.0,
                       Nb:    int=2,
                       units: str='meters'):
        self.__name                = name
        self.__pitch               = pitch
        self.__dia                 = dia
        self.__Radius              = dia/2.0
        self.__Nb                  = Nb
        self.__units               = units
        
        
##---------------------------------------------------------------------------##
    def __RPM_to_rps(self, rpm: float=0.0):
        '''
        
        This method converts rotation per minute RPM to 
        RPS (radians per second)
        
        '''

        return 2*pi*rpm/60.
    
    
##---------------------------------------------------------------------------##
    def __unitsHandle(self, parameter: float=0.0, convertTo: str='meters'):
        '''
        
        This method converts units of distance between one and another
        parameter: parameter to convert units for
        convertTo: string to indicate to which units parameter should be
        converted to
        
        '''
        if self.__units == 'meters' and convertTo == 'inches':
            return 39.3701*parameter
        
        elif self.__units == 'inches' and convertTo == 'meters':
            return parameter/39.3701
        
        elif self.__units == 'inches' and convertTo == 'inches':
            return parameter
        
        elif self.__units == 'meters' and convertTo == 'meters':
            return parameter
    
    
##---------------------------------------------------------------------------##
    def __computeTheta(self):
        '''
        
        This method computes angle relating pitch and rotation
        
        '''
        
        return arctan2(pi*self.__dia, self.__pitch)
    
    
##---------------------------------------------------------------------------##
    def __get_k(self):
        '''
        
        This method computes k

        '''
        
        d                          = self.__unitsHandle(parameter=self.__dia, convertTo='inches')
        
        if d == 4.0:
            c_by_d                 = 0.09
            
        elif d>=5. and d<=6.:
            c_by_d                 = 0.1
            
        elif d>=7. and d<=9.:
            c_by_d                 = 0.11
            
        elif d>=10. and d<=12.:
            c_by_d                 = 0.12
            
        elif d>=13. and d<=14.:
            c_by_d                 = 0.13
            
        elif d>=15. and d<=16.:
            c_by_d                 = 0.14
            
        return 0.5*self.__Nb*c_by_d
    
    
##---------------------------------------------------------------------------##
    def __get_ed(self):
        '''
        
        This method estimates the effective propeller diameter

        '''
        
        p_by_d                     = self.__pitch/self.__dia
        
        if p_by_d < 0.4:
            ed                     = 0.91
        
        elif p_by_d >= 0.4 and p_by_d < 0.8:
            ed                     = 0.88
            
        elif p_by_d >=0.8 and p_by_d<0.9:
            ed                     = 0.86
            
        elif p_by_d>=0.9:
            ed                     = 0.8
            
        else:
            ed                     = None
            
        return ed
        
    
##---------------------------------------------------------------------------##
    def __getCT(self):
        '''
        
        This method computes the thrust coefficient CT
        rpm: Rotations per minute
        
        '''
        
        tht                        = self.__computeTheta()
        k                          = self.__get_k()
        ed                         = self.__get_ed()
        
        if not ed:
            raise Exception('ed is not evaluated. Check the calculation')
        
        CT                         = ((4./3.)*k*tht*(1-(1-ed)**3)) + (k*(sqrt(k*(1+k)) - sqrt(k))*(1-(1-ed)**2))
                                     
        return CT, ed
    

##---------------------------------------------------------------------------##    
    def getThrust(self, rpm: float=0.0,
                        rho: float=1.225):
        '''
        
        This method computes thrust given RPM
        rpm: Rotations per minute of the propeller
        rho: Fluid density in which the propeller is operating
        
        '''
        
        ## get the angular velocity in radians per sec
        self.__w                   = self.__RPM_to_rps(rpm)
        
        CT , ed                    = self.__getCT()
        R                          = self.__unitsHandle(parameter=self.__Radius, convertTo='meters')
        
        Thrust                     = (1./6.)*rho*pi*(R**4)*(ed**4)*CT*self.__w**2
        
        return Thrust
        
        
## ************************************************************************* ##
##           Thrust using Blade element and momentum theory 2                ##
## ************************************************************************* ##
class PropellerThrust_BET_MT_Handler2(object):
    
    '''
    
    This class is a model for propeller using blade element and moment theory
    Reference: https://charlestytler.com/modeling-vehicle-dynamics-6dof-nonlinear-simulation/
    
    name:  name of the propeller which identifies its location on the UAM
    pitch: pitch of the propeller
    dia:   diameter of the propeller
    Nb:    Number of propellers
    c:     Mean chord length
    a:     Lift curve slope
    eta:   Propeller efficiency
    dx:    x Position of motor propeller thrust axis from CG of UAM
    dy:    y position of motor propeller thrust axis from CG of UAM
    units: Whatever unit in which pitch and dia are mentioned
    
    '''
    
    def __init__(self, name:  str='',
                       pitch: float=0.0,
                       dia:   float=0.0,
                       Nb:    int=2,
                       c:     float=0.0,
                       a:     float=5.7,
                       eta:   float=1.0,
                       dx:    float=0.0,
                       dy:    float=0.0,
                       units: str='meters'):
        self.__name                = name
        self.__pitch               = pitch
        self.__dia                 = dia
        self.__Radius              = dia/2.0
        self.__Nb                  = Nb
        self.__chord               = c
        self.__liftslope           = a
        self.__eta                 = eta
        self.__dx                  = dx
        self.__dy                  = dy
        self.__units               = units
        self.__area                = self.__computeArea()
    
    
##---------------------------------------------------------------------------##
    def __computeArea(self):
        '''
        
        Returns the Area of the propeller given the radius

        '''
        
        return pi*self.__Radius**2
    
    
##---------------------------------------------------------------------------##
    def __RPM_to_rps(self, rpm: float=0.0):
        '''
        
        This method converts rotation per minute RPM to 
        RPS (radians per second)
        
        '''

        return 2*pi*rpm/60.
    
    
##---------------------------------------------------------------------------##
    def __thrustEqn(self, vi, *prop_params):
    
        '''
        
        This method computes the thrust equation given the propeller parameters
        and induced velocity
        
        '''
        
        ## Unpack parameters
        rho,theta0,theta1,U,V,W    = prop_params
        
        ## Calculate local airflow velocity at propeller with vi, V'
        Vprime                     = sqrt(U**2 + V**2 + (W - vi)**2)
        
        ## Calculate Thrust averaged over one revolution of propeller using vi
        Thrust                     = 1/4 * rho * self.__liftslope * self.__Nb * self.__chord * self.__Radius * \
                                    ( (W - vi) * self.__w * self.__Radius + 2/3 * (self.__w * self.__Radius)**2 * (theta0 + 3/4 * theta1) + \
                                    (U**2 + V**2) * (theta0 + 1/2 * theta1) )
        
        ## Calculate residual for equation: Thrust = mass flow rate * delta Velocity
        residual                   = self.__eta * 2 * vi * rho * self.__area * Vprime - Thrust
        
        return residual
    
    
##---------------------------------------------------------------------------##
    def getThrust(self, rpm: float=0.0,
                        rho: float=1.225,
                        ub:  float=0.0,
                        vb:  float=0.0,
                        wb:  float=0.0,
                        p:   float=0.0,
                        q:   float=0.0,
                        r:   float=0.0,
                        vi0: float=0.1):
        '''
        
        This method computes thrust given RPM
        rpm: Rotations per minute of the propeller
        rho: Fluid density in which the propeller is operating
        ub:  UAM linear velocity in body frame x axis
        vb:  UAM linear velocity in body frame y axis
        wb:  UAM linear velocity in body frame z axis
        p:   UAM angular velocity about body frame x axis
        q:   UAM angular velocity about body frame y axis
        r:   UAM angular velocity about body frame z axis
        vi0: Initial guess for the propeller induced velocity
        
        '''
        
        ## compute the fluid velocity vector components in body frame
        U,V,W                     = array([ub,vb,wb]) + cross(array([p,q,r]), array([self.__dx,self.__dy,0]))
        
        ## get the angular velocity in radians per sec
        self.__w                   = self.__RPM_to_rps(rpm)
        
        ## compute theta0 and theta1
        theta0                     = 2*arctan2(self.__pitch, (2 * pi * 3/4 * self.__dia/2))
        theta1                     = -4 / 3 * arctan2(self.__pitch, 2 * pi * 3/4 * self.__dia/2)
        
        ## pack the propeller params
        prop_params                = (rho,theta0,theta1,U,V,W)
        
        ## Numerically solve for propeller induced velocity, vi
        ## using nonlinear root finder, fsolve, and prop_params
        vi                         = fsolve(self.__thrustEqn, vi0, args=prop_params)
        
        ## obtain the thrust from the converged induced velocity
        Vprime                     = sqrt(U**2 + V**2 + (W - vi)**2)
        Thrust                     = self.__eta * 2 * vi * rho * self.__area * Vprime
        
        return Thrust
        
        
    





