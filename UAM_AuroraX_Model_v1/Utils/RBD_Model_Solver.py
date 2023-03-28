# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:47:17 2022

@author: i0001386

This script solves or integrates the equations of motion 

https://www.wilselby.com/research/arducopter/controller-design/

"""


## ************************************************************************* ##
##                          Import statements                                ##
## ************************************************************************* ## 
from numpy import arange, zeros, size, pi, sqrt, arctan2, cos, sin, array, where, radians, empty, append, transpose
from numpy.linalg import inv
from matplotlib.pyplot import plot, xlabel, ylabel, grid, axhline
from warnings import warn
import sys
from PyQt5 import QtCore
from PyQt5.QtWidgets import QHBoxLayout, QGroupBox, QDialog, QVBoxLayout, QLCDNumber, QLabel
from OCC.Display.backend import load_backend
load_backend('qt-pyqt5')
from OCC.Display.qtDisplay import qtViewer3d
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Pnt, gp_Dir, gp_Trsf, gp_Vec
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Display.SimpleGui import init_display


## ************************************************************************* ##
##                         RPM to rps converter                              ##
## ************************************************************************* ## 
def RPM_to_RPS(rpm: float=0.0):
    return 2*pi*rpm/60.


## ************************************************************************* ##
##                         RPS to rpm converter                              ##
## ************************************************************************* ## 
def RPS_to_RPM(rps: float=0.0):
    return (rps*60.)/(2*pi)


## ************************************************************************* ##
##                          RK4 STEP Algorithm                               ##
## ************************************************************************* ## 
def RK4_step1(func, y0, t0, dt, params):
    k1   = dt * func(t0, y0, params)
    k2   = dt * func(t0 + 0.5*dt, y0 + 0.5*k1, params)
    k3   = dt * func(t0 + 0.5*dt, y0 + 0.5*k2, params)
    k4   = dt * func(t0 + dt, y0 + k3, params)
    
    ynew = y0 + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
    
    return ynew


def pid_force1(t, err, derr_dt, int_err, kp, kd, ki):
    return (kp * err) + (kd * derr_dt) + (ki * int_err)


## ************************************************************************* ##
##       RBD Solver with Controller (Linear equations in inertial frame)     ##
## ************************************************************************* ##
class RBDSolverHandler1(object):
    
    '''
    
    This class integrates equations of motion with PID control
    for manuevers specified.
    
    The linear equations of motion are in the inertial frame
    The angular equations of motion are also in the body frame
    
    INPUTS :  w1, w2, w3, w4
    OUTPUTS:  ue, ve, we, p, q, r, phi, tht, psi, xe, ye, ze
    
    '''

    def __init__(self, rbd_model,
                       params: dict={}):
        self.model                         = rbd_model
        self.params                        = params
        self.__getRPMtoLoadsMatrix()
        self.__getTimeParams()
        self.__initializePlottingDict()
        self.n_states                      = 12                # Number of states
        self.n_inputs                      = 4                 # Number of inputs
        self.x                             = zeros((size(self.params['t']), self.n_states))      # time history of state vectors
        self.__hover_alt_tol               = 0.1                               # 10 cm
        self.inp                           = zeros((size(self.params['t']), self.n_inputs))      # time history of input vectors
        self.plotNo                        = 1


##---------------------------------------------------------------------------##
    def __getTimeParams(self):
        self.params['t']                   = arange(0, self.params['tf'], self.params['dt'])
        
        
##---------------------------------------------------------------------------##
    def __initializePlottingDict(self):
        self.__plotDict                    = {
            'u'      :          (0,  'x velocity(m/s)'),
            'v'      :          (1,  'y velocity(m/s)'),
            'w'      :          (2,  'z velocity(m/s)'),
            'p'      :          (3,  'roll rate(rad/s)'),
            'q'      :          (4,  'pitch rate(rad/s)'),
            'r'      :          (5,  'yaw rate(rad/s)'),
            'phi'    :          (6,  'Roll angle(deg)'),
            'tht'    :          (7,  'pitch Angle(deg)'),
            'psi'    :          (8,  'yaw Angle(deg)'),
            'xe'     :          (9,  'x (m)'),
            'ye'     :          (10, 'y (m'),
            'ze'     :          (11, 'Altitude (m'),
            'rpm'    :          (12, 'RPM')
            }
        
        
##---------------------------------------------------------------------------##
    def __getRPMtoLoadsMatrix(self):
        km                                 = self.params['km']
        bm                                 = self.params['bm']
        self.w_to_Loads_matrix             = array([
            [km, km, km, km],
            [km*self.params['d1y'], km*self.params['d2y'], km*self.params['d3y'], km*self.params['d4y']],
            [-km*self.params['d1x'], -km*self.params['d2x'], -km*self.params['d3x'], -km*self.params['d4x']],
            [bm, -bm, bm, -bm]
            ])
        
        
##---------------------------------------------------------------------------##
    def getRPMgivenLoads(self, Ft, L, M, N):
        
        '''
        
        This function computes RPM for a specific Thrust, L, M, and N

        '''
        
        self.loads_to_w_matrix             = inv(self.w_to_Loads_matrix)
        w1_sq, w2_sq, w3_sq, w4_sq         = self.loads_to_w_matrix @ array([Ft, L, M, N])
        return sqrt(w1_sq)*60/(2*pi), sqrt(w2_sq)*60/(2*pi), sqrt(w3_sq)*60/(2*pi), sqrt(w4_sq)*60/(2*pi)
    
    
##---------------------------------------------------------------------------##
    def getStaticEqRPM(self):
        
        '''

        This function returns the RPM at which the Thrust force exactly 
        balances the weight of the UAM for the given configuration

        '''
        
        w1, w1, w1, w1                     = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0., 0., 0.)
        return w1
        
        
##---------------------------------------------------------------------------##
    def __eval_eq_of_mot(self, t, y, params, raw_loads=False):
        
        '''
        
        This method evaluates the non linear equations of motions given all the 
        parameters of the system
        
        '''  
        
        ## expand the params
        mv                                 = params['m']
        Ixxv                               = params['Ixx']
        Iyyv                               = params['Iyy']
        Izzv                               = params['Izz']
        gv                                 = params['g']
        
        
        if not raw_loads:
            w1                             = RPM_to_RPS(params['w1'])
            w2                             = RPM_to_RPS(params['w2'])
            w3                             = RPM_to_RPS(params['w3'])
            w4                             = RPM_to_RPS(params['w4'])
            
            Ftv, Lv, Mv, Nv                = self.w_to_Loads_matrix @ array([w1**2, w2**2, w3**2, w4**2])
        else:
            Ftv                            = params['Ft']
            Lv                             = params['L']
            Mv                             = params['M']
            Nv                             = params['N']
        #print(Ftv/(mv*gv))
        
        uev                                = y[0]
        vev                                = y[1]
        wev                                = y[2]
        pv                                 = y[3]
        qv                                 = y[4]
        rv                                 = y[5]
        phiv                               = y[6]
        thtv                               = y[7]
        psiv                               = y[8]
        #xEv                                = y[9]
        #yEv                                = y[10]
        #zEv                                = y[11]
        
        output                             = zeros(12)
        for i in range(12):
            output[i]                      = self.model.eq_of_mot_func[i](uev,vev,wev,pv,qv,rv,phiv,thtv,psiv,mv,Ixxv,Iyyv,Izzv,gv,Ftv,Lv,Mv,Nv)
        
        return output


##---------------------------------------------------------------------------##
    def getEulerAngleRate(self, euler_param:str='thtd', p:float=0., q:float=0.,
                          r:float=0., phi:float=0., tht:float=0., psi:float=0.):
        
        '''
        
        This method evaluates phid, thtd, psid given p, q and r, phi, tht, psi
        
        '''
        
        if euler_param == 'phid':
            return self.model.eq_of_mot_func[6](0.,0.,0.,p,q,r,phi,tht,psi,0.,1.,1.,1.,9.81,0.,0.,0.,0.)
        
        elif euler_param == 'thtd':
            return self.model.eq_of_mot_func[7](0.,0.,0.,p,q,r,phi,tht,psi,0.,1.,1.,1.,9.81,0.,0.,0.,0.)
        
        elif euler_param == 'psid':
            return self.model.eq_of_mot_func[8](0.,0.,0.,p,q,r,phi,tht,psi,0.,1.,1.,1.,9.81,0.,0.,0.,0.)
        
        
##---------------------------------------------------------------------------##
    def getAngularVelocityRate(self, ang_vel_param:str='pd', p:float=0., q:float=0.,
                          r:float=0., L:float=0., M:float=0., N:float=0.):
        
        '''
        
        This method evaluates pd, qd, rd from the dynamic equations of motion
        
        '''
        
        if ang_vel_param == 'pd':
            return self.model.eq_of_mot_func[3](0.,0.,0.,p,q,r,0.,0.,0.,self.params['m'],self.params['Ixx'],self.params['Iyy'],self.params['Izz'],self.params['g'],0.,L,0.,0.)
        
        elif ang_vel_param == 'qd':
            return self.model.eq_of_mot_func[4](0.,0.,0.,p,q,r,0.,0.,0.,self.params['m'],self.params['Ixx'],self.params['Iyy'],self.params['Izz'],self.params['g'],0.,0.,M,0.)
        
        elif ang_vel_param == 'rd':
            return self.model.eq_of_mot_func[5](0.,0.,0.,p,q,r,0.,0.,0.,self.params['m'],self.params['Ixx'],self.params['Iyy'],self.params['Izz'],self.params['g'],0.,0.,0.,N)

    
##---------------------------------------------------------------------------##
    def solve(self):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        self.x[0,:]                        = X0
        
        for t_idx in range(len(self.params['t'])-1):
            
            y_tmp                          = RK4_step1(self.__eval_eq_of_mot, self.x[t_idx,:], self.params['t'][t_idx], self.params['dt'], self.params)
            
            ## restricting altitude to physical solutions
            if y_tmp[11] < 0:
                break
            
            self.x[t_idx+1, :]             = y_tmp
        
        print('*******************Solver ran successfully********************')


##---------------------------------------------------------------------------##
    def hover_at(self, hover_alt: float=0.0, ti: float=0.0, tf:float=None):
        ti_idx                             = where(self.params['t'] == ti)[0]
        if ti_idx.size == 0:
            raise Exception('requested start time is invalid')
        else:
            ti_idx                         = ti_idx[0]
            
        tf_idx                             = where(self.params['t'] == tf)[0]
        if tf_idx.size == 0:
            tf_idx                         = len(self.params['t'])
        else:
            tf_idx                         = tf_idx[0]
            
        if ti_idx == 0:
            if 'X0' in list(self.params.keys()):
                X0                         = self.params['X0']
            else:
                X0                         = zeros(12)
            self.x[0,:]                    = X0
            y_tmp, y_new                   = X0.copy(), X0.copy()
        else:
            y_tmp, y_new                   = self.x[ti_idx,:].copy(), self.x[ti_idx,:].copy()
            
        if abs(y_tmp[11] - hover_alt) > self.__hover_alt_tol:
            warn('Altitude to hover is not close to current altitude within tolerance')
            return
        
        for t_idx in range(ti_idx, tf_idx):
            
            t0                 = self.params['t'][t_idx]
            self.x[t_idx, :]   = y_new
            err_h              = y_tmp[11] - hover_alt
            derr_dt_h          = y_tmp[2]
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
            #print(Ft)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
            k1                 = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params, raw_loads=False)
            
            err_h              = (y_tmp + 0.5*k1)[11] - hover_alt
            derr_dt_h          = (y_tmp + 0.5*k1)[2]
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
            #print(Ft)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
            k2                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
            
            err_h              = (y_tmp + 0.5*k2)[11] - hover_alt
            derr_dt_h          = (y_tmp + 0.5*k2)[2]
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
            #print(Ft)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
            k3                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
            
            err_h              = (y_tmp + k3)[11] - hover_alt
            derr_dt_h          = (y_tmp + k3)[2]
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
            #print(Ft)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
            k4                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
            
            y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp              = y_new.copy()
            
            self.inp[t_idx,:]  = array([self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']])
            
            
##---------------------------------------------------------------------------##
    def climb_and_hover_at(self, hover_alt:float = 0.0, climb_speed:float=0.5):
       # ws                     = self.getStaticEqRPM()
        
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        y_tmp, y_new                       = X0.copy(), X0.copy()
        self.x[0,:]                        = X0
        
        ti_hover                           = None
        
        for t_idx in range(len(self.params['t'])-1):
            
            t0                 = self.params['t'][t_idx]
            self.x[t_idx, :]   = y_new
            err_h              = y_tmp[11] - hover_alt
            derr_dt_h          = y_tmp[2] - climb_speed
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
            #print(Ft)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
            k1                 = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params, raw_loads=False)
            
            err_h              = (y_tmp + 0.5*k1)[11] - hover_alt
            derr_dt_h          = (y_tmp + 0.5*k1)[2] - climb_speed
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
            #print(Ft)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
            k2                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
            
            err_h              = (y_tmp + 0.5*k2)[11] - hover_alt
            derr_dt_h          = (y_tmp + 0.5*k2)[2] - climb_speed
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
            #print(Ft)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
            k3                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
            
            err_h              = (y_tmp + k3)[11] - hover_alt
            derr_dt_h          = (y_tmp + k3)[2] - climb_speed
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
            #print(Ft)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
            k4                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
            
            y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp              = y_new.copy()
            
            ## restricting altitude to physical solutions
            if y_tmp[11] < 0:
                break
            
            self.x[t_idx+1, :]             = y_tmp
            self.inp[t_idx,:]              = array([self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']])
            
            ## if altitude reaches close to hover altitude requested, also break and start closed loop  control
            if abs(y_tmp[11] - hover_alt) <= self.__hover_alt_tol:
                ti_hover                   = self.params['t'][t_idx+1]
                break
            
        self.hover_at(hover_alt=hover_alt, ti=ti_hover)
        
        print('*******************Solver ran successfully********************')
            
            
##---------------------------------------------------------------------------##
    def stabilize_roll_at(self, phi: float=0.0):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        y_tmp, y_new                       = X0.copy(), X0.copy()
        self.x[0,:]                        = X0
        
        phi_rad                          = radians(phi)
        
        for t_idx in range(len(self.params['t'])-1):
            
            t0                 = self.params['t'][t_idx]
            self.x[t_idx, :]   = y_new
            err_phi            = y_tmp[6] - phi_rad
            derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
            L                  = -pid_force1(t0, err_phi, derr_dt_phi, 0.0, self.params['kp_phi'], self.params['kd_phi'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], L, 0., 0.)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0., 0.
            k1                 = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params, raw_loads=False)
            
            err_phi            = (y_tmp + 0.5*k1)[6] - phi_rad
            derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
            L                  = -pid_force1(t0, err_phi, derr_dt_phi, 0.0, self.params['kp_phi'], self.params['kd_phi'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], L, 0, 0.)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k2                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
            
            err_phi            = (y_tmp + 0.5*k2)[6] - phi_rad
            derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
            L                  = -pid_force1(t0, err_phi, derr_dt_phi, 0.0, self.params['kp_phi'], self.params['kd_phi'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], L, 0., 0.)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k3                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
            
            err_phi            = (y_tmp + k3)[6] - phi_rad
            derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
            L                  = -pid_force1(t0, err_phi, derr_dt_phi, 0.0, self.params['kp_phi'], self.params['kd_phi'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], L, 0, 0.)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k4                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
            
            y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp              = y_new.copy()
            
            ## restricting altitude to physical solutions
            if y_tmp[11] < 0:
                break
            
            self.x[t_idx+1, :]             = y_tmp
            self.inp[t_idx,:]              = array([self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']])
        
        
##---------------------------------------------------------------------------##
    def stabilize_pitch_at(self, theta: float=0.0):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        y_tmp, y_new                       = X0.copy(), X0.copy()
        self.x[0,:]                        = X0
        
        theta_rad                          = radians(theta)
        
        for t_idx in range(len(self.params['t'])-1):
            
            t0                 = self.params['t'][t_idx]
            self.x[t_idx, :]   = y_new
            err_tht            = y_tmp[7] - theta_rad
            derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
            M                  = -pid_force1(t0, err_tht, derr_dt_tht, 0.0, self.params['kp_tht'], self.params['kd_tht'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0., M, 0.)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), 0, M, 0.
            k1                 = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params, raw_loads=False)
            
            err_tht            = (y_tmp + 0.5*k1)[7] - theta_rad
            derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
            M                  = -pid_force1(t0, err_tht, derr_dt_tht, 0.0, self.params['kp_tht'], self.params['kd_tht'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0., M, 0.)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), 0, M, 0.
            k2                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
            
            err_tht            = (y_tmp + 0.5*k2)[7] - theta_rad
            derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
            M                  = -pid_force1(t0, err_tht, derr_dt_tht, 0.0, self.params['kp_tht'], self.params['kd_tht'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0., M, 0.)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), 0, M, 0.
            k3                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
            
            err_tht            = (y_tmp + k3)[7] - theta_rad
            derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
            M                  = -pid_force1(t0, err_tht, derr_dt_tht, 0.0, self.params['kp_tht'], self.params['kd_tht'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0., M, 0.)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), 0, M, 0.
            k4                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
            
            y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp              = y_new.copy()
            
            ## restricting altitude to physical solutions
            if y_tmp[11] < 0:
                break
            
            self.x[t_idx+1, :]             = y_tmp
            self.inp[t_idx,:]              = array([self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']])
        
        
##---------------------------------------------------------------------------##
    def stabilize_yaw_at(self, psi: float=0.0):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        y_tmp, y_new                       = X0.copy(), X0.copy()
        self.x[0,:]                        = X0
        
        psi_rad                          = radians(psi)
        
        for t_idx in range(len(self.params['t'])-1):
            
            t0                 = self.params['t'][t_idx]
            self.x[t_idx, :]   = y_new
            err_psi            = y_tmp[8] - psi_rad
            derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
            N                  = -pid_force1(t0, err_psi, derr_dt_psi, 0.0, self.params['kp_psi'], self.params['kd_psi'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0., 0., N)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0., 0.
            k1                 = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params, raw_loads=False)
            
            err_psi            = (y_tmp + 0.5*k1)[8] - psi_rad
            derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
            N                  = -pid_force1(t0, err_psi, derr_dt_psi, 0.0, self.params['kp_psi'], self.params['kd_psi'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0, 0, N)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k2                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
            
            err_psi            = (y_tmp + 0.5*k2)[8] - psi_rad
            derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
            N                  = -pid_force1(t0, err_psi, derr_dt_psi, 0.0, self.params['kp_psi'], self.params['kd_psi'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0, 0., N)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k3                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
            
            err_psi            = (y_tmp + k3)[8] - psi_rad
            derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
            N                  = -pid_force1(t0, err_psi, derr_dt_psi, 0.0, self.params['kp_psi'], self.params['kd_psi'], 0.0)
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0, 0, N)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k4                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
            
            y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp              = y_new.copy()
            
            ## restricting altitude to physical solutions
            if y_tmp[11] < 0:
                break
            
            self.x[t_idx+1, :]             = y_tmp
            self.inp[t_idx,:]              = array([self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']])
        
        
##---------------------------------------------------------------------------##
    def stabilize_roll_pitch_yaw_alt_at(self, phi: float=0.0, theta: float=0., 
                                        psi:float=0., hover_alt: float=0.):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        y_tmp, y_new                       = X0.copy(), X0.copy()
        self.x[0,:]                        = X0
        
        phi_rad                            = radians(phi)
        theta_rad                          = radians(theta)
        psi_rad                            = radians(psi)
        
        if abs(y_tmp[11] - hover_alt) > self.__hover_alt_tol:
            warn('Altitude to hover is not close to current altitude within tolerance')
            return
        
        int_err_h              = 0.0
        int_err_phi            = 0.0
        int_err_tht            = 0.0
        int_err_psi            = 0.0
        
        for t_idx in range(len(self.params['t'])-1):
            
            t0                 = self.params['t'][t_idx]
            #self.x[t_idx, :]   = y_new
            err_phi            = y_tmp[6] - phi_rad
            derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
            int_err_phi        = int_err_phi + (self.params['dt']/self.params['tf'])*err_phi
            L                  = -pid_force1(t0, err_phi, derr_dt_phi, int_err_phi, self.params['kp_phi'], self.params['kd_phi'], self.params['ki_phi'])
            err_tht            = y_tmp[7] - theta_rad
            derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
            int_err_tht        = int_err_tht + (self.params['dt']/self.params['tf'])*err_tht
            M                  = -pid_force1(t0, err_tht, derr_dt_tht, int_err_tht, self.params['kp_tht'], self.params['kd_tht'], self.params['ki_tht'])
            err_psi            = y_tmp[8] - psi_rad
            derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
            int_err_psi        = int_err_psi + (self.params['dt']/self.params['tf'])*err_psi
            N                  = -pid_force1(t0, err_psi, derr_dt_psi, int_err_psi, self.params['kp_psi'], self.params['kd_psi'], self.params['ki_psi'])
            err_h              = y_tmp[11] - hover_alt
            derr_dt_h          = y_tmp[2]
            int_err_h          = int_err_h + (self.params['dt']/self.params['tf'])*err_h
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, int_err_h, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, L, M, N)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0., 0.
            k1                 = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params, raw_loads=False)
            
            err_phi            = (y_tmp + 0.5*k1)[6] - phi_rad
            derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
            int_err_phi        = int_err_phi + (self.params['dt']/self.params['tf'])*err_phi
            L                  = -pid_force1(t0, err_phi, derr_dt_phi, int_err_phi, self.params['kp_phi'], self.params['kd_phi'], self.params['ki_phi'])
            err_tht            = (y_tmp + 0.5*k1)[7] - theta_rad
            derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
            int_err_tht        = int_err_tht + (self.params['dt']/self.params['tf'])*err_tht
            M                  = -pid_force1(t0, err_tht, derr_dt_tht, int_err_tht, self.params['kp_tht'], self.params['kd_tht'], self.params['ki_tht'])
            err_psi            = (y_tmp + 0.5*k1)[8] - psi_rad
            derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
            int_err_psi        = int_err_psi + (self.params['dt']/self.params['tf'])*err_psi
            N                  = -pid_force1(t0, err_psi, derr_dt_psi, int_err_psi, self.params['kp_psi'], self.params['kd_psi'], self.params['ki_psi'])
            err_h              = (y_tmp + 0.5*k1)[11] - hover_alt
            derr_dt_h          = (y_tmp + 0.5*k1)[2]
            int_err_h          = int_err_h + (self.params['dt']/self.params['tf'])*err_h
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, int_err_h, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, L, M, N)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k2                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
            
            err_phi            = (y_tmp + 0.5*k2)[6] - phi_rad
            derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
            int_err_phi        = int_err_phi + (self.params['dt']/self.params['tf'])*err_phi
            L                  = -pid_force1(t0, err_phi, derr_dt_phi, int_err_phi, self.params['kp_phi'], self.params['kd_phi'], self.params['ki_phi'])
            err_tht            = (y_tmp + 0.5*k2)[7] - theta_rad
            derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
            int_err_tht        = int_err_tht + (self.params['dt']/self.params['tf'])*err_tht
            M                  = -pid_force1(t0, err_tht, derr_dt_tht, int_err_tht, self.params['kp_tht'], self.params['kd_tht'], self.params['ki_tht'])
            err_psi            = (y_tmp + 0.5*k2)[8] - psi_rad
            derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
            int_err_psi        = int_err_psi + (self.params['dt']/self.params['tf'])*err_psi
            N                  = -pid_force1(t0, err_psi, derr_dt_psi, int_err_psi, self.params['kp_psi'], self.params['kd_psi'], self.params['ki_psi'])
            err_h              = (y_tmp + 0.5*k2)[11] - hover_alt
            derr_dt_h          = (y_tmp + 0.5*k2)[2]
            int_err_h          = int_err_h + (self.params['dt']/self.params['tf'])*err_h
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, int_err_h, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, L, M, N)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k3                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
            
            err_phi            = (y_tmp + k3)[6] - phi_rad
            derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
            int_err_phi        = int_err_phi + (self.params['dt']/self.params['tf'])*err_phi
            L                  = -pid_force1(t0, err_phi, derr_dt_phi, int_err_phi, self.params['kp_phi'], self.params['kd_phi'], self.params['ki_phi'])
            err_tht            = (y_tmp + k3)[7] - theta_rad
            derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
            int_err_tht        = int_err_tht + (self.params['dt']/self.params['tf'])*err_tht
            M                  = -pid_force1(t0, err_tht, derr_dt_tht, int_err_tht, self.params['kp_tht'], self.params['kd_tht'], self.params['ki_tht'])
            err_psi            = (y_tmp + k3)[8] - psi_rad
            derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
            int_err_psi        = int_err_psi + (self.params['dt']/self.params['tf'])*err_psi
            N                  = -pid_force1(t0, err_psi, derr_dt_psi, int_err_psi, self.params['kp_psi'], self.params['kd_psi'], self.params['ki_psi'])
            err_h              = (y_tmp + k3)[11] - hover_alt
            derr_dt_h          = (y_tmp + k3)[2]
            int_err_h          = int_err_h + (self.params['dt']/self.params['tf'])*err_h
            Ft                 = (self.params['m']*self.params['g']) - pid_force1(t0, err_h, derr_dt_h, int_err_h, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
            self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, L, M, N)
            #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
            k4                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
            
            y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp              = y_new.copy()
            
            ## restricting altitude to physical solutions
            if y_tmp[11] < 0:
                break
            
            self.x[t_idx+1, :]             = y_tmp
            self.inp[t_idx,:]              = array([self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']])
        
        
##---------------------------------------------------------------------------##
    def stabilize_pqr_alt_at(self, p_des:float=0., q_des:float=0., r_des:float=0., zd:float=0.):
        
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        y_tmp, y_new                       = X0.copy(), X0.copy()
        self.x[0,:]                        = X0
        
        int_err_z                          = 0.
        int_err_p                          = 0.
        int_err_q                          = 0.
        int_err_r                          = 0.
        L                                  = 0.
        M                                  = 0.
        N                                  = 0.
        
        for t_idx in range(len(self.params['t'])-1):
            
            t0                             = self.params['t'][t_idx]
            
            ## Outer loop runs first to get desired attitude for target pos
            curr_w                         = y_tmp[2]
            curr_p                         = y_tmp[3]
            curr_q                         = y_tmp[4]
            curr_r                         = y_tmp[5]
            curr_z                         = y_tmp[11]
            
            err_z                          = curr_z - zd
            derr_dt_z                      = curr_w
            int_err_z                      = int_err_z + (self.params['dt']/self.params['tf'])*err_z
            Ft                             = (self.params['m']*self.params['g']) - pid_force1(t0, err_z, derr_dt_z, int_err_z, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
            
            err_p                          = curr_p - p_des
            derr_dt_p                      = self.getAngularVelocityRate(ang_vel_param='pd', p=curr_p, q=curr_q, r=curr_r, L=L)
            int_err_p                      = int_err_p + (self.params['dt']/self.params['tf'])*err_p
            L                              = -pid_force1(t0, err_p, derr_dt_p, int_err_p, self.params['kp_p'], self.params['kd_p'], self.params['ki_p'])
            
            err_q                          = curr_q - q_des
            derr_dt_q                      = self.getAngularVelocityRate(ang_vel_param='qd', p=curr_p, q=curr_q, r=curr_r, M=M)
            int_err_q                      = int_err_q + (self.params['dt']/self.params['tf'])*err_q
            M                              = -pid_force1(t0, err_q, derr_dt_q, int_err_q, self.params['kp_q'], self.params['kd_q'], self.params['ki_q'])
            
            err_r                          = curr_r - r_des
            derr_dt_r                      = self.getAngularVelocityRate(ang_vel_param='rd', p=curr_p, q=curr_q, r=curr_r, N=N)
            int_err_r                      = int_err_r + (self.params['dt']/self.params['tf'])*err_r
            N                              = -pid_force1(t0, err_r, derr_dt_r, int_err_r, self.params['kp_r'], self.params['kd_r'], self.params['ki_r'])
            
            self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, L, M, N
            k1                             = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params, raw_loads=True)
            
            curr_w                         = (y_tmp + 0.5*k1)[2]
            curr_p                         = (y_tmp + 0.5*k1)[3]
            curr_q                         = (y_tmp + 0.5*k1)[4]
            curr_r                         = (y_tmp + 0.5*k1)[5]
            curr_z                         = (y_tmp + 0.5*k1)[11]
            
            err_z                          = curr_z - zd
            derr_dt_z                      = curr_w
            int_err_z                      = int_err_z + (self.params['dt']/self.params['tf'])*err_z
            Ft                             = (self.params['m']*self.params['g']) - pid_force1(t0 + 0.5*self.params['dt'], err_z, derr_dt_z, int_err_z, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
            
            err_p                          = curr_p - p_des
            derr_dt_p                      = self.getAngularVelocityRate(ang_vel_param='pd', p=curr_p, q=curr_q, r=curr_r, L=L)
            int_err_p                      = int_err_p + (self.params['dt']/self.params['tf'])*err_p
            L                              = -pid_force1(t0 + 0.5*self.params['dt'], err_p, derr_dt_p, int_err_p, self.params['kp_p'], self.params['kd_p'], self.params['ki_p'])
            
            err_q                          = curr_q - q_des
            derr_dt_q                      = self.getAngularVelocityRate(ang_vel_param='qd', p=curr_p, q=curr_q, r=curr_r, M=M)
            int_err_q                      = int_err_q + (self.params['dt']/self.params['tf'])*err_q
            M                              = -pid_force1(t0 + 0.5*self.params['dt'], err_q, derr_dt_q, int_err_q, self.params['kp_q'], self.params['kd_q'], self.params['ki_q'])
            
            err_r                          = curr_r - r_des
            derr_dt_r                      = self.getAngularVelocityRate(ang_vel_param='rd', p=curr_p, q=curr_q, r=curr_r, N=N)
            int_err_r                      = int_err_r + (self.params['dt']/self.params['tf'])*err_r
            N                              = -pid_force1(t0 + 0.5*self.params['dt'], err_r, derr_dt_r, int_err_r, self.params['kp_r'], self.params['kd_r'], self.params['ki_r'])
            
            self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, L, M, N
            k2                             = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=True)
            
            curr_w                         = (y_tmp + 0.5*k2)[2]
            curr_p                         = (y_tmp + 0.5*k2)[3]
            curr_q                         = (y_tmp + 0.5*k2)[4]
            curr_r                         = (y_tmp + 0.5*k2)[5]
            curr_z                         = (y_tmp + 0.5*k2)[11]
            
            err_z                          = curr_z - zd
            derr_dt_z                      = curr_w
            int_err_z                      = int_err_z + (self.params['dt']/self.params['tf'])*err_z
            Ft                             = (self.params['m']*self.params['g']) - pid_force1(t0 + 0.5*self.params['dt'], err_z, derr_dt_z, int_err_z, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
            
            err_p                          = curr_p - p_des
            derr_dt_p                      = self.getAngularVelocityRate(ang_vel_param='pd', p=curr_p, q=curr_q, r=curr_r, L=L)
            int_err_p                      = int_err_p + (self.params['dt']/self.params['tf'])*err_p
            L                              = -pid_force1(t0 + 0.5*self.params['dt'], err_p, derr_dt_p, int_err_p, self.params['kp_p'], self.params['kd_p'], self.params['ki_p'])
            
            err_q                          = curr_q - q_des
            derr_dt_q                      = self.getAngularVelocityRate(ang_vel_param='qd', p=curr_p, q=curr_q, r=curr_r, M=M)
            int_err_q                      = int_err_q + (self.params['dt']/self.params['tf'])*err_q
            M                              = -pid_force1(t0 + 0.5*self.params['dt'], err_q, derr_dt_q, int_err_q, self.params['kp_q'], self.params['kd_q'], self.params['ki_q'])
            
            err_r                          = curr_r - r_des
            derr_dt_r                      = self.getAngularVelocityRate(ang_vel_param='rd', p=curr_p, q=curr_q, r=curr_r, N=N)
            int_err_r                      = int_err_r + (self.params['dt']/self.params['tf'])*err_r
            N                              = -pid_force1(t0 + 0.5*self.params['dt'], err_r, derr_dt_r, int_err_r, self.params['kp_r'], self.params['kd_r'], self.params['ki_r'])
            
            self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, L, M, N
            k3                             = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=True)
            
            curr_w                         = (y_tmp + k3)[2]
            curr_p                         = (y_tmp + k3)[3]
            curr_q                         = (y_tmp + k3)[4]
            curr_r                         = (y_tmp + k3)[5]
            curr_z                         = (y_tmp + k3)[11]
    
            err_z                          = curr_z - zd
            derr_dt_z                      = curr_w
            int_err_z                      = int_err_z + (self.params['dt']/self.params['tf'])*err_z
            Ft                             = (self.params['m']*self.params['g']) - pid_force1(t0 + 0.5*self.params['dt'], err_z, derr_dt_z, int_err_z, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
            
            err_p                          = curr_p - p_des
            derr_dt_p                      = self.getAngularVelocityRate(ang_vel_param='pd', p=curr_p, q=curr_q, r=curr_r, L=L)
            int_err_p                      = int_err_p + (self.params['dt']/self.params['tf'])*err_p
            L                              = -pid_force1(t0 + 0.5*self.params['dt'], err_p, derr_dt_p, int_err_p, self.params['kp_p'], self.params['kd_p'], self.params['ki_p'])
            
            err_q                          = curr_q - q_des
            derr_dt_q                      = self.getAngularVelocityRate(ang_vel_param='qd', p=curr_p, q=curr_q, r=curr_r, M=M)
            int_err_q                      = int_err_q + (self.params['dt']/self.params['tf'])*err_q
            M                              = -pid_force1(t0 + 0.5*self.params['dt'], err_q, derr_dt_q, int_err_q, self.params['kp_q'], self.params['kd_q'], self.params['ki_q'])
            
            err_r                          = curr_r - r_des
            derr_dt_r                      = self.getAngularVelocityRate(ang_vel_param='rd', p=curr_p, q=curr_q, r=curr_r, N=N)
            int_err_r                      = int_err_r + (self.params['dt']/self.params['tf'])*err_r
            N                              = -pid_force1(t0 + 0.5*self.params['dt'], err_r, derr_dt_r, int_err_r, self.params['kp_r'], self.params['kd_r'], self.params['ki_r'])
            
            self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, L, M, N
            k4                             = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + k3, self.params, raw_loads=True)
            
            y_new                          = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp                          = y_new.copy()
            
            ## restricting altitude to physical solutions
            # if y_tmp[11] < 0:
            #     break
            
            self.x[t_idx+1, :]             = y_tmp
            
            
##---------------------------------------------------------------------------##
    def move_to_XYZ(self, xd:float=0., yd:float=0., zd:float=0., psi_des:float=0., conv_tol:float=0.01):
        
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        y_tmp, y_new                       = X0.copy(), X0.copy()
        self.x[0,:]                        = X0
        
        int_err_x                          = 0.
        int_err_y                          = 0.
        int_err_z                          = 0.
        int_err_phi                        = 0.
        int_err_tht                        = 0.
        int_err_psi                        = 0.
        psi_des                            = psi_des * pi/180
        
        for t_idx in range(len(self.params['t'])-1):
            
            t0                             = self.params['t'][t_idx]
            
            ## Outer loop runs first to get desired attitude for target pos
            curr_u                         = y_tmp[0]
            curr_v                         = y_tmp[1]
            curr_w                         = y_tmp[2]
            curr_p                         = y_tmp[3]
            curr_q                         = y_tmp[4]
            curr_r                         = y_tmp[5]
            curr_x                         = y_tmp[9]
            curr_y                         = y_tmp[10]
            curr_z                         = y_tmp[11]
            curr_phi                       = y_tmp[6]
            curr_tht                       = y_tmp[7]
            curr_psi                       = y_tmp[8]
            
            err_x                          = curr_x - xd
            derr_x_dt                      = curr_u
            int_err_x                      = int_err_x + (self.params['dt']/self.params['tf'])*err_x
            theta_des                      = -pid_force1(t0, err_x, derr_x_dt, int_err_x, self.params['kp_x'], self.params['kd_x'], self.params['ki_x'])
            
            err_y                          = curr_y - yd
            derr_y_dt                      = curr_v
            int_err_y                      = int_err_y + (self.params['dt']/self.params['tf'])*err_y
            phi_des                        = -pid_force1(t0, err_y, derr_y_dt, int_err_y, self.params['kp_y'], self.params['kd_y'], self.params['ki_y'])
        
            err_z                          = curr_z - zd
            derr_dt_z                      = curr_w
            int_err_z                      = int_err_z + (self.params['dt']/self.params['tf'])*err_z
            Ft                             = (self.params['m']*self.params['g']) - pid_force1(t0, err_z, derr_dt_z, int_err_z, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
    
            err_phi                        = curr_phi - phi_des
            derr_dt_phi                    = self.getEulerAngleRate(euler_param='phid', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_phi                    = int_err_phi + (self.params['dt']/self.params['tf'])*err_phi
            L                              = -pid_force1(t0, err_phi, derr_dt_phi, int_err_phi, self.params['kp_phi'], self.params['kd_phi'], self.params['ki_phi'])
            
            err_tht                        = curr_tht - theta_des
            derr_dt_tht                    = self.getEulerAngleRate(euler_param='thtd', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_tht                    = int_err_tht + (self.params['dt']/self.params['tf'])*err_tht
            M                              = -pid_force1(t0, err_tht, derr_dt_tht, int_err_tht, self.params['kp_tht'], self.params['kd_tht'], self.params['ki_tht'])
            
            err_psi                        = curr_psi - psi_des
            derr_dt_psi                    = self.getEulerAngleRate(euler_param='psid', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_psi                    = int_err_psi + (self.params['dt']/self.params['tf'])*err_psi
            N                              = -pid_force1(t0, err_psi, derr_dt_psi, int_err_psi, self.params['kp_psi'], self.params['kd_psi'], self.params['ki_psi'])
            
            self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, L, M, N
            k1                             = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params, raw_loads=True)
            
            curr_u                         = (y_tmp + 0.5*k1)[0]
            curr_v                         = (y_tmp + 0.5*k1)[1]
            curr_w                         = (y_tmp + 0.5*k1)[2]
            curr_p                         = (y_tmp + 0.5*k1)[3]
            curr_q                         = (y_tmp + 0.5*k1)[4]
            curr_r                         = (y_tmp + 0.5*k1)[5]
            curr_x                         = (y_tmp + 0.5*k1)[9]
            curr_y                         = (y_tmp + 0.5*k1)[10]
            curr_z                         = (y_tmp + 0.5*k1)[11]
            curr_phi                       = (y_tmp + 0.5*k1)[6]
            curr_tht                       = (y_tmp + 0.5*k1)[7]
            curr_psi                       = (y_tmp + 0.5*k1)[8]
            
            err_x                          = curr_x - xd
            derr_x_dt                      = curr_u
            int_err_x                      = int_err_x + (self.params['dt']/self.params['tf'])*err_x
            theta_des                      = -pid_force1(t0 + 0.5*self.params['dt'], err_x, derr_x_dt, int_err_x, self.params['kp_x'], self.params['kd_x'], self.params['ki_x'])
            
            err_y                          = curr_y - yd
            derr_y_dt                      = curr_v
            int_err_y                      = int_err_y + (self.params['dt']/self.params['tf'])*err_y
            phi_des                        = -pid_force1(t0 + 0.5*self.params['dt'], err_y, derr_y_dt, int_err_y, self.params['kp_y'], self.params['kd_y'], self.params['ki_y'])
        
            err_z                          = curr_z - zd
            derr_dt_z                      = curr_w
            int_err_z                      = int_err_z + (self.params['dt']/self.params['tf'])*err_z
            Ft                             = (self.params['m']*self.params['g']) - pid_force1(t0 + 0.5*self.params['dt'], err_z, derr_dt_z, int_err_z, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
    
            err_phi                        = curr_phi - phi_des
            derr_dt_phi                    = self.getEulerAngleRate(euler_param='phid', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_phi                    = int_err_phi + (self.params['dt']/self.params['tf'])*err_phi
            L                              = -pid_force1(t0 + 0.5*self.params['dt'], err_phi, derr_dt_phi, int_err_phi, self.params['kp_phi'], self.params['kd_phi'], self.params['ki_phi'])
            
            err_tht                        = curr_tht - theta_des
            derr_dt_tht                    = self.getEulerAngleRate(euler_param='thtd', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_tht                    = int_err_tht + (self.params['dt']/self.params['tf'])*err_tht
            M                              = -pid_force1(t0 + 0.5*self.params['dt'], err_tht, derr_dt_tht, int_err_tht, self.params['kp_tht'], self.params['kd_tht'], self.params['ki_tht'])
            
            err_psi                        = curr_psi - psi_des
            derr_dt_psi                    = self.getEulerAngleRate(euler_param='psid', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_psi                    = int_err_psi + (self.params['dt']/self.params['tf'])*err_psi
            N                              = -pid_force1(t0 + 0.5*self.params['dt'], err_psi, derr_dt_psi, int_err_psi, self.params['kp_psi'], self.params['kd_psi'], self.params['ki_psi'])
            
            self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, L, M, N
            k2                             = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=True)
            
            curr_u                         = (y_tmp + 0.5*k2)[0]
            curr_v                         = (y_tmp + 0.5*k2)[1]
            curr_w                         = (y_tmp + 0.5*k2)[2]
            curr_p                         = (y_tmp + 0.5*k2)[3]
            curr_q                         = (y_tmp + 0.5*k2)[4]
            curr_r                         = (y_tmp + 0.5*k2)[5]
            curr_x                         = (y_tmp + 0.5*k2)[9]
            curr_y                         = (y_tmp + 0.5*k2)[10]
            curr_z                         = (y_tmp + 0.5*k2)[11]
            curr_phi                       = (y_tmp + 0.5*k2)[6]
            curr_tht                       = (y_tmp + 0.5*k2)[7]
            curr_psi                       = (y_tmp + 0.5*k2)[8]
            
            err_x                          = curr_x - xd
            derr_x_dt                      = curr_u
            int_err_x                      = int_err_x + (self.params['dt']/self.params['tf'])*err_x
            theta_des                      = -pid_force1(t0 + 0.5*self.params['dt'], err_x, derr_x_dt, int_err_x, self.params['kp_x'], self.params['kd_x'], self.params['ki_x'])
            
            err_y                          = curr_y - yd
            derr_y_dt                      = curr_v
            int_err_y                      = int_err_y + (self.params['dt']/self.params['tf'])*err_y
            phi_des                        = -pid_force1(t0 + 0.5*self.params['dt'], err_y, derr_y_dt, int_err_y, self.params['kp_y'], self.params['kd_y'], self.params['ki_y'])
        
            err_z                          = curr_z - zd
            derr_dt_z                      = curr_w
            int_err_z                      = int_err_z + (self.params['dt']/self.params['tf'])*err_z
            Ft                             = (self.params['m']*self.params['g']) - pid_force1(t0 + 0.5*self.params['dt'], err_z, derr_dt_z, int_err_z, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
    
            err_phi                        = curr_phi - phi_des
            derr_dt_phi                    = self.getEulerAngleRate(euler_param='phid', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_phi                    = int_err_phi + (self.params['dt']/self.params['tf'])*err_phi
            L                              = -pid_force1(t0 + 0.5*self.params['dt'], err_phi, derr_dt_phi, int_err_phi, self.params['kp_phi'], self.params['kd_phi'], self.params['ki_phi'])
            
            err_tht                        = curr_tht - theta_des
            derr_dt_tht                    = self.getEulerAngleRate(euler_param='thtd', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_tht                    = int_err_tht + (self.params['dt']/self.params['tf'])*err_tht
            M                              = -pid_force1(t0 + 0.5*self.params['dt'], err_tht, derr_dt_tht, int_err_tht, self.params['kp_tht'], self.params['kd_tht'], self.params['ki_tht'])
            
            err_psi                        = curr_psi - psi_des
            derr_dt_psi                    = self.getEulerAngleRate(euler_param='psid', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_psi                    = int_err_psi + (self.params['dt']/self.params['tf'])*err_psi
            N                              = -pid_force1(t0 + 0.5*self.params['dt'], err_psi, derr_dt_psi, int_err_psi, self.params['kp_psi'], self.params['kd_psi'], self.params['ki_psi'])
            
            self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, L, M, N
            k3                             = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=True)
            
            curr_u                         = (y_tmp + k3)[0]
            curr_v                         = (y_tmp + k3)[1]
            curr_w                         = (y_tmp + k3)[2]
            curr_p                         = (y_tmp + k3)[3]
            curr_q                         = (y_tmp + k3)[4]
            curr_r                         = (y_tmp + k3)[5]
            curr_x                         = (y_tmp + k3)[9]
            curr_y                         = (y_tmp + k3)[10]
            curr_z                         = (y_tmp + k3)[11]
            curr_phi                       = (y_tmp + k3)[6]
            curr_tht                       = (y_tmp + k3)[7]
            curr_psi                       = (y_tmp + k3)[8]
            
            err_x                          = curr_x - xd
            derr_x_dt                      = curr_u
            int_err_x                      = int_err_x + (self.params['dt']/self.params['tf'])*err_x
            theta_des                      = -pid_force1(t0 + 0.5*self.params['dt'], err_x, derr_x_dt, int_err_x, self.params['kp_x'], self.params['kd_x'], self.params['ki_x'])
            
            err_y                          = curr_y - yd
            derr_y_dt                      = curr_v
            int_err_y                      = int_err_y + (self.params['dt']/self.params['tf'])*err_y
            phi_des                        = -pid_force1(t0 + 0.5*self.params['dt'], err_y, derr_y_dt, int_err_y, self.params['kp_y'], self.params['kd_y'], self.params['ki_y'])
        
            err_z                          = curr_z - zd
            derr_dt_z                      = curr_w
            int_err_z                      = int_err_z + (self.params['dt']/self.params['tf'])*err_z
            Ft                             = (self.params['m']*self.params['g']) - pid_force1(t0 + 0.5*self.params['dt'], err_z, derr_dt_z, int_err_z, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
    
            err_phi                        = curr_phi - phi_des
            derr_dt_phi                    = self.getEulerAngleRate(euler_param='phid', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_phi                    = int_err_phi + (self.params['dt']/self.params['tf'])*err_phi
            L                              = -pid_force1(t0 + 0.5*self.params['dt'], err_phi, derr_dt_phi, int_err_phi, self.params['kp_phi'], self.params['kd_phi'], self.params['ki_phi'])
            
            err_tht                        = curr_tht - theta_des
            derr_dt_tht                    = self.getEulerAngleRate(euler_param='thtd', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_tht                    = int_err_tht + (self.params['dt']/self.params['tf'])*err_tht
            M                              = -pid_force1(t0 + 0.5*self.params['dt'], err_tht, derr_dt_tht, int_err_tht, self.params['kp_tht'], self.params['kd_tht'], self.params['ki_tht'])
            
            err_psi                        = curr_psi - psi_des
            derr_dt_psi                    = self.getEulerAngleRate(euler_param='psid', p=curr_p, q=curr_q, r=curr_r, phi=curr_phi, tht=curr_tht, psi=curr_psi)
            int_err_psi                    = int_err_psi + (self.params['dt']/self.params['tf'])*err_psi
            N                              = -pid_force1(t0 + 0.5*self.params['dt'], err_psi, derr_dt_psi, int_err_psi, self.params['kp_psi'], self.params['kd_psi'], self.params['ki_psi'])
            
            self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, L, M, N
            k4                             = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + k3, self.params, raw_loads=True)
            
            y_new                          = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp                          = y_new.copy()
            
            ## restricting altitude to physical solutions
            # if y_tmp[11] < 0:
            #     break
            
            self.x[t_idx+1, :]             = y_tmp
            
            x_per                          = abs(xd - y_tmp[9])/xd * 100
            y_per                          = abs(yd - y_tmp[10])/yd * 100
            z_per                          = abs(zd - y_tmp[11])/zd * 100
            
            #print('time = ', self.params['t'][t_idx], ' x_error = ', x_per, ' y_error = ', y_per, ' z_error = ', z_per)
    
            if (x_per < conv_tol) and (y_per < conv_tol) and (z_per < conv_tol):
                print('converged on position at time = ', self.params['t'][t_idx])
                for tt_idx in range(len(self.params['t'][t_idx:])-1):
                    self.x[t_idx+tt_idx+1, :]              = y_tmp
                    #print(self.params['t'][t_idx+tt_idx])
                break
            
            
##---------------------------------------------------------------------------##
    #def follow_trajectory(self, points_arr)
    
    
##---------------------------------------------------------------------------##
    def plotter_with_time(self, yvar:str='ze'):
        
        self.plotNo = self.plotNo + 1
        
        if yvar == 'phi' or yvar == 'tht' or yvar == 'psi':
            plot(self.params['t'], self.x[:, self.__plotDict[yvar][0]]*180/pi)
        elif yvar == 'rpm':
            plot(self.params['t'], self.inp[:,0], color='red')
            plot(self.params['t'], self.inp[:,1], color='blue')
            plot(self.params['t'], self.inp[:,2], color='green')
            plot(self.params['t'], self.inp[:,3], color='black')
            axhline(y=self.getStaticEqRPM(), color='red')
        else:
            plot(self.params['t'], self.x[:, self.__plotDict[yvar][0]])
        grid('on')
        xlabel('time(seconds)')
        ylabel(self.__plotDict[yvar][1])


##---------------------------------------------------------------------------##
    def __build_shape(self, grnd_size,arm_len,arm_ang,display):

        ground          = BRepPrimAPI_MakeBox(gp_Pnt(-grnd_size/2., -grnd_size/2., -0.2), grnd_size, grnd_size, 0.01).Shape()
        arm1            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-self.params['dx'], -self.params['dy'], 0.0), gp_Dir(sin(arm_ang), cos(arm_ang), 0.0)), 0.01*arm_len, arm_len, 2*pi).Shape()
        arm2            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-self.params['dx'], self.params['dy'], 0.0), gp_Dir(sin(arm_ang), -cos(arm_ang), 0.0)), 0.01*arm_len, arm_len, 2*pi).Shape()
        
        mot1            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(self.params['dx'], self.params['dy'], 0.0), gp_Dir(0.0, 0.0, 1.0)), 0.05*arm_len, 0.1*arm_len, 2*pi).Shape()
        mot2            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-self.params['dx'], -self.params['dy'], 0.0), gp_Dir(0.0, 0.0, 1.0)), 0.05*arm_len, 0.1*arm_len, 2*pi).Shape()
        mot3            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(self.params['dx'], -self.params['dy'], 0.0), gp_Dir(0.0, 0.0, 1.0)), 0.05*arm_len, 0.1*arm_len, 2*pi).Shape()
        mot4            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-self.params['dx'], self.params['dy'], 0.0), gp_Dir(0.0, 0.0, 1.0)), 0.05*arm_len, 0.1*arm_len, 2*pi).Shape()

        prop1           = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(self.params['dx']/2.0, self.params['dy'], 0.11*arm_len), gp_Dir(1.0, 0.0, 0.0)), 0.01*arm_len, self.params['dx'], 2*pi).Shape()
        prop2           = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-1.5*self.params['dx'], -self.params['dy'], 0.11*arm_len), gp_Dir(1.0, 0.0, 0.0)), 0.01*arm_len, self.params['dx'], 2*pi).Shape()
        prop3           = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(self.params['dx']/2.0, -self.params['dy'], 0.11*arm_len), gp_Dir(1.0, 0.0, 0.0)), 0.01*arm_len, self.params['dx'], 2*pi).Shape()
        prop4           = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-1.5*self.params['dx'], self.params['dy'], 0.11*arm_len), gp_Dir(1.0, 0.0, 0.0)), 0.01*arm_len, self.params['dx'], 2*pi).Shape()

        uam             = BRepAlgoAPI_Fuse(arm1, arm2).Shape()
        uam             = BRepAlgoAPI_Fuse(uam, mot1).Shape()
        uam             = BRepAlgoAPI_Fuse(uam, mot2).Shape()
        uam             = BRepAlgoAPI_Fuse(uam, mot3).Shape()
        uam             = BRepAlgoAPI_Fuse(uam, mot4).Shape()

        ground_shp      = display.DisplayColoredShape(ground, 'BLACK', update=False)[0]
        uam_shp         = display.DisplayColoredShape(uam, 'RED', update=True)[0]
        prop1_shp       = display.DisplayColoredShape(prop1, 'RED', update=True)[0]
        prop2_shp       = display.DisplayColoredShape(prop2, 'RED', update=True)[0]
        prop3_shp       = display.DisplayColoredShape(prop3, 'RED', update=True)[0]
        prop4_shp       = display.DisplayColoredShape(prop4, 'RED', update=True)[0]
        
        return ground_shp, uam_shp, prop1_shp, prop2_shp, prop3_shp, prop4_shp


##---------------------------------------------------------------------------##
    def animate(self, grnd_size:float=3.0,
                      ang_fac:float=5.0,
                      fitall:bool=True,
                      w1: float=3500,
                      w2: float=3500,
                      w3: float=3500,
                      w4: float=3500):
        arm_len             = 2*sqrt(self.params['dx']**2 + self.params['dy']**2)
        arm_ang             = arctan2(self.params['dx'], self.params['dy'])
        display, start_display, add_menu, add_function_to_menu = init_display()
        
        ground_shp, uam_shp, prop1_shp, prop2_shp, prop3_shp, prop4_shp          = self.__build_shape(grnd_size,arm_len,arm_ang,display)
        
        display.FitAll()

        prop1_trsf_spin           = gp_Trsf()
        prop2_trsf_spin           = gp_Trsf()
        prop3_trsf_spin           = gp_Trsf()
        prop4_trsf_spin           = gp_Trsf()

        prop1_trsf_trns           = gp_Trsf()
        prop2_trsf_trns           = gp_Trsf()
        prop3_trsf_trns           = gp_Trsf()
        prop4_trsf_trns           = gp_Trsf()

        prop1_trsf_rotpsi         = gp_Trsf()
        prop2_trsf_rotpsi         = gp_Trsf()
        prop3_trsf_rotpsi         = gp_Trsf()
        prop4_trsf_rotpsi         = gp_Trsf()

        prop1_trsf_rottht         = gp_Trsf()
        prop2_trsf_rottht         = gp_Trsf()
        prop3_trsf_rottht         = gp_Trsf()
        prop4_trsf_rottht         = gp_Trsf()

        prop1_trsf_rotphi         = gp_Trsf()
        prop2_trsf_rotphi         = gp_Trsf()
        prop3_trsf_rotphi         = gp_Trsf()
        prop4_trsf_rotphi         = gp_Trsf()

        uam_trsf_trns             = gp_Trsf()
        uam_trsf_rotpsi           = gp_Trsf()
        uam_trsf_rottht           = gp_Trsf()
        uam_trsf_rotphi           = gp_Trsf()

        ang1                 = 0.0
        ang2                 = 0.0
        ang3                 = 0.0
        ang4                 = 0.0
        
        for i in range(len(self.params['t'])):
            
            ## altitude animation parameter
            xpos             = self.x[i,9]
            ypos             = self.x[i,10]
            alt              = self.x[i,11]
            phival           = self.x[i,6]
            thtval           = self.x[i,7]
            psival           = self.x[i,8]
            
            #prop_rot_axis_c  = self.model.kinematics.R_B_I_func(phival, thtval, psival) @ array([0.,0.,1.])
            #prop_rot_axis_ac = self.model.kinematics.R_B_I_func(phival, thtval, psival) @ array([0.,0.,-1.])
            
            ## propeller rotation animation
            ang1             = ang1 + w1*self.params['dt']
            ang2             = ang2 + w2*self.params['dt']
            ang3             = ang3 + w3*self.params['dt']
            ang4             = ang4 + w4*self.params['dt']
            
            prop1_trsf_spin.SetRotation(gp_Ax1(gp_Pnt(self.params['dx'], self.params['dy'], 0.11*arm_len), gp_Dir(0.,0.,1.)), ang1)
            prop1_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
            prop1_trsf_rotpsi.Multiply(prop1_trsf_spin)
            prop1_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
            prop1_trsf_rottht.Multiply(prop1_trsf_rotpsi)
            prop1_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
            prop1_trsf_rotphi.Multiply(prop1_trsf_rottht)
            prop1_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
            prop1_trsf_trns.Multiply(prop1_trsf_rotphi)
            prop1Toploc    = TopLoc_Location(prop1_trsf_trns)
            display.Context.SetLocation(prop1_shp, prop1Toploc)
            
            prop2_trsf_spin.SetRotation(gp_Ax1(gp_Pnt(-self.params['dx'], -self.params['dy'], 0.11*arm_len), gp_Dir(0.,0.,1.)), ang2)
            prop2_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
            prop2_trsf_rotpsi.Multiply(prop2_trsf_spin)
            prop2_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
            prop2_trsf_rottht.Multiply(prop2_trsf_rotpsi)
            prop2_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
            prop2_trsf_rotphi.Multiply(prop2_trsf_rottht)
            prop2_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
            prop2_trsf_trns.Multiply(prop2_trsf_rotphi)
            prop2Toploc    = TopLoc_Location(prop2_trsf_trns)
            display.Context.SetLocation(prop2_shp, prop2Toploc)
            
            prop3_trsf_spin.SetRotation(gp_Ax1(gp_Pnt(self.params['dx'], -self.params['dy'], 0.11*arm_len), gp_Dir(0.,0.,-1.)), ang3)
            prop3_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
            prop3_trsf_rotpsi.Multiply(prop3_trsf_spin)
            prop3_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
            prop3_trsf_rottht.Multiply(prop3_trsf_rotpsi)
            prop3_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
            prop3_trsf_rotphi.Multiply(prop3_trsf_rottht)
            prop3_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
            prop3_trsf_trns.Multiply(prop3_trsf_rotphi)
            prop3Toploc    = TopLoc_Location(prop3_trsf_trns)
            display.Context.SetLocation(prop3_shp, prop3Toploc)
            
            prop4_trsf_spin.SetRotation(gp_Ax1(gp_Pnt(-self.params['dx'], self.params['dy'], 0.11*arm_len), gp_Dir(0.,0.,-1.)), ang4)
            prop4_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
            prop4_trsf_rotpsi.Multiply(prop4_trsf_spin)
            prop4_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
            prop4_trsf_rottht.Multiply(prop4_trsf_rotpsi)
            prop4_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
            prop4_trsf_rotphi.Multiply(prop4_trsf_rottht)
            prop4_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
            prop4_trsf_trns.Multiply(prop4_trsf_rotphi)
            prop4Toploc    = TopLoc_Location(prop4_trsf_trns)
            display.Context.SetLocation(prop4_shp, prop4Toploc)
            
            uam_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
            uam_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
            uam_trsf_rottht.Multiply(uam_trsf_rotpsi)
            uam_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
            uam_trsf_rotphi.Multiply(uam_trsf_rottht)
            uam_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
            uam_trsf_trns.Multiply(uam_trsf_rotphi)
            uamToploc     = TopLoc_Location(uam_trsf_trns)
            display.Context.SetLocation(uam_shp, uamToploc)
            
            display.Context.UpdateCurrentViewer()
            
            if fitall:
                display.FitAll()

        start_display()


## ************************************************************************* ##
##     RBD Solver with Controller/pyqt (Linear equations in inertial frame)  ##
## ************************************************************************* ##
class RBDSolverHandler2(QDialog):
    
    '''
    
    This class is a live flight controller of the quadrotor in python
    Steps to control are as follows
    
    1.  S          --> Start the quadcopter, propellers rotate at start rpm
    2.  C          --> climb, altitude of quadcopter increases, while attitude 
                         remains level
    3.  D          --> descend, altitude of quadcopter decreases, while attitude
                         remains level
    4.  uparrow    --> moves forward at pitch angle of 3 deg at current yaw
    5.  downarrow  --> moves backward at pitch angle of -3 deg at current yaw
    6.  <          --> moves left at roll angle of 3 deg at current yaw
    7.  >          --> moves right at roll angle of -3 deg at current yaw
    8.  L          --> yaw left with increment of 3 deg
    9.  R          --> yaw right with increment of -3 deg
    10. Q          --> land, by maintaining a low descend velocity such that
                         it does not crash
    11. H          --> Hover at a given altitude, while attitude remains level
    12. E          --> Exit the simulation
    
    '''

    def __init__(self, rbd_model,
                       app,
                       params: dict={}):
        super().__init__()
        self.title                         = 'UAM_AuroraX_FlightSimulator'
        self.left                          = 100
        self.top                           = 100
        self.width                         = 1800
        self.height                        = 900
        self.app                           = app
        
        self.model                         = rbd_model
        self.params                        = params
        self.__getRPMtoLoadsMatrix()
        self.__getTimeParams()
        self.__initializePlottingDict()
        self.n_states                      = 12                # Number of states
        self.n_inputs                      = 4                 # Number of inputs
        self.x                             = empty((0, self.n_states), float)      # time history of state vectors
        self.__currAlt                     = 0
        self.__hover_alt_tol               = 0.1                               # 10 cm
        self.inp                           = empty((0, self.n_inputs), float)      # time history of input vectors
        self.plotNo                        = 1
        
        self.__climbSpeed                  = 1.
        self.__descendSpeed                = -0.5
        self.__int_err_h                   = 0.0
        
        self.__S_cnt                       = 0
        self.__C_cnt                       = 0
        self.__D_cnt                       = 0
        self.__up_cnt                      = 0
        self.__down_cnt                    = 0
        self.__left_cnt                    = 0
        self.__right_cnt                   = 0
        self.__L_cnt                       = 0
        self.__R_cnt                       = 0
        self.__Q_cnt                       = 0
        self.__H_cnt                       = 0
        self.__E_cnt                       = 0
        
        self.__startSimulator()
        qtViewer3d.keyPressEvent  = self.__keyPressEvent
        
        
##---------------------------------------------------------------------------##
    def __initialize_int_err_alt(self):
        self.__int_err_h                   = 0.0
        
        
##---------------------------------------------------------------------------##
    def __initializePlottingDict(self):
        self.__plotDict                    = {
            'u'      :          (0,  'x velocity(m/s)'),
            'v'      :          (1,  'y velocity(m/s)'),
            'w'      :          (2,  'z velocity(m/s)'),
            'p'      :          (3,  'roll rate(rad/s)'),
            'q'      :          (4,  'pitch rate(rad/s)'),
            'r'      :          (5,  'yaw rate(rad/s)'),
            'phi'    :          (6,  'Roll angle(deg)'),
            'tht'    :          (7,  'pitch Angle(deg)'),
            'psi'    :          (8,  'yaw Angle(deg)'),
            'xe'     :          (9,  'x (m)'),
            'ye'     :          (10, 'y (m'),
            'ze'     :          (11, 'Altitude (m'),
            'rpm'    :          (12, 'RPM')
            }
        
        
##---------------------------------------------------------------------------##
    def __getRPMtoLoadsMatrix(self):
        km                                 = self.params['km']
        bm                                 = self.params['bm']
        self.w_to_Loads_matrix             = array([
            [km, km, km, km],
            [km*self.params['d1y'], km*self.params['d2y'], km*self.params['d3y'], km*self.params['d4y']],
            [-km*self.params['d1x'], -km*self.params['d2x'], -km*self.params['d3x'], -km*self.params['d4x']],
            [bm, -bm, bm, -bm]
            ])
        
        
##---------------------------------------------------------------------------##
    def getRPMgivenLoads(self, Ft, L, M, N):
        
        '''
        
        This function computes RPM for a specific Thrust, L, M, and N

        '''
        
        self.loads_to_w_matrix             = inv(self.w_to_Loads_matrix)
        w1_sq, w2_sq, w3_sq, w4_sq         = self.loads_to_w_matrix @ array([Ft, L, M, N])
        return sqrt(w1_sq)*60/(2*pi), sqrt(w2_sq)*60/(2*pi), sqrt(w3_sq)*60/(2*pi), sqrt(w4_sq)*60/(2*pi)
    
    
##---------------------------------------------------------------------------##
    def getStaticEqRPM(self):
        
        '''

        This function returns the RPM at which the Thrust force exactly 
        balances the weight of the UAM for the given configuration

        '''
        
        w1, w1, w1, w1                     = self.getRPMgivenLoads(self.params['m']*self.params['g'], 0., 0., 0.)
        return w1
        
        
##---------------------------------------------------------------------------##
    def __eval_eq_of_mot(self, t, y, params, raw_loads=False):
        
        '''
        
        This method evaluates the non linear equations of motions given all the 
        parameters of the system
        
        '''  
        
        ## expand the params
        mv                                 = params['m']
        Ixxv                               = params['Ixx']
        Iyyv                               = params['Iyy']
        Izzv                               = params['Izz']
        gv                                 = params['g']
        
        
        if not raw_loads:
            w1                             = RPM_to_RPS(params['w1'])
            w2                             = RPM_to_RPS(params['w2'])
            w3                             = RPM_to_RPS(params['w3'])
            w4                             = RPM_to_RPS(params['w4'])
            
            Ftv, Lv, Mv, Nv                = self.w_to_Loads_matrix @ array([w1**2, w2**2, w3**2, w4**2])
        else:
            Ftv                            = params['Ft']
            Lv                             = params['L']
            Mv                             = params['M']
            Nv                             = params['N']
        #print(Ftv/(mv*gv))
        
        uev                                = y[0]
        vev                                = y[1]
        wev                                = y[2]
        pv                                 = y[3]
        qv                                 = y[4]
        rv                                 = y[5]
        phiv                               = y[6]
        thtv                               = y[7]
        psiv                               = y[8]
        #xEv                                = y[9]
        #yEv                                = y[10]
        #zEv                                = y[11]
        
        output                             = zeros(12)
        for i in range(12):
            output[i]                      = self.model.eq_of_mot_func[i](uev,vev,wev,pv,qv,rv,phiv,thtv,psiv,mv,Ixxv,Iyyv,Izzv,gv,Ftv,Lv,Mv,Nv)
        
        return output


##---------------------------------------------------------------------------##
    def getEulerAngleRate(self, euler_param:str='thtd', p:float=0., q:float=0.,
                          r:float=0., phi:float=0., tht:float=0., psi:float=0.):
        
        '''
        
        This method evaluates phid, thtd, psid given p, q and r, phi, tht, psi
        
        '''
        
        if euler_param == 'phid':
            return self.model.eq_of_mot_func[6](0.,0.,0.,p,q,r,phi,tht,psi,0.,1.,1.,1.,9.81,0.,0.,0.,0.)
        
        elif euler_param == 'thtd':
            return self.model.eq_of_mot_func[7](0.,0.,0.,p,q,r,phi,tht,psi,0.,1.,1.,1.,9.81,0.,0.,0.,0.)
        
        elif euler_param == 'psid':
            return self.model.eq_of_mot_func[8](0.,0.,0.,p,q,r,phi,tht,psi,0.,1.,1.,1.,9.81,0.,0.,0.,0.)
        
        
##---------------------------------------------------------------------------##
    def __getTimeParams(self):
        self.__cIdx                        = 0
        self.__currTime                    = 0.
        self.params['t']                   = empty(0, float)
        
        
##---------------------------------------------------------------------------##
    def __keyPressEvent(self, event):
        
        if event.key() == 83:
            
            '''
            
            Key S is pressed, Start the quadcopter and rotate the propellers
            at constant RPM
            
            '''
            
            self.__executeStartManeuver()
            self.__S_cnt                   = self.__S_cnt + 1
            
        elif event.key() == 67:
            
            '''
            
            Key C is pressed, climb with a constant climb speed of 1m/s
            while maintaining level attitude
            
            '''
            
            self.__executeClimbManuever()
            self.__C_cnt                   = self.__C_cnt + 1
            self.__currAlt                 = self.x[self.__cIdx,11]
            
        elif event.key() == 68:
            
            '''
            
            Key D is pressed, descend down with a constant descend velocity 
            of -0.5m/s, while maintaining level attitude
            
            '''
            
            self.__executeDescendManuever()
            self.__D_cnt                   = self.__D_cnt + 1
            self.__currAlt                 = self.x[self.__cIdx,11]
            
        elif event.key() == 16777235:
            
            '''
            
            Key uparrow is pressed, move forward with constant pitch angle
            of 3 deg, while maintaining constant altitude
            
            '''
            
            #self.__initialize_int_err_alt()
            self.__stabilize_roll_pitch_yaw_alt_at(theta=3.0, hover_alt=self.__currAlt)
            self.__up_cnt                  = self.__up_cnt + 1
            
        elif event.key() == 16777237:
            
            '''
            
            Key downarrow is pressed, move backward with constant pitch angle
            of -3 deg, while maintaining constant altitude
            
            '''
            
            self.__stabilize_roll_pitch_yaw_alt_at(theta=-3.0, hover_alt=self.__currAlt)
            self.__down_cnt                = self.__down_cnt + 1
            
        elif event.key() == 16777234:
            
            '''
            
            Key leftarrow is pressed, move left with constant roll angle of
            3 deg, while maintaining constant altitude
            
            '''
            
            self.__stabilize_roll_pitch_yaw_alt_at(phi=3.0, hover_alt=self.__currAlt)
            self.__left_cnt                = self.__left_cnt + 1
            
        if event.key() == 16777236:
            
            '''
            
            Key rightarrow is pressed, move right with consant roll angle of
            -3 deg, while maintaining constant altitude
            
            '''
            
            self.__stabilize_roll_pitch_yaw_alt_at(phi=-3.0, hover_alt=self.__currAlt)
            self.__right_cnt               = self.__right_cnt + 1
            
        elif event.key() == 76:
            
            '''
            
            Key L is pressed, yaw left to current yaw + 3 deg, while 
            maintaining constant altitude
            
            '''
            
            self.__stabilize_roll_pitch_yaw_alt_at(psi=3.0, hover_alt=self.__currAlt)
            self.__L_cnt                   = self.__L_cnt + 1
            
        elif event.key() == 82:
            
            '''
            
            Key R is pressed, yaw right to current yaw - 3 deg, while
            maintaining constant altitude
            
            '''
            
            self.__stabilize_roll_pitch_yaw_alt_at(psi=-3.0, hover_alt=self.__currAlt)
            self.__R_cnt                   = self.__R_cnt + 1
            
        elif event.key() == 81:
            
            '''
            
            Key Q is pressed, land with small descend velocity, while 
            maintaining a level attitutde, to avoid crash
            
            '''
            self.__Q_cnt                   = self.__Q_cnt + 1
            
        elif event.key() == 72:
            
            '''
            
            Key H is pressed, hover at current altitude while maintaining 
            a level attitude
            
            '''
            
            self.__stabilize_roll_pitch_yaw_alt_at(hover_alt=self.__currAlt)
            self.__H_cnt                   = self.__H_cnt + 1
            
        elif event.key() == 69:
            sys.exit(self.app.exec_())
            
        else:
            print(event.key())
            
            
##---------------------------------------------------------------------------##
    def __executeStartManeuver(self):
        
        if self.__currTime == 0:
            if 'X0' in list(self.params.keys()):
                X0                         = self.params['X0']
            else:
                X0                         = zeros(12)
            self.x                         = append(self.x, array([X0]), axis=0)
            y_tmp, y_new                   = X0.copy(), X0.copy()
        else:
            y_tmp, y_new                   = self.x[self.__cIdx,:].copy(), self.x[self.__cIdx,:].copy()
            
        self.__flightTimeDisplay.display(self.__currTime)
        
        self.x             = append(self.x, array([y_new]), axis=0)
        Ft                 = 0.5 * (self.params['m']*self.params['g'])
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k1                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime, y_tmp, self.params, raw_loads=False)
        
        Ft                 = 0.5 * (self.params['m']*self.params['g'])
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k2                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
        
        Ft                 = 0.5 * (self.params['m']*self.params['g'])
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k3                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
        
        Ft                 = 0.5 * (self.params['m']*self.params['g'])
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k4                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
        
        y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
        
        y_tmp              = y_new.copy()
        
        self.inp           = append(self.inp, array([[self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']]]), axis=0)
        
        self.__updateAnimation()
        
        self.__currTime    = self.__currTime + self.params['dt']
        self.__cIdx        = self.__cIdx + 1
        
        print('starting')
        
        
##---------------------------------------------------------------------------##
    def __executeClimbManuever(self):
        
        if self.__currTime == 0:
            if 'X0' in list(self.params.keys()):
                X0                         = self.params['X0']
            else:
                X0                         = zeros(12)
            self.x                         = append(self.x, array([X0]), axis=0)
            y_tmp, y_new                   = X0.copy(), X0.copy()
        else:
            y_tmp, y_new                   = self.x[self.__cIdx,:].copy(), self.x[self.__cIdx,:].copy()
            
        self.__flightTimeDisplay.display(self.__currTime)
        self.__altDisplay.display(y_tmp[11])
        
        err_h              = y_tmp[11] - 0.
        derr_dt_h          = y_tmp[2] - self.__climbSpeed
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k1                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime, y_tmp, self.params, raw_loads=False)
        
        err_h              = (y_tmp + 0.5*k1)[11] - 0.
        derr_dt_h          = (y_tmp + 0.5*k1)[2] - self.__climbSpeed
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k2                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
        
        err_h              = (y_tmp + 0.5*k2)[11] - 0.
        derr_dt_h          = (y_tmp + 0.5*k2)[2] - self.__climbSpeed
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k3                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
        
        err_h              = (y_tmp + k3)[11] - 0.
        derr_dt_h          = (y_tmp + k3)[2] - self.__climbSpeed
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k4                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
        
        y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
        
        y_tmp              = y_new.copy()
        
        self.x             = append(self.x, array([y_tmp]), axis=0)
        self.inp           = append(self.inp, array([[self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']]]), axis=0)
            
        self.__updateAnimation()
        
        self.__currTime    = self.__currTime + self.params['dt']
        self.__cIdx        = self.__cIdx + 1
        
        print('u = ', y_tmp[0], ' v = ', y_tmp[1], ' w = ', y_tmp[2],
              ' p = ', y_tmp[3], ' q = ', y_tmp[4], ' r = ', y_tmp[5],
              ' phi = ', y_tmp[6], ' tht = ', y_tmp[7], ' psi = ', y_tmp[8],
              ' x = ', y_tmp[9], ' y = ', y_tmp[10], ' z = ', y_tmp[11])
        
        
##---------------------------------------------------------------------------##
    def __executeDescendManuever(self):
        
        if self.__currTime == 0:
            if 'X0' in list(self.params.keys()):
                X0                         = self.params['X0']
            else:
                X0                         = zeros(12)
            self.x                         = append(self.x, array([X0]), axis=0)
            y_tmp, y_new                   = X0.copy(), X0.copy()
        else:
            y_tmp, y_new                   = self.x[self.__cIdx,:].copy(), self.x[self.__cIdx,:].copy()
            
        self.__flightTimeDisplay.display(self.__currTime)
        self.__altDisplay.display(y_tmp[11])
        
        err_h              = y_tmp[11] - 0.
        derr_dt_h          = y_tmp[2] - self.__descendSpeed
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k1                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime, y_tmp, self.params, raw_loads=False)
        
        err_h              = (y_tmp + 0.5*k1)[11] - 0.
        derr_dt_h          = (y_tmp + 0.5*k1)[2] - self.__descendSpeed
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k2                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
        
        err_h              = (y_tmp + 0.5*k2)[11] - 0.
        derr_dt_h          = (y_tmp + 0.5*k2)[2] - self.__descendSpeed
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k3                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
        
        err_h              = (y_tmp + k3)[11] - 0.
        derr_dt_h          = (y_tmp + k3)[2] - self.__descendSpeed
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, 0.0, 0.0, self.params['kd_h'], 0.0)
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, 0., 0., 0.)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = Ft, 0., 0., 0.
        k4                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
        
        y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
        
        y_tmp              = y_new.copy()
        
        self.x             = append(self.x, array([y_tmp]), axis=0)
        self.inp           = append(self.inp, array([[self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']]]), axis=0)
            
        self.__updateAnimation()
        
        self.__currTime    = self.__currTime + self.params['dt']
        self.__cIdx        = self.__cIdx + 1
        
        print('Descending')
        
        
##---------------------------------------------------------------------------##
    def __stabilize_roll_pitch_yaw_alt_at(self, phi: float=0.0, theta: float=0., 
                                        psi:float=0., hover_alt:float=0.):
        
        if self.__currTime == 0:
            if 'X0' in list(self.params.keys()):
                X0                         = self.params['X0']
            else:
                X0                         = zeros(12)
            self.x                         = append(self.x, array([X0]), axis=0)
            y_tmp, y_new                   = X0.copy(), X0.copy()
        else:
            y_tmp, y_new                   = self.x[self.__cIdx,:].copy(), self.x[self.__cIdx,:].copy()
        
        self.__flightTimeDisplay.display(self.__currTime)
        self.__altDisplay.display(y_tmp[11])
        
        phi_rad                            = radians(phi)
        theta_rad                          = radians(theta)
        psi_rad                            = radians(psi)
        
        # if abs(y_tmp[11] - hover_alt) > self.__hover_alt_tol:
        #     warn('Altitude to hover is not close to current altitude within tolerance')
        #     return
            
        err_phi            = y_tmp[6] - phi_rad
        derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
        L                  = -pid_force1(self.__currTime, err_phi, derr_dt_phi, 0.0, self.params['kp_phi'], self.params['kd_phi'], 0.0)
        err_tht            = y_tmp[7] - theta_rad
        derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
        M                  = -pid_force1(self.__currTime, err_tht, derr_dt_tht, 0.0, self.params['kp_tht'], self.params['kd_tht'], 0.0)
        err_psi            = y_tmp[8] - psi_rad
        derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=y_tmp[3], q=y_tmp[4], r=y_tmp[5], phi=y_tmp[6], tht=y_tmp[7], psi=y_tmp[8]) - 0.
        N                  = -pid_force1(self.__currTime, err_psi, derr_dt_psi, 0.0, self.params['kp_psi'], self.params['kd_psi'], 0.0)
        err_h              = y_tmp[11] - hover_alt
        derr_dt_h          = y_tmp[2]
        self.__int_err_h   = self.__int_err_h + (self.params['dt']/self.params['tf'])*err_h
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, self.__int_err_h, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, L, M, N)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0., 0.
        k1                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime, y_tmp, self.params, raw_loads=False)
        
        err_phi            = (y_tmp + 0.5*k1)[6] - phi_rad
        derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
        L                  = -pid_force1(self.__currTime, err_phi, derr_dt_phi, 0.0, self.params['kp_phi'], self.params['kd_phi'], 0.0)
        err_tht            = (y_tmp + 0.5*k1)[7] - theta_rad
        derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
        M                  = -pid_force1(self.__currTime, err_tht, derr_dt_tht, 0.0, self.params['kp_tht'], self.params['kd_tht'], 0.0)
        err_psi            = (y_tmp + 0.5*k1)[8] - psi_rad
        derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + 0.5*k1)[3], q=(y_tmp + 0.5*k1)[4], r=(y_tmp + 0.5*k1)[5], phi=(y_tmp + 0.5*k1)[6], tht=(y_tmp + 0.5*k1)[7], psi=(y_tmp + 0.5*k1)[8]) - 0.
        N                  = -pid_force1(self.__currTime, err_psi, derr_dt_psi, 0.0, self.params['kp_psi'], self.params['kd_psi'], 0.0)
        err_h              = (y_tmp + 0.5*k1)[11] - hover_alt
        derr_dt_h          = (y_tmp + 0.5*k1)[2]
        self.__int_err_h   = self.__int_err_h + (self.params['dt']/self.params['tf'])*err_h
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, self.__int_err_h, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, L, M, N)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
        k2                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params, raw_loads=False)
        
        err_phi            = (y_tmp + 0.5*k2)[6] - phi_rad
        derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
        L                  = -pid_force1(self.__currTime, err_phi, derr_dt_phi, 0.0, self.params['kp_phi'], self.params['kd_phi'], 0.0)
        err_tht            = (y_tmp + 0.5*k2)[7] - theta_rad
        derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
        M                  = -pid_force1(self.__currTime, err_tht, derr_dt_tht, 0.0, self.params['kp_tht'], self.params['kd_tht'], 0.0)
        err_psi            = (y_tmp + 0.5*k2)[8] - psi_rad
        derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + 0.5*k2)[3], q=(y_tmp + 0.5*k2)[4], r=(y_tmp + 0.5*k2)[5], phi=(y_tmp + 0.5*k2)[6], tht=(y_tmp + 0.5*k2)[7], psi=(y_tmp + 0.5*k2)[8]) - 0.
        N                  = -pid_force1(self.__currTime, err_psi, derr_dt_psi, 0.0, self.params['kp_psi'], self.params['kd_psi'], 0.0)
        err_h              = (y_tmp + 0.5*k2)[11] - hover_alt
        derr_dt_h          = (y_tmp + 0.5*k2)[2]
        self.__int_err_h   = self.__int_err_h + (self.params['dt']/self.params['tf'])*err_h
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, self.__int_err_h, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, L, M, N)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
        k3                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params, raw_loads=False)
        
        err_phi            = (y_tmp + k3)[6] - phi_rad
        derr_dt_phi        = self.getEulerAngleRate(euler_param='phid', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
        L                  = -pid_force1(self.__currTime, err_phi, derr_dt_phi, 0.0, self.params['kp_phi'], self.params['kd_phi'], 0.0)
        err_tht            = (y_tmp + k3)[7] - theta_rad
        derr_dt_tht        = self.getEulerAngleRate(euler_param='thtd', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
        M                  = -pid_force1(self.__currTime, err_tht, derr_dt_tht, 0.0, self.params['kp_tht'], self.params['kd_tht'], 0.0)
        err_psi            = (y_tmp + k3)[8] - psi_rad
        derr_dt_psi        = self.getEulerAngleRate(euler_param='psid', p=(y_tmp + k3)[3], q=(y_tmp + k3)[4], r=(y_tmp + k3)[5], phi=(y_tmp + k3)[6], tht=(y_tmp + k3)[7], psi=(y_tmp + k3)[8]) - 0.
        N                  = -pid_force1(self.__currTime, err_psi, derr_dt_psi, 0.0, self.params['kp_psi'], self.params['kd_psi'], 0.0)
        err_h              = (y_tmp + k3)[11] - hover_alt
        derr_dt_h          = (y_tmp + k3)[2]
        self.__int_err_h   = self.__int_err_h + (self.params['dt']/self.params['tf'])*err_h
        Ft                 = (self.params['m']*self.params['g']) - pid_force1(self.__currTime, err_h, derr_dt_h, self.__int_err_h, self.params['kp_h'], self.params['kd_h'], self.params['ki_h'])
        self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4'] = self.getRPMgivenLoads(Ft, L, M, N)
        #self.params['Ft'], self.params['L'], self.params['M'], self.params['N'] = (self.params['m']*self.params['g']), L, 0, 0.
        k4                 = self.params['dt'] * self.__eval_eq_of_mot(self.__currTime + self.params['dt'], y_tmp + k3, self.params, raw_loads=False)
        
        y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
        
        y_tmp              = y_new.copy()
        
        self.x             = append(self.x, array([y_tmp]), axis=0)
        self.inp           = append(self.inp, array([[self.params['w1'], self.params['w2'], self.params['w3'], self.params['w4']]]), axis=0)
            
        self.__updateAnimation()
        
        self.__currTime    = self.__currTime + self.params['dt']
        self.__cIdx        = self.__cIdx + 1
        
        print('u = ', y_tmp[0], ' v = ', y_tmp[1], ' w = ', y_tmp[2],
              ' p = ', y_tmp[3], ' q = ', y_tmp[4], ' r = ', y_tmp[5],
              ' phi = ', y_tmp[6], ' tht = ', y_tmp[7], ' psi = ', y_tmp[8],
              ' x = ', y_tmp[9], ' y = ', y_tmp[10], ' z = ', y_tmp[11])
        
        
##---------------------------------------------------------------------------##
    def __build_shape(self, grnd_size,arm_len,arm_ang):

        ground          = BRepPrimAPI_MakeBox(gp_Pnt(-grnd_size/2., -grnd_size/2., -0.2), grnd_size, grnd_size, 0.01).Shape()
        arm1            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-self.params['dx'], -self.params['dy'], 0.0), gp_Dir(sin(arm_ang), cos(arm_ang), 0.0)), 0.01*arm_len, arm_len, 2*pi).Shape()
        arm2            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-self.params['dx'], self.params['dy'], 0.0), gp_Dir(sin(arm_ang), -cos(arm_ang), 0.0)), 0.01*arm_len, arm_len, 2*pi).Shape()
        
        mot1            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(self.params['dx'], self.params['dy'], 0.0), gp_Dir(0.0, 0.0, 1.0)), 0.05*arm_len, 0.1*arm_len, 2*pi).Shape()
        mot2            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-self.params['dx'], -self.params['dy'], 0.0), gp_Dir(0.0, 0.0, 1.0)), 0.05*arm_len, 0.1*arm_len, 2*pi).Shape()
        mot3            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(self.params['dx'], -self.params['dy'], 0.0), gp_Dir(0.0, 0.0, 1.0)), 0.05*arm_len, 0.1*arm_len, 2*pi).Shape()
        mot4            = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-self.params['dx'], self.params['dy'], 0.0), gp_Dir(0.0, 0.0, 1.0)), 0.05*arm_len, 0.1*arm_len, 2*pi).Shape()

        prop1           = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(self.params['dx']/2.0, self.params['dy'], 0.11*arm_len), gp_Dir(1.0, 0.0, 0.0)), 0.01*arm_len, self.params['dx'], 2*pi).Shape()
        prop2           = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-1.5*self.params['dx'], -self.params['dy'], 0.11*arm_len), gp_Dir(1.0, 0.0, 0.0)), 0.01*arm_len, self.params['dx'], 2*pi).Shape()
        prop3           = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(self.params['dx']/2.0, -self.params['dy'], 0.11*arm_len), gp_Dir(1.0, 0.0, 0.0)), 0.01*arm_len, self.params['dx'], 2*pi).Shape()
        prop4           = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-1.5*self.params['dx'], self.params['dy'], 0.11*arm_len), gp_Dir(1.0, 0.0, 0.0)), 0.01*arm_len, self.params['dx'], 2*pi).Shape()

        uam             = BRepAlgoAPI_Fuse(arm1, arm2).Shape()
        uam             = BRepAlgoAPI_Fuse(uam, mot1).Shape()
        uam             = BRepAlgoAPI_Fuse(uam, mot2).Shape()
        uam             = BRepAlgoAPI_Fuse(uam, mot3).Shape()
        uam             = BRepAlgoAPI_Fuse(uam, mot4).Shape()

        self.ground_shp      = self.display.DisplayColoredShape(ground, 'BLACK', update=False)[0]
        self.uam_shp         = self.display.DisplayColoredShape(uam, 'RED', update=True)[0]
        self.prop1_shp       = self.display.DisplayColoredShape(prop1, 'RED', update=True)[0]
        self.prop2_shp       = self.display.DisplayColoredShape(prop2, 'RED', update=True)[0]
        self.prop3_shp       = self.display.DisplayColoredShape(prop3, 'RED', update=True)[0]
        self.prop4_shp       = self.display.DisplayColoredShape(prop4, 'RED', update=True)[0]
        
        self.prop1_trsf_spin           = gp_Trsf()
        self.prop2_trsf_spin           = gp_Trsf()
        self.prop3_trsf_spin           = gp_Trsf()
        self.prop4_trsf_spin           = gp_Trsf()

        self.prop1_trsf_trns           = gp_Trsf()
        self.prop2_trsf_trns           = gp_Trsf()
        self.prop3_trsf_trns           = gp_Trsf()
        self.prop4_trsf_trns           = gp_Trsf()

        self.prop1_trsf_rotpsi         = gp_Trsf()
        self.prop2_trsf_rotpsi         = gp_Trsf()
        self.prop3_trsf_rotpsi         = gp_Trsf()
        self.prop4_trsf_rotpsi         = gp_Trsf()

        self.prop1_trsf_rottht         = gp_Trsf()
        self.prop2_trsf_rottht         = gp_Trsf()
        self.prop3_trsf_rottht         = gp_Trsf()
        self.prop4_trsf_rottht         = gp_Trsf()

        self.prop1_trsf_rotphi         = gp_Trsf()
        self.prop2_trsf_rotphi         = gp_Trsf()
        self.prop3_trsf_rotphi         = gp_Trsf()
        self.prop4_trsf_rotphi         = gp_Trsf()

        self.uam_trsf_trns             = gp_Trsf()
        self.uam_trsf_rotpsi           = gp_Trsf()
        self.uam_trsf_rottht           = gp_Trsf()
        self.uam_trsf_rotphi           = gp_Trsf()
        
        self.ang1                      = 0.0
        self.ang2                      = 0.0
        self.ang3                      = 0.0
        self.ang4                      = 0.0


##---------------------------------------------------------------------------##
    def __buildDRONE(self, grnd_size:float=3.0,
                      ang_fac:float=5.0,
                      fitall:bool=True,
                      w1: float=3500,
                      w2: float=3500,
                      w3: float=3500,
                      w4: float=3500):
        self.arm_len             = 2*sqrt(self.params['dx']**2 + self.params['dy']**2)
        self.arm_ang             = arctan2(self.params['dx'], self.params['dy'])
        
        self.__build_shape(grnd_size,self.arm_len,self.arm_ang)
        self.display.FitAll()
        
        
##---------------------------------------------------------------------------##
    def __updateAnimation(self):
        
        ang_fac          = 5.0
        
        xpos             = self.x[self.__cIdx,9]
        ypos             = self.x[self.__cIdx,10]
        alt              = self.x[self.__cIdx,11]
        phival           = self.x[self.__cIdx,6]
        thtval           = self.x[self.__cIdx,7]
        psival           = self.x[self.__cIdx,8]
        
        #prop_rot_axis_c  = self.model.kinematics.R_B_I_func(phival, thtval, psival) @ array([0.,0.,1.])
        #prop_rot_axis_ac = self.model.kinematics.R_B_I_func(phival, thtval, psival) @ array([0.,0.,-1.])
        
        # w1               = RPM_to_RPS(self.inp[self.__cIdx, 0])
        # w2               = RPM_to_RPS(self.inp[self.__cIdx, 1])
        # w3               = RPM_to_RPS(self.inp[self.__cIdx, 2])
        # w4               = RPM_to_RPS(self.inp[self.__cIdx, 3])
        
        w1               = RPM_to_RPS(3500)
        w2               = RPM_to_RPS(3500)
        w3               = RPM_to_RPS(3500)
        w4               = RPM_to_RPS(3500)
        
        ## propeller rotation animation
        self.ang1        = self.ang1 + w1*self.params['dt']
        self.ang2        = self.ang2 + w2*self.params['dt']
        self.ang3        = self.ang3 + w3*self.params['dt']
        self.ang4        = self.ang4 + w4*self.params['dt']
        
        self.prop1_trsf_spin.SetRotation(gp_Ax1(gp_Pnt(self.params['dx'], self.params['dy'], 0.11*self.arm_len), gp_Dir(0.,0.,1.)), self.ang1)
        self.prop1_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
        self.prop1_trsf_rotpsi.Multiply(self.prop1_trsf_spin)
        self.prop1_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
        self.prop1_trsf_rottht.Multiply(self.prop1_trsf_rotpsi)
        self.prop1_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
        self.prop1_trsf_rotphi.Multiply(self.prop1_trsf_rottht)
        self.prop1_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
        self.prop1_trsf_trns.Multiply(self.prop1_trsf_rotphi)
        self.prop1Toploc    = TopLoc_Location(self.prop1_trsf_trns)
        self.display.Context.SetLocation(self.prop1_shp, self.prop1Toploc)
        
        self.prop2_trsf_spin.SetRotation(gp_Ax1(gp_Pnt(-self.params['dx'], -self.params['dy'], 0.11*self.arm_len), gp_Dir(0.,0.,1.)), self.ang2)
        self.prop2_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
        self.prop2_trsf_rotpsi.Multiply(self.prop2_trsf_spin)
        self.prop2_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
        self.prop2_trsf_rottht.Multiply(self.prop2_trsf_rotpsi)
        self.prop2_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
        self.prop2_trsf_rotphi.Multiply(self.prop2_trsf_rottht)
        self.prop2_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
        self.prop2_trsf_trns.Multiply(self.prop2_trsf_rotphi)
        self.prop2Toploc    = TopLoc_Location(self.prop2_trsf_trns)
        self.display.Context.SetLocation(self.prop2_shp, self.prop2Toploc)
        
        self.prop3_trsf_spin.SetRotation(gp_Ax1(gp_Pnt(self.params['dx'], -self.params['dy'], 0.11*self.arm_len), gp_Dir(0.,0.,-1.)), self.ang3)
        self.prop3_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
        self.prop3_trsf_rotpsi.Multiply(self.prop3_trsf_spin)
        self.prop3_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
        self.prop3_trsf_rottht.Multiply(self.prop3_trsf_rotpsi)
        self.prop3_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
        self.prop3_trsf_rotphi.Multiply(self.prop3_trsf_rottht)
        self.prop3_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
        self.prop3_trsf_trns.Multiply(self.prop3_trsf_rotphi)
        self.prop3Toploc    = TopLoc_Location(self.prop3_trsf_trns)
        self.display.Context.SetLocation(self.prop3_shp, self.prop3Toploc)
        
        self.prop4_trsf_spin.SetRotation(gp_Ax1(gp_Pnt(-self.params['dx'], self.params['dy'], 0.11*self.arm_len), gp_Dir(0.,0.,-1.)), self.ang4)
        self.prop4_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
        self.prop4_trsf_rotpsi.Multiply(self.prop4_trsf_spin)
        self.prop4_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
        self.prop4_trsf_rottht.Multiply(self.prop4_trsf_rotpsi)
        self.prop4_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
        self.prop4_trsf_rotphi.Multiply(self.prop4_trsf_rottht)
        self.prop4_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
        self.prop4_trsf_trns.Multiply(self.prop4_trsf_rotphi)
        self.prop4Toploc    = TopLoc_Location(self.prop4_trsf_trns)
        self.display.Context.SetLocation(self.prop4_shp, self.prop4Toploc)
        
        self.uam_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
        self.uam_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), thtval*ang_fac)
        self.uam_trsf_rottht.Multiply(self.uam_trsf_rotpsi)
        self.uam_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
        self.uam_trsf_rotphi.Multiply(self.uam_trsf_rottht)
        self.uam_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
        self.uam_trsf_trns.Multiply(self.uam_trsf_rotphi)
        self.uamToploc     = TopLoc_Location(self.uam_trsf_trns)
        self.display.Context.SetLocation(self.uam_shp, self.uamToploc)
        
        self.display.Context.UpdateCurrentViewer()
        
        self.display.FitAll()
            
##---------------------------------------------------------------------------##  
    def __createHorizontalLayout(self):
        self.horizontalGroupBox   = QGroupBox("UAM_AuroraX")
        
        self.__altDisplay       = QLCDNumber(3)
        self.__altDisplay.setStyleSheet("QLCDNumber {color: rgb(0,0,0); background-color:rgb(255,255,255);}")
        self.__altDisplay.setSegmentStyle(QLCDNumber.Flat)
        self.__altDisplay.display(0)
        self.__altLabel         = QLabel('ALTITUDE')
        self.__altLabel.setAlignment(QtCore.Qt.AlignCenter)
        
        self.__rollDisplay        = QLCDNumber(3)
        self.__rollDisplay.setStyleSheet("QLCDNumber {color: rgb(0,0,0); background-color:rgb(255,255,255);}")
        self.__rollDisplay.setSegmentStyle(QLCDNumber.Flat)
        self.__rollDisplay.display(0)
        self.__rollLabel          = QLabel('Roll Angle')
        self.__rollLabel.setAlignment(QtCore.Qt.AlignCenter)
        
        self.__pitchDisplay        = QLCDNumber(3)
        self.__pitchDisplay.setStyleSheet("QLCDNumber {color: rgb(0,0,0); background-color:rgb(255,255,255);}")
        self.__pitchDisplay.setSegmentStyle(QLCDNumber.Flat)
        self.__pitchDisplay.display(0)
        self.__pitchLabel          = QLabel('Pitch Angle')
        self.__pitchLabel.setAlignment(QtCore.Qt.AlignCenter)
        
        self.__yawDisplay         = QLCDNumber(3)
        self.__yawDisplay.setStyleSheet("QLCDNumber {color: rgb(0,0,0); background-color:rgb(255,255,255);}")
        self.__yawDisplay.setSegmentStyle(QLCDNumber.Flat)
        self.__yawDisplay.display(0)
        self.__yawLabel           = QLabel('Yaw Angle')
        self.__yawLabel.setAlignment(QtCore.Qt.AlignCenter)
        
        self.__flightTimeDisplay  = QLCDNumber(3)
        self.__flightTimeDisplay.setStyleSheet("QLCDNumber {color: rgb(0,0,0); background-color:rgb(255,255,255);}")
        self.__flightTimeDisplay.setSegmentStyle(QLCDNumber.Flat)
        self.__flightTimeDisplay.display(self.__currTime)
        self.__timeLabel           = QLabel('FLIGHT TIME (SEC)')
        self.__timeLabel.setAlignment(QtCore.Qt.AlignCenter)
        
        vbox                      = QVBoxLayout()
        vbox.addWidget(self.__altLabel)
        vbox.addWidget(self.__altDisplay)
        vbox.addWidget(self.__rollLabel)
        vbox.addWidget(self.__rollDisplay)
        vbox.addWidget(self.__pitchLabel)
        vbox.addWidget(self.__pitchDisplay)
        vbox.addWidget(self.__yawLabel)
        vbox.addWidget(self.__yawDisplay)
        vbox.addWidget(self.__timeLabel)
        vbox.addWidget(self.__flightTimeDisplay)
        vbox.addStretch()
        
        layout                    = QHBoxLayout()
        
        self.canvas               = qtViewer3d(self)
        layout.addWidget(self.canvas)
        
        self.canvas.resize(800, 800)
        self.canvas.InitDriver()
        self.display              = self.canvas._display
        self.display.SetSelectionModeVertex()
        
        self.__buildDRONE()
        
        layout.addLayout(vbox)
        self.horizontalGroupBox.setLayout(layout)
        
            
##---------------------------------------------------------------------------##  
    def __startSimulator(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        windowLayout = QVBoxLayout()
        self.__createHorizontalLayout()
        
        windowLayout.addWidget(self.horizontalGroupBox)
        self.setLayout(windowLayout)
        self.show()






















































