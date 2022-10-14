# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:47:17 2022

@author: i0001386

This script solves or integrates the equations of motion 

"""


## ************************************************************************* ##
##                          Import statements                                ##
## ************************************************************************* ## 
from numpy import arange, zeros, size, pi, sqrt, arctan2, cos, sin, array
from matplotlib.pyplot import figure, plot, xlabel, ylabel, grid, show
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
    return kp * err + kd * derr_dt


## ************************************************************************* ##
##           RBD Solver Class (Linear equations in body frame)               ##
## ************************************************************************* ##
class RBD_Solver_Handler_No_Control_Body_Frame(object):
    
    '''
    
    This class integrates equations of motion without any controller
    or uses open source controller
    
    The linear equations of motion are in the body frame
    The angular equations of motion are also in the body frame
    Positional equations transforms position vector from body frame to 
    inertial frame
    
    '''

    def __init__(self, rbd_model,
                       params: dict={}):
        self.model                         = rbd_model
        self.params                        = params
        self.__getTimeParams()
        self.n_states                      = 12                # Number of states
        self.n_inputs                      = 4                 # Number of inputs
        self.x                             = zeros((size(self.params['t']), self.n_states))      # time history of state vectors
        self.inp                           = zeros((size(self.params['t']), self.n_inputs))      # time history of input vectors
        
        
##---------------------------------------------------------------------------##
    def __getTimeParams(self):
        self.params['t']                   = arange(0, self.params['tf'], self.params['dt'])
        
        
##---------------------------------------------------------------------------##
    def __eval_eq_of_mot(self, t, y, params):
        ## expand the params
        mv                                 = params['m']
        Ixxv                               = params['Ixx']
        Iyyv                               = params['Iyy']
        Izzv                               = params['Izz']
        gv                                 = params['g']
        Ftv                                = params['Ft/mg'] * params['m'] * params['g']
        Lv                                 = params['L']
        Mv                                 = params['M']
        Nv                                 = params['N']
        
        uv                                 = y[0]
        vv                                 = y[1]
        wv                                 = y[2]
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
            output[i]                      = self.model.eq_of_mot_func[i](uv,vv,wv,pv,qv,rv,phiv,thtv,psiv,mv,Ixxv,Iyyv,Izzv,gv,Ftv,Lv,Mv,Nv)
        
        return output
        
    
##---------------------------------------------------------------------------##
    def solve(self):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        self.x[0,:]                        = X0
        for t_idx in range(len(self.params['t'])-1):
            y_tmp                          = RK4_step1(self.__eval_eq_of_mot, self.x[t_idx,:], self.params['t'][t_idx], self.params['dt'], self.params)
            self.x[t_idx+1, :]             = y_tmp
        

## ************************************************************************* ##
##           RBD Solver Class (Linear equations in inertial frame)           ##
## ************************************************************************* ##
class RBD_Solver_Handler_No_Control_Inertial_Frame(object):
    
    '''
    
    This class integrates equations of motion without any controller
    
    The linear equations of motion are in the inertial frame
    The angular equations of motion are also in the body frame
    
    '''

    def __init__(self, rbd_model,
                       params: dict={}):
        self.model                         = rbd_model
        self.params                        = params
        self.__getTimeParams()
        self.n_states                      = 12                # Number of states
        self.n_inputs                      = 4                 # Number of inputs
        self.x                             = zeros((size(self.params['t']), self.n_states))      # time history of state vectors
        self.inp                           = zeros((size(self.params['t']), self.n_inputs))      # time history of input vectors
        
        
##---------------------------------------------------------------------------##
    def __getTimeParams(self):
        self.params['t']                   = arange(0, self.params['tf'], self.params['dt'])
        
        
##---------------------------------------------------------------------------##
    def __eval_eq_of_mot(self, t, y, params):
        ## expand the params
        mv                                 = params['m']
        Ixxv                               = params['Ixx']
        Iyyv                               = params['Iyy']
        Izzv                               = params['Izz']
        gv                                 = params['g']
        Ftv                                = params['Ft/mg'] * params['m'] * params['g']
        Lv                                 = params['L']
        Mv                                 = params['M']
        Nv                                 = params['N']
        
        uve                                = y[0]
        vve                                = y[1]
        wve                                = y[2]
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
            output[i]                      = self.model.eq_of_mot_func[i](uve,vve,wve,pv,qv,rv,phiv,thtv,psiv,mv,Ixxv,Iyyv,Izzv,gv,Ftv,Lv,Mv,Nv)
        
        return output
        
    
##---------------------------------------------------------------------------##
    def solve(self):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        self.x[0,:]                        = X0
        for t_idx in range(len(self.params['t'])-1):
            
            y_tmp                          = RK4_step1(self.__eval_eq_of_mot, self.x[t_idx,:], self.params['t'][t_idx], self.params['dt'], self.params)
            self.x[t_idx+1, :]             = y_tmp    

        print('*******************Solver ran successfully********************')


##---------------------------------------------------------------------------##
    def plotter_with_time(self, yvar:str='h'):
        
        if yvar == 'h':
            figure(1)
            plot(self.params['t'], self.x[:, 11])
            grid()
            xlabel('time(seconds)')
            ylabel('Altitude(meters)')
            show()
            
        if yvar == 'x':
            figure(2)
            plot(self.params['t'], self.x[:, 9])
            grid()
            xlabel('time(seconds)')
            ylabel('x translation(meters)')
            show()
            
        if yvar == 'y':
            figure(3)
            plot(self.params['t'], self.x[:, 10])
            grid()
            xlabel('time(seconds)')
            ylabel('y translation(meters)')
            show()
            
            
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
            prop1_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
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
            prop2_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
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
            prop3_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
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
            prop4_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
            prop4_trsf_rottht.Multiply(prop4_trsf_rotpsi)
            prop4_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
            prop4_trsf_rotphi.Multiply(prop4_trsf_rottht)
            prop4_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
            prop4_trsf_trns.Multiply(prop4_trsf_rotphi)
            prop4Toploc    = TopLoc_Location(prop4_trsf_trns)
            display.Context.SetLocation(prop4_shp, prop4Toploc)
            
            uam_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
            uam_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
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
##       RBD Solver with Controller (Linear equations in inertial frame)     ##
## ************************************************************************* ##
class RBD_Solver_Handler_PID_Control_Inertial_Frame1(object):
    
    '''
    
    This class integrates equations of motion with PID control
    for manuevers specified.
    
    The linear equations of motion are in the inertial frame
    The angular equations of motion are also in the body frame
    
    '''

    def __init__(self, rbd_model,
                       rbd_control: list=[],
                       params: dict={}):
        self.model                         = rbd_model
        self.control                       = rbd_control
        self.params                        = params
        self.__initializePlottingDict()
        self.__getTimeParams()
        self.n_states                      = 12                # Number of states
        self.n_inputs                      = 4                 # Number of inputs
        self.x                             = zeros((size(self.params['t']), self.n_states))      # time history of state vectors
        self.inp                           = zeros((size(self.params['t']), self.n_inputs))      # time history of input vectors
        self.plotNo                        = 1


##---------------------------------------------------------------------------##
    def __getTimeParams(self):
        self.params['t']                   = arange(0, self.params['tf'], self.params['dt'])
        
        
##---------------------------------------------------------------------------##
    def __eval_eq_of_mot(self, t, y, params):
        ## expand the params
        mv                                 = params['m']
        Ixxv                               = params['Ixx']
        Iyyv                               = params['Iyy']
        Izzv                               = params['Izz']
        gv                                 = params['g']
        Ftv                                = params['Ft/mg'] * params['m'] * params['g']
        Lv                                 = params['L']
        Mv                                 = params['M']
        Nv                                 = params['N']
        
        uv                                 = y[0]
        vv                                 = y[1]
        wv                                 = y[2]
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
            output[i]                      = self.model.eq_of_mot_func[i](uv,vv,wv,pv,qv,rv,phiv,thtv,psiv,mv,Ixxv,Iyyv,Izzv,gv,Ftv,Lv,Mv,Nv)
        
        return output
        
    
##---------------------------------------------------------------------------##
    def solve(self):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        self.x[0,:]                        = X0
        
        for t_idx in range(len(self.params['t'])-1):
            
            #self.params,control_identified = self.control.update_control(self.params['t'][t_idx], self.x[t_idx,:], self.params)
            
            y_tmp                          = RK4_step1(self.__eval_eq_of_mot, self.x[t_idx,:], self.params['t'][t_idx], self.params['dt'], self.params)
            
            ## restricting altitude to physical solutions
            if y_tmp[11] < 0:
                break
            
            self.x[t_idx+1, :]             = y_tmp
        
        print('*******************Solver ran successfully********************')
        
        
##---------------------------------------------------------------------------##
    def altitude_control_test1(self, hover_alt=0.0):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        self.x[0,:]                        = X0
        
        prev_err_h                         = 0.0
        y_tmp, y_new                       = X0.copy(), X0.copy()
        
        for t_idx in range(len(self.params['t'])):
            
            t0                 = self.params['t'][t_idx]
            self.x[t_idx, :]   = y_new
            err_h              = y_tmp[11] - hover_alt
            derr_dt_h          = (err_h - prev_err_h)/self.params['dt']
            self.params['Ft/mg']  = 1.0 - (pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0))/(self.params['m']*self.params['g'])
            k1                 = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params)
            
            err_h              = (y_tmp + 0.5*k1)[11] - hover_alt
            derr_dt_h          = (err_h - prev_err_h)/self.params['dt']
            self.params['Ft/mg']  = 1.0 - (pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0))/(self.params['m']*self.params['g'])
            k2                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params)
            
            err_h              = (y_tmp + 0.5*k2)[11] - hover_alt
            derr_dt_h          = (err_h - prev_err_h)/self.params['dt']
            self.params['Ft/mg']  = 1.0 - (pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0))/(self.params['m']*self.params['g'])
            k3                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params)
            
            err_h              = (y_tmp + k3)[11] - hover_alt
            derr_dt_h          = (err_h - prev_err_h)/self.params['dt']
            self.params['Ft/mg']  = 1.0 - (pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0))/(self.params['m']*self.params['g'])
            k4                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + self.params['dt'], y_tmp + k3, self.params)
            
            y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp              = y_new.copy()
            
            prev_err_h         = err_h
            
            
##---------------------------------------------------------------------------##
    def altitude_control_test2(self, hover_alt=0.0):
        if 'X0' in list(self.params.keys()):
            X0                             = self.params['X0']
        else:
            X0                             = zeros(12)
        self.x[0,:]                        = X0
        
        y_tmp, y_new                       = X0.copy(), X0.copy()
        
        for t_idx in range(len(self.params['t'])):
            
            t0                 = self.params['t'][t_idx]
            self.x[t_idx, :]   = y_new
            err_h              = y_tmp[11] - hover_alt
            derr_dt_h          = y_tmp[2]
            self.params['Ft/mg']  = 1.0 - (pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0))/(self.params['m']*self.params['g'])
            k1                 = self.params['dt'] * self.__eval_eq_of_mot(t0, y_tmp, self.params)
            
            err_h              = (y_tmp + 0.5*k1)[11] - hover_alt
            derr_dt_h          = (y_tmp + 0.5*k1)[2]
            self.params['Ft/mg']  = 1.0 - (pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0))/(self.params['m']*self.params['g'])
            k2                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k1, self.params)
            
            err_h              = (y_tmp + 0.5*k2)[11] - hover_alt
            derr_dt_h          = (y_tmp + 0.5*k2)[2]
            self.params['Ft/mg']  = 1.0 - (pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0))/(self.params['m']*self.params['g'])
            k3                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + 0.5*self.params['dt'], y_tmp + 0.5*k2, self.params)
            
            err_h              = (y_tmp + k3)[11] - hover_alt
            derr_dt_h          = (y_tmp + k3)[2]
            self.params['Ft/mg']  = 1.0 - (pid_force1(t0, err_h, derr_dt_h, 0.0, self.params['kp_h'], self.params['kd_h'], 0.0))/(self.params['m']*self.params['g'])
            k4                 = self.params['dt'] * self.__eval_eq_of_mot(t0 + self.params['dt'], y_tmp + k3, self.params)
            
            y_new              = y_tmp + (1./6.)* (k1 + 2*k2 + 2*k3 + k4)
            
            y_tmp              = y_new.copy()
                
            
##---------------------------------------------------------------------------##
    def __initializePlottingDict(self):
        self.__plotDict        = {
            'u'     :      (0,  'x velocity(m/s)'),
            'v'     :      (1,  'y velocity(m/s)'),
            'w'     :      (2,  'climb velocity(m/s)'),
            'p'     :      (3,  'roll rate(rad/s)'),
            'q'     :      (4,  'pitch rate(rad/s)'),
            'r'     :      (5,  'yaw rate(rad/s)'),
            'phi'   :      (6,  'roll Angle(deg)'),
            'tht'   :      (7,  'pitch Angle(deg)'),
            'psi'   :      (8,  'yaw Angle(deg)'),
            'x'     :      (9,  'x (m)'),
            'y'     :      (10, 'y(m)'),
            'h'     :      (11, 'Altitude(m)')
            }
            

##---------------------------------------------------------------------------##
    def plotter_with_time(self, yvar:str='h'):
        
        self.plotNo = self.plotNo + 1
        plot(self.params['t'], self.x[:, self.__plotDict[yvar][0]])
        grid()
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
            prop1_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
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
            prop2_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
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
            prop3_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
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
            prop4_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
            prop4_trsf_rottht.Multiply(prop4_trsf_rotpsi)
            prop4_trsf_rotphi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.)), phival*ang_fac)
            prop4_trsf_rotphi.Multiply(prop4_trsf_rottht)
            prop4_trsf_trns.SetTranslation(gp_Vec(xpos, ypos, alt))
            prop4_trsf_trns.Multiply(prop4_trsf_rotphi)
            prop4Toploc    = TopLoc_Location(prop4_trsf_trns)
            display.Context.SetLocation(prop4_shp, prop4Toploc)
            
            uam_trsf_rotpsi.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), psival*ang_fac)
            uam_trsf_rottht.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.)), -thtval*ang_fac)
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








