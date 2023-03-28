# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:36:26 2022

@author: I0001386

This script has the classes to develop and share RBD model of a quadrotor

RBDModelHandler1

"""


## ************************************************************************* ##
##                          Import statements                                ##
## ************************************************************************* ## 
from sympy import symbols, Matrix, sin, cos, solve, lambdify
from sympy.physics.mechanics import ReferenceFrame, dynamicsymbols, Point, inertia, dot, cross


## ************************************************************************* ##
##                    ReferenceFrame Creator class                           ##
## ************************************************************************* ##
class FrameHandler(object):
    
    '''
    
    This class creates and handles the frames required for RBD classes
    
    I: Inertial Frame
    B: Body Frame
    
    '''
    
    def __init__(self):
        self.__createFrames()
        
        
##---------------------------------------------------------------------------##
    def __createFrames(self):
        self.I                = ReferenceFrame('I')                            # inertial frame
        self.B                = ReferenceFrame('B')                            # Body frame
        
        
## ************************************************************************* ##
##                    Geometry parameters Creator class                      ##
## ************************************************************************* ##
class GeometryHandler(object):
    
    '''
    
    This class defines and handles all geometry paramaters
    
    d1x, d1y --> coordinates of thrust axis 1 in body frame
    d2x, d2y --> coordinates of thrust axis 2 in body frame
    d3x, d3y --> coordinates of thrust axis 3 in body frame
    d4x, d4y --> coordinates of thrust axis 4 in body frame
    
    '''
    
    def __init__(self, obj):
        self.__createMotorPositionParams()
        self.__createMotorPositionVectors(obj)
        
        
##---------------------------------------------------------------------------##
    def __createMotorPositionParams(self):
        self.d1x,self.d1y     = symbols('d_{1x}, d_{1y}', real=True)
        self.d2x,self.d2y     = symbols('d_{2x}, d_{2y}', real=True)
        self.d3x,self.d3y     = symbols('d_{3x}, d_{3y}', real=True)
        self.d4x,self.d4y     = symbols('d_{4x}, d_{4y}', real=True)
        
        
##---------------------------------------------------------------------------##
    def __createMotorPositionVectors(self, obj):
        self.rm1C_B           = self.d1x*obj.frames.B.x  + self.d1y*obj.frames.B.y
        self.rm2C_B           = self.d2x*obj.frames.B.x  + self.d2y*obj.frames.B.y
        self.rm3C_B           = self.d3x*obj.frames.B.x  + self.d3y*obj.frames.B.y
        self.rm4C_B           = self.d4x*obj.frames.B.x  + self.d4y*obj.frames.B.y
        
        
## ************************************************************************* ##
##                        State Vector Creator class                         ##
## ************************************************************************* ##
class StateVectorHandler(object):
    
    '''
    
    This class defines and handles all states in the statevector
    
    ub        --> x velocity of CG in body frame
    vb        --> y velocity of CG in body frame
    wb        --> z velocity of CG in body frame
    ue        --> x velocity of CG in inertial frame
    ve        --> y velocity of CG in inertial frame
    we        --> z velocity of CG in inertial frame
    p         --> roll rate in body frame
    q         --> pitch rate in body frame
    r         --> yaw rate in body frame
    phi       --> roll angle in euler zyx rotation
    theta     --> pitch angle in euler zyx rotation
    psi       --> yaw angle in euler zyx rotation
    
    '''
    
    def __init__(self):
        self.__createStateVector()
        
        
##---------------------------------------------------------------------------##
    def __createStateVector(self):
        self.phi, self.tht, self.psi          = dynamicsymbols('phi, theta, psi')
        self.ub, self.vb,self.wb              = dynamicsymbols('u_b, v_b, w_b')
        self.ue, self.ve,self.we              = dynamicsymbols('u_e, v_e, w_e')
        self.p, self.q, self.r                = dynamicsymbols('p, q, r')
        
        
## ************************************************************************* ##
##                      Points of interest Creator class                     ##
## ************************************************************************* ##
class PointsHandler(object):
    
    '''
    
    This class defines and handles all points of interest for kinematics
    kinetics
    
    CG
    
    '''
    
    def __init__(self):
        self.__createPoints()
        
        
##---------------------------------------------------------------------------##
    def __createPoints(self):
        self.CG               = Point('CG')   
        
        
## ************************************************************************* ##
##                          Kinematics Handler class                         ##
## ************************************************************************* ##
class KinematicsHandler(object):
    
    '''
    
    This class establishes all kinematics relationship between all kinematic
    parameters
    
    
    
    '''
    
    def __init__(self, obj, lin_vel_in_inertial: bool=True):
        self.__orientFrames(obj, method='EULER_ZYX')
        self.__computeLinearKinematics(obj, lin_vel_in_inertial)
        self.__computeAngularKinematics(obj)
        self.__initializeRPM()
                
        
##---------------------------------------------------------------------------##
    def __orientFrames(self, obj, method: str='EULER_ZYX'):
        if method == 'EULER_ZYX':
            obj.frames.B.orient_body_fixed(obj.frames.I, (obj.state.psi, obj.state.tht, obj.state.phi), 'ZYX')
            self.R_B_I                  = obj.frames.I.dcm(obj.frames.B)
            self.R_B_I_func             = lambdify([obj.state.phi, obj.state.tht, obj.state.psi], self.R_B_I)
        
        
##---------------------------------------------------------------------------##
    def __computeLinearKinematics(self, obj, lin_vel_in_inertial):
        if lin_vel_in_inertial:
            obj.points.CG.set_vel(obj.frames.I, obj.state.ue*obj.frames.I.x + 
                                                obj.state.ve*obj.frames.I.y +
                                                obj.state.we*obj.frames.I.z)
            self.I_vCG_I                = obj.points.CG.vel(obj.frames.I)
        else:
            obj.points.CG.set_vel(obj.frames.I, obj.state.u*obj.frames.B.x + 
                                                obj.state.v*obj.frames.B.y +
                                                obj.state.w*obj.frames.B.z)
            self.I_vCG_B                = obj.points.CG.vel(obj.frames.I)
            
        
##---------------------------------------------------------------------------##
    def __computeAngularKinematics(self, obj):
        obj.frames.B.set_ang_vel(obj.frames.I, obj.state.p*obj.frames.B.x + obj.state.q*obj.frames.B.y + obj.state.r*obj.frames.B.z)
        self.I_w_BI_B           = obj.frames.B.ang_vel_in(obj.frames.I)
        
        
##---------------------------------------------------------------------------##
    def __initializeRPM(self):
        self.w1, self.w2, self.w3, self.w4        = dynamicsymbols('omega_1, omega_2, omega_3, omega_4')        
        
        
## ************************************************************************* ##
##                        Mass parameters Creator class                      ##
## ************************************************************************* ##
class MassHandler(object):
    
    '''
    
    This class defines and handles all mass paramaters
    
    mtow,         --> takeoff mass of the quadrotor
    Ixx,          --> moment of inertia constants in body frame
    Iyy,
    Izz
    
    '''
    
    def __init__(self, obj):
        self.__createLinearInertiaParams(obj)
        self.__createRotationalInertiaParams(obj)
        
        
##---------------------------------------------------------------------------##
    def __createLinearInertiaParams(self, obj):
        self.mtow                          = symbols('m')
        
        
##---------------------------------------------------------------------------##
    def __createRotationalInertiaParams(self, obj):
        self.Ixx, self.Iyy, self.Izz       = symbols('I_{xx}, I_{yy}, I_{zz}')
        self.ICG_B                         = inertia(obj.frames.B, self.Ixx, self.Iyy, self.Izz)


## ************************************************************************* ##
##                             Loads Creator class                           ##
## ************************************************************************* ##
class LoadsHandler(object):
    
    '''
    
    This class defines various loads applied, acting on the UAM
    Fx, Fy, Fz      --> External applied loads
    L, M, N         --> External applied moments
    
    '''
    
    def __init__(self, obj, lin_vel_in_inertial: bool=True):
        self.__createAppliedForces(obj, lin_vel_in_inertial)
        self.__createAppliedMoments(obj)
        self.__createInertialForces(obj)
        self.__createInertialMoments(obj)
        self.__createGravitationalForces(obj)
        self.__createPropellerModelConstants()
        
        
##---------------------------------------------------------------------------##
    def __createAppliedForces(self, obj, lin_vel_in_inertial):
        self.Fx, self.Fy, self.Fz          = dynamicsymbols('F_x, F_y, F_z')
        self.Ft                            = dynamicsymbols('F_t')
        
        if lin_vel_in_inertial:
            self.FCG_I                     = self.Fx*obj.frames.I.x + self.Fy*obj.frames.I.y + self.Fz*obj.frames.I.z
        else:
            self.FCG_B                     = self.Fx*obj.frames.B.x + self.Fy*obj.frames.B.y + self.Fz*obj.frames.B.z
            
            
##---------------------------------------------------------------------------##
    def __createAppliedMoments(self, obj):
        self.L, self.M, self.N             = dynamicsymbols('L, M, N')
        self.MCG_B                         = self.L*obj.frames.B.x + self.M*obj.frames.B.y + self.N*obj.frames.B.z
        
        
##---------------------------------------------------------------------------##
    def __createInertialForces(self, obj):
        self.Fi                            = -obj.mass.mtow*obj.points.CG.acc(obj.frames.I)
        
        
##---------------------------------------------------------------------------##
    def __createInertialMoments(self, obj):
        self.Mi                           = -(dot(obj.frames.B.ang_acc_in(obj.frames.I), obj.mass.ICG_B) + dot(cross(obj.frames.B.ang_vel_in(obj.frames.I), obj.mass.ICG_B), obj.frames.B.ang_vel_in(obj.frames.I))) 


##---------------------------------------------------------------------------##
    def __createGravitationalForces(self, obj):
        self.g                            = symbols('g', real=True)        
        
        
##---------------------------------------------------------------------------##
    def __createPropellerModelConstants(self):
        self.km                           = symbols('k_m', real=True)
        self.bm                           = symbols('b_m', real=True)
        
        
## ************************************************************************* ##
##                        RBD Model Handler 1                                ##
## ************************************************************************* ## 
class RBDModelHandler1(object):
    
    '''
    
    This class creates sympy model for Rigid Body Dynamics of quadrotor
    The linear dynamics equations are in inertial frame
    The rotational dynamics equations are in body frame
    
    INPUTS:   
    ------
    (
        1. Sum of thrusts
        2. rolling moment
        3. pitching moment
        4. yawing moment
        
    )
        
    OUTPUTS: 
    -------
        1. Inertial x velocity in Body Frame
        2. Inertial y velocity in Body Frame
        3. Inertial z velocity in Body Frame
        4. Inertial roll rate in Body Frame
        5. Inertial picth rate in Body Frame
        6. Inertial yaw rate in Body Frame
        7. roll angle in euler ZYX convention
        8. pitch angle in euler ZYX convention
        9. yaw angle in euler ZYX convention
        
    METHOD:
    -------
        It derives the equations of motion using the Kane's dynamics
    
    '''
    
    def __init__(self):
        self.frames       = FrameHandler()
        self.geometry     = GeometryHandler(self)
        self.state        = StateVectorHandler()
        self.points       = PointsHandler()
        self.kinematics   = KinematicsHandler(self)
        self.mass         = MassHandler(self)
        self.loads        = LoadsHandler(self)
        self.__deriveKinematicEquationsOfMotion()
        self.__deriveDynamicEquationsOfMotion()
        self.__derivePositionalEquationsOfMotion()
        self.__concatEquationsOfMotion()       
        
        
##---------------------------------------------------------------------------##
    def __deriveKinematicEquationsOfMotion(self):
        mat                       = Matrix([[1, 0, -sin(self.state.tht)], 
                                            [0, cos(self.state.phi), sin(self.state.phi)*cos(self.state.tht)], 
                                            [0, -sin(self.state.phi), cos(self.state.phi)*cos(self.state.tht)]])
        mat1                      = mat.inv()
        self.kin_eq_of_mot        = mat1 * Matrix([self.state.p, self.state.q, self.state.r])
        
        
##---------------------------------------------------------------------------##
    def __deriveDynamicEquationsOfMotion(self):
        ## generalized linear velocities
        v_C_1                  = self.kinematics.I_vCG_I.diff(self.state.ue, self.frames.I, var_in_dcm=False)
        v_C_2                  = self.kinematics.I_vCG_I.diff(self.state.ve, self.frames.I, var_in_dcm=False)
        v_C_3                  = self.kinematics.I_vCG_I.diff(self.state.we, self.frames.I, var_in_dcm=False)
        v_C_4                  = self.kinematics.I_vCG_I.diff(self.state.p, self.frames.I, var_in_dcm=False)
        v_C_5                  = self.kinematics.I_vCG_I.diff(self.state.q, self.frames.I, var_in_dcm=False)
        v_C_6                  = self.kinematics.I_vCG_I.diff(self.state.r, self.frames.I, var_in_dcm=False)
        
        ## generalized angular velocites
        w_B_1                  = self.kinematics.I_w_BI_B.diff(self.state.ue, self.frames.I, var_in_dcm=False)
        w_B_2                  = self.kinematics.I_w_BI_B.diff(self.state.ve, self.frames.I, var_in_dcm=False)
        w_B_3                  = self.kinematics.I_w_BI_B.diff(self.state.we, self.frames.I, var_in_dcm=False)
        w_B_4                  = self.kinematics.I_w_BI_B.diff(self.state.p, self.frames.I, var_in_dcm=False)
        w_B_5                  = self.kinematics.I_w_BI_B.diff(self.state.q, self.frames.I, var_in_dcm=False)
        w_B_6                  = self.kinematics.I_w_BI_B.diff(self.state.r, self.frames.I, var_in_dcm=False)
        
        F1r                    = v_C_1.dot(self.loads.FCG_I) + w_B_1.dot(self.loads.MCG_B)
        F2r                    = v_C_2.dot(self.loads.FCG_I) + w_B_2.dot(self.loads.MCG_B)
        F3r                    = v_C_3.dot(self.loads.FCG_I) + w_B_3.dot(self.loads.MCG_B)
        F4r                    = v_C_4.dot(self.loads.FCG_I) + w_B_4.dot(self.loads.MCG_B)
        F5r                    = v_C_5.dot(self.loads.FCG_I) + w_B_5.dot(self.loads.MCG_B)
        F6r                    = v_C_6.dot(self.loads.FCG_I) + w_B_6.dot(self.loads.MCG_B)

        ## generalized applied forces
        Fr                     = Matrix([F1r, F2r, F3r, F4r, F5r, F6r])

        ## generalized inertia forces
        F1s                    = v_C_1.dot(self.loads.Fi) + w_B_1.dot(self.loads.Mi)
        F2s                    = v_C_2.dot(self.loads.Fi) + w_B_2.dot(self.loads.Mi)
        F3s                    = v_C_3.dot(self.loads.Fi) + w_B_3.dot(self.loads.Mi)
        F4s                    = v_C_4.dot(self.loads.Fi) + w_B_4.dot(self.loads.Mi)
        F5s                    = v_C_5.dot(self.loads.Fi) + w_B_5.dot(self.loads.Mi)
        F6s                    = v_C_6.dot(self.loads.Fi) + w_B_6.dot(self.loads.Mi)

        Frs                    = Matrix([F1s, F2s, F3s, F4s, F5s, F6s])
        
        self.dyn_eq_of_mot_raw     = Fr+Frs

        Fg_I                       = -self.mass.mtow*self.loads.g*self.frames.I.z
        
        Fg_I                       = Fg_I.to_matrix(self.frames.I)
        Ft_I                       = (self.loads.Ft*self.frames.B.z).to_matrix(self.frames.I)
        F_app                      = Fg_I + Ft_I
        self.dyn_eq_of_mot_raw     = self.dyn_eq_of_mot_raw.subs([
                                                        (self.loads.Fx, F_app[0]),
                                                        (self.loads.Fy, F_app[1]),
                                                        (self.loads.Fz, F_app[2])
        ])
        self.dyn_eq_of_mot         = solve(self.dyn_eq_of_mot_raw, (self.state.ue.diff(), 
                                                                     self.state.ve.diff(), 
                                                                     self.state.we.diff(), 
                                                                     self.state.p.diff(), 
                                                                     self.state.q.diff(), 
                                                                     self.state.r.diff()))
        

##---------------------------------------------------------------------------##
    def __derivePositionalEquationsOfMotion(self):
        self.pos_eq_of_mot         = Matrix([self.state.ue,self.state.ve,self.state.we])        
        
        
##---------------------------------------------------------------------------##
    def __concatEquationsOfMotion(self):
        self.statevector           = [
            self.state.ue,
            self.state.ve,
            self.state.we,
            self.state.p,
            self.state.q,
            self.state.r,
            self.state.phi,
            self.state.tht,
            self.state.psi
        ]
        
        self.eq_of_mot             = [
        self.dyn_eq_of_mot[self.state.ue.diff()],
        self.dyn_eq_of_mot[self.state.ve.diff()],
        self.dyn_eq_of_mot[self.state.we.diff()],
        self.dyn_eq_of_mot[self.state.p.diff()],
        self.dyn_eq_of_mot[self.state.q.diff()],
        self.dyn_eq_of_mot[self.state.r.diff()],
        self.kin_eq_of_mot[0],
        self.kin_eq_of_mot[1],
        self.kin_eq_of_mot[2],
        self.pos_eq_of_mot[0],
        self.pos_eq_of_mot[1],
        self.pos_eq_of_mot[2]
        ]
        
        self.paramslst             = [self.state.ue,
                                      self.state.ve,
                                      self.state.we,
                                      self.state.p,
                                      self.state.q,
                                      self.state.r,
                                      self.state.phi,
                                      self.state.tht,
                                      self.state.psi,
                                      self.mass.mtow,
                                      self.mass.Ixx,
                                      self.mass.Iyy,
                                      self.mass.Izz,
                                      self.loads.g,
                                      self.loads.Ft,
                                      self.loads.L,
                                      self.loads.M,
                                      self.loads.N]
            
        self.eq_of_mot_func        = []
        for eq in self.eq_of_mot:
            self.eq_of_mot_func.append(lambdify(self.paramslst, eq))
            

        print('*******************Model created successfully*****************')        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        





