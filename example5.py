# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 00:26:21 2022

@author: I0001386
"""


import sys, os
sys.path.append(r'C:\Users\i0001386\Desktop\UAM_AURORA')

from PyQt5.QtWidgets import QApplication
from numpy import array, radians
from matplotlib.pyplot import figure, gca, show
from matplotlib.ticker import StrMethodFormatter
import logging
import UAM_AuroraX_Model_v1 as uamxv1

logging.getLogger('matplotlib.font_manager').disabled = True

rbd_model                     = uamxv1.RBDModelHandler1()


params                        = {
    
    ## mass parameters
    'm'                       : 0.468,             # takeoff mass in kg
    'Ixx'                     : 4.856e-3,          # Moment of inertia kgm^2
    'Iyy'                     : 4.856e-3,          # Moment of inertia
    'Izz'                     : 8.801e-3,          # Moment of inertia (Assume nearly flat object, z=0)
    
    'g'                       : 9.81,              # acceleration due to gravity
    
    ## geometry parameters
    'd1x'                     : 0.225,             # x offset of thrust axis
    'd1y'                     : 0.0000,            # y offset of thrust axis
    'd2x'                     : 0.0000,            # x offset of thrust axis
    'd2y'                     : -0.225,            # y offset of thrust axis
    'd3x'                     : -0.225,            # x offset of thrust axis
    'd3y'                     : 0.0000,            # y offset of thrust axis
    'd4x'                     : 0.000 ,            # x offset of thrust axis
    'd4y'                     : 0.225,             # y offset of thrust axis
    
    'dx'                      : 0.225,
    'dy'                      : 0.225,
    
    ## loads parameters
    'km'                      : 2.98e-6,           # thrust coeff of motor propeller
    'bm'                      : 0.114e-6,          # yawing moment coeff of motor propeller
    
    ## motor rpms
    'w1'                      : 6000,
    'w2'                      : 6000 + 10,
    'w3'                      : 6000,
    'w4'                      : 6000 - 10,
    
    ## simulation parameters
    'dt'                      : 0.02,              # Sampling time (sec)
    'tf'                      : 10.0,              # Length of time to run simulation (sec),
    'X0'                      : array([0,            # u0
                                       0,            # v0
                                       0.,           # w0
                                       0,            # p0
                                       0,            # q0
                                       0,            # r0
                                       radians(0.),  # phi0
                                       radians(0.),  # tht0
                                       radians(0.),  # psi0
                                       0.,            # x0
                                       0.,            # y0
                                       0.]),          # z0 
    
    ## control parameters
    'kp_h'                    : 20,                 # proportional constant for altitude control
    'kd_h'                    : 5,                  # derivative constant for altitude control   
    'ki_h'                    : 50,

    'kp_tht'                  : 1,
    'kd_tht'                  : 0.12,
    'ki_tht'                  : 2.,
    
    'kp_phi'                  : 3,
    'kd_phi'                  : 0.12,
    'ki_phi'                  : 0.,
    
    'kp_psi'                  : 1,
    'kd_psi'                  : 0.12,
    'ki_psi'                  : 0.,
    
    'kp_x'                    : 0.1,
    'kd_x'                    : 0.2,
    'ki_x'                    : 0.,
    
    'kp_y'                    : -0.25,
    'kd_y'                    : -0.2,
    'ki_y'                    : 0.,
    
    'kp_p'                    : 0.25,
    'kd_p'                    : 0.001,
    'ki_p'                    : 0.3,
    
    'kp_q'                    : 0.25,
    'kd_q'                    : 0.001,
    'ki_q'                    : 0.3,
    
    'kp_r'                    : 0.25,
    'kd_r'                    : 0.001,
    'ki_r'                    : 0.3
    
}

run_option = 1


if run_option == 1:
    rbd_solver                   = uamxv1.RBDSolverHandler1(rbd_model, params=params)
    # #rbd_solver.solve()
    # rbd_solver.hover_at(hover_alt=5.0, ti=0.0, tf=10.0)
    # rbd_solver.climb_and_hover_at(hover_alt=5, climb_speed=1.)
    # rbd_solver.stabilize_pitch_at(theta=3.)
    # #rbd_solver.stabilize_roll_at(phi=2.)
    # #rbd_solver.stabilize_yaw_at(psi=0.)
    #rbd_solver.stabilize_roll_pitch_yaw_alt_at(phi=0., theta=-10, psi=0., hover_alt=1.)
    # rbd_solver.stabilize_pqr_alt_at(p_des=0.0174533, q_des=0.0174533, r_des=0, zd=1.)
    rbd_solver.move_to_XYZ(xd=1., yd=2., zd=0.1, psi_des=30., conv_tol=0.5)
    
    
    # print('steady state error in Altitude = ', ((5.0 - rbd_solver.x[-1,11])/5.0)*100)
    
    # vars_to_plot = ['phi', 'tht', 'psi', 'xe', 'ye', 'ze', 'u', 'v', 'w', 'p', 'q', 'r']
    # #vars_to_plot = ['psi', 'r']
    # for var in vars_to_plot:
    #     figure(rbd_solver.plotNo)
    #     rbd_solver.plotter_with_time(yvar=var)
    #     # if var == 'ze':
    #     #     gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.6f}'))
        
    # show()
    
    rbd_solver.animate(grnd_size = 5, ang_fac=1.0, fitall=True)

elif run_option == 2:
    if __name__ == '__main__':
        app = QApplication(sys.argv)
        ex = uamxv1.RBDSolverHandler2(rbd_model, app, params=params)
        if os.getenv('APPVEYOR') is None:
            sys.exit(app.exec_())


# python C:\Users\kiran\Desktop\UAM_AuroraX\example5.py














