'''

Python package tool chain for complete model of UAM Aurora X version 1.0

                            /\ Bx
    propeller 3             |      (Propeller 1)
                            |       / 
                        ----|----- / 
                        |   |   | / 
                        |   |   |/ 
                        |    ---|\---------> By
                        |       | \
                        ---------  \
                                    \
    propeller 2                      (Propeller 4)

'''

print('**************************UAM AuroraX*********************************')
print('''
      The Package supports creating Models for the following Keys
      -----------------------------------------------------------
      1. ReferenceFrame
      2. Geometry
      3. Kinematics
      4. Mass
      5. Loads
      6. Propeller
      7. Aerodynamics
      8. Controller
      9. RBD
      10. Model_3D
      11. Simulator
      ''')
      
print('\n')


from UAM_AuroraX_Model_v1.Utils.Model_Handler import Model
from UAM_AuroraX_Model_v1.Utils.Rigid_Body_Model import RBD_Model_Handler1
from UAM_AuroraX_Model_v1.Utils.RBD_Model_Solver import RPM_to_RPS, RPS_to_RPM, RK4_step1, RBD_Solver_Handler_No_Control_Body_Frame, RBD_Solver_Handler_PID_Control_Inertial_Frame1, RBD_Solver_Handler_No_Control_Inertial_Frame
from UAM_AuroraX_Model_v1.Utils.RBD_Model_Controller import ControlHandler
from UAM_AuroraX_Model_v1.Utils.Propeller_Thrust_estimator import PropellerThrust_BET_MT_Handler1, PropellerThrust_BET_MT_Handler2


__all__                               = (
                                         'Model',
                                         'RBD_Model_Handler1',
                                         'RPM_to_RPS',
                                         'RPS_to_RPM',
                                         'RK4_step1',
                                         'RBD_Solver_Handler_No_Control_Body_Frame',
                                         'RBD_Solver_Handler_No_Control_Inertial_Frame',
                                         'RBD_Solver_Handler_PID_Control_Inertial_Frame1',
                                         'ControlHandler',
                                         'PropellerThrust_BET_MT_Handler1',
                                         'PropellerThrust_BET_MT_Handler2'
    )







