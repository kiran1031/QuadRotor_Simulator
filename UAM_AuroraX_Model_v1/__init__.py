'''

Python package tool chain for complete model of UAM Aurora X version 1.0

                            /\ Bx
    propeller 4             |      (Propeller 1)
                            |       / 
                        ----|----- / 
                        |   |   | / 
                        |   |   |/ 
        By     <----------- |   |
                        |       | \
                        ---------  \
                                    \
    propeller 3                      (Propeller 2)

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
from UAM_AuroraX_Model_v1.Utils.Rigid_Body_Model import RBDModelHandler1
from UAM_AuroraX_Model_v1.Utils.RBD_Model_Solver import RBDSolverHandler1, RBDSolverHandler2
from UAM_AuroraX_Model_v1.Utils.Propeller_Thrust_estimator import PropellerThrust_BET_MT_Handler1, PropellerThrust_BET_MT_Handler2


__all__                               = (
                                         'Model',
                                         'RBDModelHandler1',
                                         'RBDSolverHandler1',
                                         'RBDSolverHandler2',
                                         'PropellerThrust_BET_MT_Handler1',
                                         'PropellerThrust_BET_MT_Handler2'
    )







