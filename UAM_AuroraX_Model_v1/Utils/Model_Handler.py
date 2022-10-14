# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 16:29:53 2022

@author: I0001386

Scripts which has handlers for creating and managing models

"""


## ************************************************************************* ##
##                          Import statements                                ##
## ************************************************************************* ## 
import sympy as sp
from sympy import symbols
from sympy.physics.mechanics import dynamicsymbols
from UAM_AuroraX_Model_v1.Utils.ReferenceFrame_Handler import ReferenceFrameModel


## ************************************************************************* ##
##                         Model Creator Class                               ##
## ************************************************************************* ## 
class Model(object):
    
    '''
    
    This class handles creation of following Models
    1. Frame
    2. Geometry
    3. Kinematics
    4. Mass
    5. Loads
    6. Propeller
    7. Aerodynamics
    8. Controller
    9. RBD
    
    '''
    
    NumModels                                     = 0
    AllModels                                     = []
    
    def __init__(self, model_tag:            str  = '',
                       model_type:           str  = '',
                       syms:                 list = [], 
                       dyn_syms:             list = [],
                       sub_models:           list = [],
                       expr:                 str  = []):
        
        '''
        
        model_tag             --> unique name for the model
        model_type            --> Type of model from the available list
        syms                  --> list of symbols to use in the model
        dyn_syms              --> list of dynamic symbols to use in the model
        sub_models            --> list of sub Models to be used in the model
        expr                  --> A sympy expression built using syms, dynsyms and submodels
        
        '''
        
        self.__list_model_types                   = [
            'FRAME',
            'GEOMETRY',
            'KINEMATICS',
            'MASS',
            'LOADS',
            'PROPELLER',
            'AERODYNAMICS',
            'CONTROLLER',
            'RBD'
            ]
        
        self.__syms                               = []
        self.__dyn_syms                           = []
        
        self.__model_tag                          = model_tag
        self.__checkModelType(model_type)
        self.__addSymsToModel(syms)
        self.__addDynSymsToModel(dyn_syms)
        self.__addSubModelsToModel(sub_models)
        self.__addExpr(expr)
        self.__createModel()
        
        Model.NumModels                           = Model.NumModels + 1
        Model.AllModels.append(self)


    ##------------------Method to check for Model type-----------------------##
    def __checkModelType(self, model_type: str=''):
        '''
        
        This Method checks if the model type entered by the user is Valid
        
        '''
        
        if model_type.upper() not in self.__list_model_types:
            raise Exception('Model Type Error --> Undefined Model type used')
            
        self.__modelType                          = model_type.upper()


    ##------------------Method to add symbols to the Model-------------------##
    def __addSymsToModel(self, syms: list=[]):
        '''
        
        This method converts all user defined dymbols into sympy symbols
        
        '''
        
        if not syms:
            return
        
        if self.__modelType                       == 'FRAME':
            self.__syms.append(syms[0])
        else:
            for sym in syms:
                if type(sym) == str:
                    sympy_sym                     = symbols(sym)
                    if sympy_sym not in self.__syms:
                        self.__syms.append(sympy_sym)
                elif type(sym) == sp.core.symbol.Symbol:
                    if sym not in self.__syms:
                        self.__syms.append(sympy_sym)


    ##------------------Method to add symbols to the Model-------------------##
    def __addDynSymsToModel(self, dyn_syms: list=[]):
        '''
        
        This method converts all user defined dymbols into sympy symbols
        
        '''
        
        if not dyn_syms:
            return
        
        for dyn_sym in dyn_syms:
            sympy_dyn_sym                         = dynamicsymbols(dyn_sym)
            if sympy_dyn_sym not in self.__dyn_syms:
                self.__dyn_syms.append(sympy_dyn_sym)


    ##------------------Method to add symbols to the Model-------------------##
    def __addSubModelsToModel(self, sub_models: list=[]):
        '''
        
        This method collects all user defined sub models
        
        '''
        
        if not sub_models:
            return


    ##---------------Method to add Expression to the Model-------------------##
    def __addExpr(self, expr: str=''):
        pass


    ##---------------------Method to create the Model------------------------##
    def __createModel(self):
        '''
        
        This Method creates the Model among the list of available models
        
        '''
        
        if self.__modelType == 'FRAME':
            self.symobj                           = ReferenceFrameModel(self.__syms[0]).createFrame()

























