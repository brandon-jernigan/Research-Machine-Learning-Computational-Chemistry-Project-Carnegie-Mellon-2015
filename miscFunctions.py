# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 05:32:38 2015

@author: brandon
"""
import os
from collections import OrderedDict

import numpy as np

PATH = os.path.join(os.path.expanduser('~'), 'Documents', 'quambo/') # Set file path
HARTREE_TO_KCAL = 627.509608


def convertToKcalPerMol(hartree):
    '''
    convertToKcalPerMol takes a ints, floats, tuples, lists, matrices, or dictionaries 
    of any of these and converts them from hartrees to kcal/mol
    '''
    if type(hartree) == dict or type(hartree) == OrderedDict:
        kcalMol = OrderedDict()
        for key in hartree:
            kcalMol[key] = convertToKcalPerMol(hartree[key])
            
<<<<<<< HEAD
    elif (type(hartree) == int or type(hartree) == float or type(hartree) == np.float64 or 
         (type(hartree) == np.ndarray and len(hartree.shape) == 0)):
=======
    elif type(hartree) == int or type(hartree) == float or type(hartree) == np.float64:
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        kcalMol = HARTREE_TO_KCAL * hartree
    else:
        kcalMol = [x*HARTREE_TO_KCAL for x in hartree] 
    
    return kcalMol
    
    
    

def zNumToLetter(zNum):
    '''zNumToLetter converts an atoms Z number to a Letter. 6 would become C'''
    
    values = {0: '', 1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F'}
    return values[zNum]




def fixMolType(molType, atomZ1, atomZ2):
    '''
    fixMolType takes a list of molTypes: ["C_unsub", "C_sub", "C_disub"] or 
    ["unsub", "one_sub", "one_two_disub", "one_one_disub"] and ensures their members are consistent with the
    atoms in question (a molecule with an N cannot logically be unsubstituted).
    '''
    if atomZ2 == 0 and (atomZ1 != 1 and atomZ1 != 6):
        try:
            molType.remove(0)
        except ValueError:
            pass
        return molType
        
    if (atomZ1 == 1 or atomZ1 == 6) and (atomZ2 == 1 or atomZ2 == 6):
        return molType
        

    if (atomZ1 != 6 and atomZ2 != 6) and (atomZ1 != 1 and atomZ2 !=1):
        try:
            molType.remove(1)
        except ValueError:
            pass
        
    for group2 in [7, 8, 9]:

        if (atomZ1 == group2 or atomZ2 == group2):
            try:
                molType.remove(0)
            except ValueError:
                pass
        if atomZ1 == group2 or atomZ2 == group2:
            try:
                molType.remove("C_unsub")
            except ValueError:
                pass

 
    return molType
<<<<<<< HEAD


=======
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
