# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 22:50:45 2015

@author: brandon
"""

<<<<<<< HEAD
from collections import OrderedDict
import itertools

import scipy.io
import os
import numpy as np

from Geometry import Geometry

PATH = os.path.join(os.path.expanduser('~'), 'Documents', 'quambo/') # Set file path
SUBS = (0, 1, 2, 3, 1, 1, 3, 1, 2, 3, 2, 3, 2, 3, 2, 3)

# The Molecule object contains a list of the geometric variants for each molecule (that is, Geometry objects).
#   a particulat theory type (LL/QU), molecule number, and set (or subset) of geometries can be specified 
class Molecule(object):
    def __init__(self, theoryType, molNum, geos = None, storeTwoElec = True):
        
        if geos is None:
            self.geos = range(1,21)
        else:
            self.geos = geos
            
        self.gSet = OrderedDict()
        self.theoryType = theoryType
        self.molNum = molNum

        self.type = SUBS[molNum -1]
        
        for i in self.geos:
            zdata = scipy.io.loadmat(PATH + 'zmats-%d-%d.mat' % (molNum, i))   
            basisSet = scipy.io.loadmat(PATH + 'basisSets-%d.mat' % (molNum)) 
            qdata = scipy.io.loadmat(PATH + '%ss-%d-%d.mat' % (theoryType, molNum, i)) 
                       
            self.gSet[i] = Geometry(zdata, basisSet, qdata, storeTwoElec = storeTwoElec)
            
            
    # getSubMatrixList creates a list of all subMatrices between two atoms in a specified matrix (oper)
    #   It is used in all getAllSubMatrices functions in the MoleculeSet class
    def getSubMatrixList(self, orbType1, orbType2, oper, atom1, atom2):
        
        subMatrixList = []
            
        # This loop loads the data, gets the shell submatrix, then
        #   adds it to a list
        for i in self.geos:
                
            list1 = self.gSet[i].onAtom(atom1, orbType1)
            list2 = self.gSet[i].onAtom(atom2, orbType2)
            subMatrixList.append(np.squeeze(self.gSet[i].getSubMatrix(oper, list1, list2)))
                       
#            print "\nsubMatrix%d:\n" % i, subMatrix
    
        return subMatrixList
        
        
    # getBondedAtoms returns a list of bonds between the specified atoms (z number) on this Molecule object
    def getBondedAtoms(self, atomZ1, atomZ2):
        bond1 = []
        bond2 = []
        theBonds = []
        zbond1 = []
        zbond2 = []
        index = []
        bonds = self.gSet[1].bond_ref[:,0]
        nums = self.gSet[1].nums[0]
        z = self.gSet[1].z[:,0]
        
        # This loop finds the positions of the z list (z numbers), nums (atom order number) and 
        #   bond.ref (what it is bonded to) lists and creates new lists of bonded atoms and thier Z numbers

        for i in range(1, len(bonds)):
            atomA = bonds[i]
            atomA_z = z[bonds[i] -1]
            atomB = nums[i]
            atomB_z = z[i]
            
            bond1.append(atomA)
            bond2.append(atomB)
            zbond1.append(atomA_z)
            zbond2.append(atomB_z)
        
        # This loop compares the list of bonded Z numbers to see if they match both atoms specified in 
        #   the initial argument
        for i in range(0, len(zbond1)):
            if (zbond1[i] == atomZ1 and zbond2[i] == atomZ2) or (zbond1[i] == atomZ2 and zbond2[i] == atomZ1):
                index.append(i)
        
        # Creates list of all matching bonds each list component is a pair of atom position numbers
        #   so for 1,2 Difluoroethane C,H (ie 6,1) bonds it would return [[1,2],[1,3],[5,6],[5,7]] 
        for i, n in enumerate(index):
            if atomZ1 == 1 and atomZ2 != 1:
                theBonds += [[bond2[n],bond1[n]]]
            else:
                theBonds += [[bond1[n],bond2[n]]]
        
        return theBonds
        
        
        
        
        
    # getNonBondedPairs returns a list of pairs between the specified atoms (z number) on this Molecule object  
    def getNonbondedPairs(self, atomZ1, atomZ2):
        index1 = []
        index2 = []
        pairs = []
        
        z = self.gSet[1].z[:,0]
        
        for i, n in enumerate(z):
            if n == atomZ1:
                index1.append(i + 1)
            elif n == atomZ2:
                index2.append(i + 1)
        

        # find pairs of the specified atoms if they are the same atoms
        if atomZ1 == atomZ2:
            pairs = list(itertools.combinations(index1, 2))

                        
        # find pairs of the specified atoms if they are different atoms
        else:
            if index1 > index2:
                temp = index1
                index1 = index2
                index2 = temp
            pairs = list(itertools.product(index1, index2))

        
        # take out any pairs that are bonds
        bonds = self.getBondedAtoms(atomZ1, atomZ2)
        end = len(pairs)
        k = 0
        
        
        for bond in bonds:
            for j in range(0, end):
                if sorted(bond) == list(sorted(pairs[j-k])):
                    del pairs[j-k]
                    end -= 1
                    k += 1

                    
        return pairs
        
    def getAtomSelves(self, atomZ):
        selves = []
        
        z = self.gSet[1].z[:,0]
        
        for i, n in enumerate(z):
            if n == atomZ:
                selves.append([i + 1,i + 1])
                
        
        return selves
        
        

#%%
#            
LL = Molecule("LL", 1)
#LL.gSet[1].printMatrix("angles", "screen")
#LL.gSet[1].printMatrix("ang_ref", "screen")
#LL.gSet[1].printMatrix("bond_ref", "screen")
#LL.gSet[1].printMatrix("bond_total", "screen")
#LL.gSet[1].printMatrix("bonds", "screen")
#LL.gSet[1].printMatrix("di_ref", "screen")
LL.gSet[1].printMatrix("KE", "file")
LL.gSet[1].printMatrix("rho", "file")
#LL.gSet[1].printMatrix("z", "screen")
LL.gSet[2].printMatrix("S", "file")
#list1 = LL.gSet[1].onAtom(1, "all")
#print list1
#list2= LL.gSet[1].onAtom(1, "all")
#print list2
#sub = LL.gSet[1].getSubMatrix("KE", list1[0,2,4], list2[0,2,4])
#print sub
#cart = LL.getCart(1)
#cart2 = []
#for i in range(0,8):
#    cart2.append(np.linalg.norm(np.array([cart[0][i], cart[1][i], cart[2][i]])))
#print cart2
    
=======
from Collection_base import Collection_base


class Molecule(Collection_base):
    '''
    Geometry objects that share a common molecular structure.
    "Common molecular structure" means the same order of chemical elements and 
    the same bondingpattern. 
    For example, Geometry.connection_order is assumed to be the same
    for all Geometry objects in the Molecule object.
    
    A Geometry is assigned a number upon addition to the Molecule object, and
    this number serves as a unique key throughout the lifetime of the
    Molecule object. 
    '''    
    
    def __init__(self, name):
        '''
        Creates a Molecule object with no geometries
        copies name into Molecule.name and initializes the internal data
        structures
        '''
        super(Molecule, self).__init__(name)
        

>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
