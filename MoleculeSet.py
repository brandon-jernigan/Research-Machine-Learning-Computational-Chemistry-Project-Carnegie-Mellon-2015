# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 22:50:45 2015

@author: brandon
"""
#def __init__(self, zdata, basisSet, qdata):

from collections import OrderedDict
import itertools
import numpy as np

from Molecule import Molecule
from miscFunctions import zNumToLetter, fixMolType
from rotate import get_rotation_matrix


# Molecule numbers corresponding to categories, where 0 represents unsubstituted, 1 is 1-substituted
# 2 is 1,2-disubstituted, and 3 is 1,1-disubstituted

SUBS = (0, 1, 2, 3, 1, 1, 3, 1, 2, 3, 2, 3, 2, 3, 2, 3)
KEYS = {0: "unsub", 1: "one_sub", 2:"one_two_disub", 3:"one_one_disub"}



class MoleculeSet(object):
    '''
    MoleculeSet takes a theory type and a set (or subset) of molecule numbers and stores them as a set of Molecule objects
    so a MoleculeSet is a set of Molecules which is a set of Geometries
    '''    
    def __init__(self, theoryType, mols = None, molType = "all", storeTwoElec = True):      
        
        if mols is None:
            self.mols = range(1,17)
        else:
            self.mols = mols
        try:
            self.mols.remove(6)
        except ValueError:
            pass
        
        self.mSet = OrderedDict({i: Molecule(theoryType, i, storeTwoElec = storeTwoElec) for i in self.mols})
        self.theoryType = theoryType
        self.molType = molType
       

            
            

    def getAllSubMatricesBonded(self, bonds, molType, orbType, oper):
        '''
        getAllSubMatricesBonded finds all the subMatrices between across the entire MoleculeSet for the particular 
        bond, orbital type and matrix type specified.
        orbType = all 1s 2s 2p 2px 2py 2pz 2s2p s p px py pz
        oper can be: 'KE' 'H1' 'H1nuc' 'J' 'K' 'S' 'c' 'H1nuc1' 'H1nuc2' ... 
        '''        
        subMatrices = OrderedDict()
        pairType = "Bonded"
        
        if bonds == "CC":
            subMatrices = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 6, 6)
            
        elif bonds == "XC":
            subMats = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 7, 6)
            keys = subMats.keys()
            for group in keys:
                subMatrices[group] = subMats[group]
                
            subMats = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 8, 6)
            keys = subMats.keys()
            for group in keys:
                subMatrices[group] = subMats[group]
                
            subMats = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 9, 6)
            keys = subMats.keys()
            for group in keys:
                subMatrices[group] = subMats[group]
        elif bonds == "CH":
            subMatrices = self.getAllSubMatrices_bondedByC(molType, orbType, oper, 1)
            
        elif bonds == "XH":
            subMats = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 7, 1)
            keys = subMats.keys()
            for group in keys:
                subMatrices[group] = subMats[group]
                
            subMats = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 8, 1)
            keys = subMats.keys()
            for group in keys:
                subMatrices[group] = subMats[group]
        else:
            print "Error, unknown bond type"
            return
    
        return subMatrices
      
    def getAllSubMatricesNonbonded(self, pairs, molType, orbType, oper):
        '''
        getAllSubMatricesNonbonded finds all the subMatrices between across the entire MoleculeSet for the particular 
        bond, orbital type and matrix type specified.
        orbType = all 1s 2s 2p 2px 2py 2pz 2s2p s p px py pz
        oper can be: 'KE' 'H1' 'H1nuc' 'J' 'K' 'S' 'c' 'H1nuc1' 'H1nuc2' ... 
        '''        
        subMatrices = OrderedDict()
        
        if pairs == "XC_1":
            distance = 1
            atoms = list(itertools.product([7,8,9], [6]))
            for n in atoms:
                subMats = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, n[0], n[1])
                keys = subMats.keys()
                for group in keys:
                    subMatrices[group] = subMats[group]
                    
        elif pairs == "XC_2":
            distance = 2
            atoms = list(itertools.product([7,8,9], [6]))
            for n in atoms:
                subMats = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, n[0], n[1])
                keys = subMats.keys()
                for group in keys:
                    subMatrices[group] = subMats[group]  
                    
        elif pairs == "XX_1":
            distance = 1
            atoms = list(itertools.combinations_with_replacement([7, 8, 9], 2))
            for n in atoms:
                subMats = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, n[0], n[1])
                keys = subMats.keys()
                for group in keys:
                    subMatrices[group] = subMats[group]
                    
        elif pairs == "XX_2":
            distance = 2
            atoms = list(itertools.combinations_with_replacement([7, 8, 9], 2))
            for n in atoms:
                subMats = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, n[0], n[1])
                keys = subMats.keys()
                for group in keys:
                    subMatrices[group] = subMats[group]      
                    
        elif pairs == "XH_1":
            distance = 1
            atoms = list(itertools.product([7,8,9], [1]))
            for n in atoms:
                subMats = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, n[0], n[1])
                keys = subMats.keys()
                for group in keys:
                    subMatrices[group] = subMats[group]
                    
        elif pairs == "XH_2":
            distance = 2
            atoms = list(itertools.product([7,8,9], [1]))
            for n in atoms:
                subMats = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, n[0], n[1])
                keys = subMats.keys()
                for group in keys:
                    subMatrices[group] = subMats[group]
        
        elif pairs == "HH_1":
            distance = 1
            subMatrices = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, 1, 1)
        
        elif pairs == "HH_2":
            distance = 2
            subMatrices = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, 1, 1)
            
        elif pairs == "CH_1":
            distance = 1
            subMatrices = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, 6, 1)
        
        elif pairs == "CH_2":
            distance = 2
            subMatrices = self.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, orbType, oper, 6, 1)    

        else:
            print "Error, unknown bond type"
            return
    
        return subMatrices
        
    def getAllSubMatricesSelf(self, atom, molType, orbType, oper):
        '''
        getAllSubMatricesSelf finds all the subMatrices between across the entire MoleculeSet for the particular 
        atom, orbital type and matrix type specified.
        orbType = all 1s 2s 2p 2px 2py 2pz 2s2p s p px py pz
        oper can be: 'KE' 'H1' 'H1nuc' 'J' 'K' 'S' 'c' 'H1nuc1' 'H1nuc2' ... 
        '''        
        subMatrices = OrderedDict()
        pairType = "Self"
        
        if atom == "H":
            subMatrices = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 1)
        elif atom == "C":
            subMatrices = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 6)
        elif atom == "N":
            subMatrices = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 7)
        elif atom == "O":
            subMatrices = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 8)            
        elif atom == "F":
            subMatrices = self.getAllSubMatrices_ByMol(pairType, molType, orbType, oper, 9)

        return subMatrices

    def makeKeys(self, molType, atomZ1, atomZ2):
        keys = []
        molType = fixMolType(molType, atomZ1, atomZ2)
        a1 = zNumToLetter(atomZ1)
        a2 = zNumToLetter(atomZ2)
        for group in molType:
            for i, group2 in enumerate(KEYS.keys()):
                if group == group2:
                    name = KEYS[i]
                    break
            start = a1 + a2
            keys.append(start + "_" + "%s" % name)
        return keys
        
    def filterFunc(self, molType):
        '''
        filterFunc takes the current MoleculeSet and filters it according to the type of molecule specified
        molType is a list (the complete list is ["unsub", "one_sub", "one_two_disub", "one_one_disub"] another is ["all"]).
        '''  
        # create empty MoleculeSet
        molSet = MoleculeSet(self.theoryType, [])
        filtMols = []

        # create list of molecule numbers of specified type    
        for i, x in enumerate(self.mols):
            if molType == "all":
                filtMols.append(x) 
            elif SUBS[x - 1] == molType:
                filtMols.append(x)
            
        # populate the empty MoleculeSet with only the specified molecules
        for group in filtMols:
            molSet.mSet[group] = self.mSet[group]
            molSet.mols.append(group)
        
        return molSet
        
        

    def getAllSubMatrices_bondedByC(self, molType, orbType, oper, atomZ):
        '''   
        getAllSubMatrices_bondsByC gets all the bonds sorted according to what kind of carbon (how many substitutions)
        its carbon is, molType is a list (the complete list is ["C_unsub", "C_sub", "C_disub"]). atomZ is the other
        atom bonded to the carbon.
        orbType = all 1s 2s 2p 2px 2py 2pz 2s2p s p px py pz
        oper can be: 'KE' 'H1' 'H1nuc' 'J' 'K' 'H1nuc1' 'H1nuc2' ... 
        '''
        subMatrices = OrderedDict()
        
        C_disub = []
        C_sub = []
        C_unsub = []
        
        a1 = "C"
        a2 = zNumToLetter(atomZ)

        molType = fixMolType(molType, 6, atomZ)
        
        for i in self.mols:
            # molecule 6 is currently redundent, must change if we update the molecules
            
            C1 = []
            C2 = []
            
            # find C-C and C-(X or H) bonds
            CC_bonds = self.mSet[i].getBondedAtoms(6, 6)[0]            
            bonds = self.mSet[i].getBondedAtoms(6, atomZ)
            
            # Find out which atom is bonded to which carbon
            for j, bond in enumerate(bonds):
                if bond[0] == CC_bonds[0]:
                    C1 += [bond]
                elif bond[0] == CC_bonds[1]:
                    C2 += [bond]
                
                    
            # assign subMatrices based on the category of carbon it is
            if len(C1) == 1:                     
                C_disub += self.mSet[i].getSubMatrixList(orbType, orbType, oper, C1[0][0], C1[0][1])
            if len(C2) == 1:                 
                C_disub += self.mSet[i].getSubMatrixList(orbType, orbType,oper, C2[0][0], C2[0][1])
            if len(C1) == 2:
                for j in range(0, 2): 
                    C_sub += self.mSet[i].getSubMatrixList(orbType, orbType, oper, C1[j][0], C1[j][1])
            if len(C2) == 2:
                for j in range(0, 2): 
                    C_sub += self.mSet[i].getSubMatrixList(orbType, orbType, oper, C2[j][0], C2[j][1])
            if len(C1) == 3:
                for j in range(0, 3): 
                    C_unsub += self.mSet[i].getSubMatrixList(orbType, orbType, oper, C1[j][0], C1[j][1])
            if len(C2) == 3:
                for j in range(0, 3): 
                    C_unsub += self.mSet[i].getSubMatrixList(orbType, orbType, oper, C2[j][0], C2[j][1])                
        
        # create orderedDict based on the molType
        start = a1 + a2
        for group in molType:
            
            if group == "C_unsub":
                subMatrices[start + "_on_C_unsub"] = C_unsub
            elif group == "C_sub":
                subMatrices[start + "_on_C_sub"] = C_sub
            elif group == "C_disub":
                subMatrices[start + "_on_C_disub"] = C_disub
            elif group == "all":
                C_all = C_unsub + C_sub + C_disub
                subMatrices[start + "_all"] = C_all
        
        return subMatrices
        
    def getAllSubMatrices_nonbondedByAtomsBetw(self, distance, molType, orbType, oper, atomZ1, atomZ2):
        '''   
        getAllSubMatrices_bondsByC gets all the bonds sorted according to what kind of carbon (how many substitutions)
        its carbon is, molType is a list (the complete list is ["C_unsub", "C_sub", "C_disub"]). atomZ is the other
        atom bonded to the Intervcarbon.
        orbType = all 1s 2s 2p 2px 2py 2pz 2s2p s p px py pz
        oper can be: 'KE' 'H1' 'H1nuc' 'J' 'K' 'H1nuc1' 'H1nuc2' ... 
        '''

        
        a1 = zNumToLetter(atomZ1)
        a2 = zNumToLetter(atomZ2)
        subMatrixList= []
        subMatrices = OrderedDict()
        molType = fixMolType(molType, atomZ1, atomZ2)
        
        for i, group in enumerate(molType):
            subMatrixList = []
            filt = self.filterFunc(group)
#            print "\n", a1 + a2
            
            for i in filt.mols:
                # molecule 6 is currently redundent, must change if we update the molecules
                C1 = []
                C2 = []
                X1_1_C = []
                X1_2_C = []
                X2_1_C= []
                X2_2_C= []
                X1_C = []
                X2_C= []
                X1_H = []
                X2_H = []
                match = []
                XH_bonds1 = []
                XH_bonds2 = []
                # find C-C and C-(X or H) bonds
                CC_bonds = filt.mSet[i].getBondedAtoms(6, 6)[0]
                for j in range(7, 10):
                    XH_bonds1 += self.mSet[i].getBondedAtoms(j, 1)
                    XH_bonds2 += self.mSet[i].getBondedAtoms(j, 6)
                bonds1 = filt.mSet[i].getBondedAtoms(6, atomZ1)
                bonds2 = filt.mSet[i].getBondedAtoms(6, atomZ2)
                # Find out which atom is bonded to which carbon
                C2.append([CC_bonds[1], CC_bonds[0]])
                
                for j, bond in enumerate(bonds1):
                    if bond[0] == CC_bonds[0]:
                        C1 += [bond]
                    elif bond[0] == CC_bonds[1]:
                        C2 += [bond]
                        
                for j, bond in enumerate(bonds2):
                    if bond[0] == CC_bonds[0]:
                        C1 += [bond]
                    elif bond[0] == CC_bonds[1]:
                        C2 += [bond]
                        
                for bond in XH_bonds2:
                    for k in C1:
                        if bond[0] == k[0]:
                            X1_C.append(bond[1])
                    for l in C2:
                        if bond[0] == l[0]:
                            X2_C.append(bond[1])  
                            
    
                X1_C = list(set(X1_C))
                X2_C = list(set(X2_C))
                X1_C = list(set(X1_C))
                X2_C = list(set(X2_C))
    
                
                if len(X1_C) == 2:
                    X1_1_C.append(X1_C[0])
                    X1_2_C.append(X1_C[1])
                    
                if len(X2_C) == 2:
                    X2_1_C.append(X2_C[0])
                    X2_2_C.append(X2_C[1])
                    
                pairs = filt.mSet[i].getNonbondedPairs(atomZ1, atomZ2)
                
                for j, n in enumerate(C1):
                    C1[j] = n[1]
                for j, m in enumerate(C2):
                    C2[j] = m[1]
                    
                for bond in XH_bonds1:
                    for k in X1_C:
                        if bond[0] == k:
                            X1_H.append(bond[1])
                    for l in X2_C:
                        if bond[0] == l:
                            X2_H.append(bond[1])
                            
                X1_H = list(set(X1_H))
                X2_H = list(set(X2_H))
             
                            
                if distance == 1:
                    
                    nonBonded = list(itertools.combinations(C1, 2))
                    nonBonded += list(itertools.combinations(C2, 2))
                    nonBonded += list(itertools.combinations(X1_C, 2))
                    nonBonded += list(itertools.combinations(X2_C, 2))
                    nonBonded += list(itertools.product(X1_1_C, X1_2_C))
                    nonBonded += list(itertools.product(X2_1_C, X2_2_C))
                    nonBonded = list(nonBonded)
                    
                if distance == 2:    

                    C1 += X1_C                
                    C2 += X2_C
                    
                    nonBonded = list(itertools.product(C1, C2))
                    nonBonded += list(itertools.product(X1_H, C1))
                    nonBonded += list(itertools.product(X2_H, C2))
                    
                if distance == 3:
                    nonBonded = list(itertools.product(X1_C, C2))
                    nonBonded += list(itertools.product(X2_C, C1))
                    
                if distance == 4:
                    nonBonded = list(itertools.product(X1_H, X2_H))
                
                nonBonded = list(set(nonBonded))
                
                for n in pairs:
                    for m in nonBonded:
                        if sorted(n) == sorted(m):
                            match.append(n)
                match = list(set(match))
                
#                print "%d" % i, match
                for j, pair in enumerate(match):
                    if len(pair) > 0:   
                        subMatrixList += self.mSet[i].getSubMatrixList(orbType, orbType, oper, pair[0], pair[1])
                        
            for i, group2 in enumerate(KEYS.keys()):
                if group == group2:
                    name = KEYS[i]
                    break          
            start = a1 + a2
            if molType[0] == "all":
                if len(subMatrices) == 0:
                    subMatrices[start + "_%s_atom( s)_betw--%s" % (distance, molType[0])] = subMatrixList
                else:
                    subMatrices[start + "_%s_atom( s)_betw--%s" % (distance, molType[0])] += subMatrixList
            else:
                subMatrices[start + "_%s_atom(s)_betw--%s" % (distance, name)] = subMatrixList
        return subMatrices

    def getAllSubMatrices_ByMol(self, pairType, molType, orbType, oper, atomZ1, atomZ2 = 0):
        '''
        getAllSubMatrices_ByMol gets all the pairs sorted according to what kind of molecule (how many substitutions)
        its carbon is and returns the submatrices for those pairs accross the entire MolSet categorized by the 
        molType list.
        pairType is "Bonded" or "Nonbonded". 
        molType is a list (the complete list is ["unsub", "one_sub", "one_two_disub", "one_one_disub"] another is ["all"]). 
        atomZ1, atomZ2 are the z numbers for the atoms in question (C-F would be 6,9).   
        orbType = all 1s 2s 2p 2px 2py 2pz 2s2p s p px py pz
        oper can be: 'KE' 'H1' 'H1nuc' 'J' 'K' 'H1nuc1' 'H1nuc2' ... 
        '''
        subMatrices = OrderedDict()      

        a1 = zNumToLetter(atomZ1)
        a2 = zNumToLetter(atomZ2)

        molType = fixMolType(molType, atomZ1, atomZ2)
                
        # create molecule set object based on molType list using filerFunc
        for i, group in enumerate(molType):
            subMatrixList = []
            filt = self.filterFunc(group)
            
            # find pairs of bonded or unbonded atoms, based on function arguments
            for x in filt.mols:

                
                if pairType == "Bonded":
                    pairs = filt.mSet[x].getBondedAtoms(atomZ1, atomZ2)
                elif pairType == "Nonbonded":
                    pairs = filt.mSet[x].getNonbondedPairs(atomZ1, atomZ2)
                elif pairType == "Self":
                    pairs = filt.mSet[x].getAtomSelves(atomZ1)
                else:
                    print "Error, unknown pair type."
                    
                    return
                    
                # use found pairs to create a list of thier submatrices based on oper type
                for j, pair in enumerate(pairs):
                    subMatrixList += filt.mSet[x].getSubMatrixList(orbType, orbType, oper, pair[0], pair[1])

            start = a1 + a2
            # create orderedDict
            for i, group2 in enumerate(KEYS.keys()):
                if group == group2:
                    name = KEYS[i]
                    break
            if molType[0] == "all":
                if len(subMatrices) == 0:
                    subMatrices[start + "_" + "%s" % molType[0]] = subMatrixList
                else:
                    subMatrices[start + "_" + "%s" % molType[0]] += subMatrixList
            else:
                subMatrices[start + "_" + "%s" % name] = subMatrixList
            
        return subMatrices
        
    def getPairsAndSelves(self, oper, atomZ1, atomZ2, orbType1 = "all", orbType2 = "all"):
        a1 = zNumToLetter(atomZ1)
        A1 = a1 + '1'
        a2 = zNumToLetter(atomZ2)
        A2 = a2 + '2'
        bond = a1 + a2
        
        subMatrices = OrderedDict()
        one_selves = []
        two_selves = []
        bond_submatrices = []
        one_submatrices = []
        two_submatrices = []
        for i, n in enumerate(self.mSet):
            
            pairs = self.mSet[n].getBondedAtoms(atomZ1, atomZ2)
            one_selves = []
            two_selves = []
            for j, m in enumerate(pairs):
                
                one_selves.append(m[0])
        
            
            for j, m in enumerate(pairs):
                
                two_selves.append(m[1])  
        
            
            for j, m in enumerate(pairs):
                
                bond_submatrices += self.mSet[n].getSubMatrixList(orbType1, orbType2, oper, m[0], m[1])
                one_submatrices += self.mSet[n].getSubMatrixList(orbType1, orbType2, oper, one_selves[j], one_selves[j])
                two_submatrices += self.mSet[n].getSubMatrixList(orbType1, orbType2, oper, two_selves[j], two_selves[j])
        
        subMatrices[bond] = bond_submatrices
        subMatrices[A1] = one_submatrices
        subMatrices[A2] = two_submatrices
        
        return subMatrices
        
        
        
        
#    def makeVectors(self, atomZ1, atomZ2):
#        if ((atomZ1 == 6 or atomZ1 == 7 or atomZ1 == 8 or atomZ1 == 9) and 
#            (atomZ1 == 6 or atomZ1 == 7 or atomZ1 == 8 or atomZ1 == 9)):
#                
#            subMatrices = self.getPairsAndSelves("S", atomZ1, atomZ2, "2p", "2s")
#        else:
#            subMatrices = self.getPairsAndSelves("S", atomZ1, atomZ2, "2p", "1s")
#            pairs = self.mSet[n].getBondedAtoms(atomZ1, atomZ2)
#        vectorList = OrderedDict()
#        for group in subMatrices.keys():
#            vectorList[group] =  OrderedDict([(1, []), (2, [])])
#            for mat in subMatrices[group][0]:
#                mat = np.squeeze(mat)
#                vectorList[group][1].append(mat)
##                vectorList[group][2].append(-1 * mat)
#        
#                    
#        return vectorList


    def makeVectors(self, atomZ1, atomZ2):

        vector_list = []
        for i in self.mols:
            mol_set = self.mSet[i]
            pairs = mol_set.getBondedAtoms(atomZ1, atomZ2)
            for j in self.mSet[i].geos:
                geo_set = mol_set.gSet[j]
                cart_list = geo_set.cartesian
                for n, pair in enumerate(pairs):
                    vector_list.append(np.array([cart_list[0][pair[1] - 1] - cart_list[0][pair[0] - 1],
                                         cart_list[1][pair[1] - 1] - cart_list[1][pair[0] - 1],
                                         cart_list[2][pair[1] - 1] - cart_list[2][pair[0] - 1]]))
                

                    
        return vector_list
        
        
    def make_rotated_matrix(self, pairType, molType, orbType, oper, atomZ1, atomZ2):
        
        a1 = zNumToLetter(atomZ1)
        A1 = a1 + '1'
        a2 = zNumToLetter(atomZ2)
        A2 = a2 + '2'
        pair = a1 + a2
        R1 = OrderedDict()
        R2 = OrderedDict()
        mat_dict = OrderedDict()
        
    #    np.set_printoptions(precision = 6, suppress = True)
        mol_set = MoleculeSet(self.theoryType, storeTwoElec = False)
        
        vectors = mol_set.makeVectors(atomZ1, atomZ2)

        
        rhos = mol_set.getPairsAndSelves("rho", atomZ1, atomZ2, "2p", "2p")

        
        R1, R2 = [], []
        
        for vec, rho1, rho2 in zip(vectors, rhos[A1], rhos[A2]):
    #        print vec, rho
            if atomZ1 == 1:
                R1.append(1)
            else:
                R1.append(get_rotation_matrix(vec, rho1))
            
            if atomZ2 == 1:
                R2.append(1)
            else:
                R2.append(get_rotation_matrix(-1 * vec, rho2)) # double check
    
        old_operators_2p = mol_set.getPairsAndSelves(oper, atomZ1, atomZ2, "2p", "2p")

        
        new_operators = mol_set.getPairsAndSelves(oper, atomZ1, atomZ2, "all")

        
        for group in [pair, A1, A2]:
            mat_list = []

            i = 0
            for ops, mat in zip(old_operators_2p[group], new_operators[group]):
                if type(R2[i]) != int:
                    r2 = R2[i].T
                else:
                    r2 = R2[i]
                    
                if mat.size == 1:
                    mat = ops
                    
                elif mat.size == 5:
                    ops = np.asmatrix(ops)
                    if atomZ2 == 1:
                        ops = ops.T

                    mat[[2,3,4]] =  R1[i] * ops * r2
                    
                elif mat.size == 25:            
                    ops = np.asmatrix(ops)
                    
                    mat[2:,2:] =  R1[i] * ops * r2
                    
                mat_list.append(mat)
                i += 1
            mat_dict[group] = mat_list
                
        return mat_dict
        
#%%
#LL = MoleculeSet("LL", [1])
#g = LL.mSet[2].gSet[1]
#g.printMatrix("rho", "file")
#g.printMatrix("KE", "file")
#g.printMatrix("H1nuc", "file")
#g.printMatrix("J", "file")
#g.printMatrix("K", "file")
#for i in range(1,17):
#    print i
#    print LL.mSet[i].gSet[1].nums
#    a = LL.mSet[i].gSet[1].z
#    print a.reshape(1,a.size)
#    print "\n\n"
#molType = ["all"]
#molType = [0, 1, 2, 3]
#atomZ1 = 6
#oper = "KE"
#orbType = "all"
#atomZ2 = 7
#
#distance = 3
#print LL.getAllSubMatrices_ByMol("Self", molType, orbType, oper, atomZ1)
#mol1  = LL.getAllSubMatrices_nonbondedByAtomsBetw(distance, molType, "2s2p", "KE", atomZ1, atomZ2)

#print mol1
#print len(mol1)
#vectors =  LL.makeVectors(6,6)
#print LL.mSet[1].gSet[1].cartesian
#print vectors
#print len(vectors)
#pairType = "Bonded"
#molType = ["all"]
#orbType = "all"
#oper = "KE"
#atomZ1 = 1
#atomZ2 = 6
#
#
#print LL.make_rotated_matrix(pairType, molType, orbType, oper, atomZ1, atomZ2)
