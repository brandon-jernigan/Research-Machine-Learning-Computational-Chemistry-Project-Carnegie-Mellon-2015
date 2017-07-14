# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 22:50:45 2015

@author: brandon
"""
#def __init__(self, zdata, basisSet, qdata):

from collections import namedtuple


from Molecule import Molecule
from Collection_base import Collection_base
#from miscFunctions import zNumToLetter, fixMolType

Atom = namedtuple('Atom',['mol','geom','atom'])
Bond = namedtuple('Bond',['mol','geom','atom1','atom2'])

class Molecule_set(Collection_base):
    '''
    Molecule objects that form a set of data for model development.
    
    A Molecule is assigned a number upon addition to the MoleculeSet object, and
    this number serves as a unique key throughout the lifetime of the
    MoleculeSet object. 
    '''    
    
    def __init__(self, name):
        '''
        Creates a MoleculeSet object with no Molecule objects
        copies name into MoleculeSet.name and initializes the internal data
        structures
        '''
        super(Molecule_set, self).__init__(name)
        
    def atoms(self, Z):
        '''
        Atoms corresponding to element Z. 
        Returns list of the namedtuple Atom
        '''
        res = []
        for imol, mol in self:
            for igeom, geom in mol:
                for iatom in geom.get_atoms(Z):
                    res.append(Atom(imol,igeom, iatom))
        return res

    def bonds(self, Z1, Z2, connection_order=1):
        '''
        Returns list of Bonds where atom1 is of type element Z1, 
        atom2 is of type Z2, and Geometry.connection_order returns the given
        connection_order. 

        Connection_order defaults to 1, meaning the atoms are bonded.

        Returns list of the namedtuple Bond
        '''
        res = []
        for imol, mol in self:
            for igeom, geom in mol:
                for bonds in geom.get_bonds(Z1,Z2,connection_order):
                    res.append(Bond(imol,igeom, bonds[0], bonds[1]))
        return res
        
    def submatrices(self, atoms_or_bonds, oper, orbs1, orbs2=None, bondAtom=None):
        '''
        Get all operator submatrices of a given list of Atom's or Bond's
        
        Required input:
        atoms_or_bonds: List of Atom or Bond namedtuples
        oper: string specifying the operator, via conventions in 
              Geometry.get_submatrix
        orbs1: orbital type associated with the Atom.atom, or Bond.atom1
               using convention of Geometry.on_atom
        orbs2: for bonds, orbital type of the second atom
        bondAtom: (0 or 1) ask for operator on first or second atom in the bond
               
        Returns:
        List of all submatrices, ordered as in atoms_or_bonds
        '''
        res = []
        if type(atoms_or_bonds[0]) is Atom:
            for atoms in atoms_or_bonds:
                geom = self[atoms.mol][atoms.geom]
                basis = geom.on_atom(atoms.atom, orbs1)
                res.append(geom.submatrix(oper,basis,basis))
        elif type(atoms_or_bonds[0]) is Bond:
            for bonds in atoms_or_bonds:
                geom = self[bonds.mol][bonds.geom]
                if bondAtom is None:
                    basis1 = geom.on_atom(bonds.atom1, orbs1)
                    basis2 = geom.on_atom(bonds.atom2, orbs2)
                    res.append(geom.submatrix(oper,basis1,basis2))
                elif bondAtom == 0:
                    basis1 = geom.on_atom(bonds.atom1, orbs1)
                    res.append(geom.submatrix(oper,basis1,basis1))
                else:
                    basis2 = geom.on_atom(bonds.atom2, orbs2)
                    res.append(geom.submatrix(oper,basis2,basis2))
        else:
            raise TypeError('Molecule_set.submatrices atoms_or_bonds type error')
        return res
