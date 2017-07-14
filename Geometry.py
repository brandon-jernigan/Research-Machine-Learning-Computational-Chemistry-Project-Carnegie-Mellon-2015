# -*- coding: utf-8 -*-

import numpy as np

<<<<<<< HEAD
import generateFiles as gen
#path = os.path.join(os.path.expanduser('~'), 'Documents', 'quambo/') # Set file path

np.set_printoptions(threshold='nan')  # Will show entire matrix in interpreter


class Geometry(object):
    """Molecule geometry and data"""

    def __init__(self, zdata, basisSet, qdata, storeTwoElec = True ):
        """Make from matlab files specifying basis set, zmatrix, and quantum data."""
        
        self.angles                = zdata['angles']
        self.ang_ref               = zdata['ang_ref']
        self.bond_ref              = zdata['bond_ref']
        self.bond_total            = zdata['bond_total']
        self.bonds                 = zdata['bonds']
        self.di_ref                = zdata['di_ref']
        self.dihedrals             = zdata['dihedrals']
        self.nums                  = zdata['nums']
        self.z                     = zdata['z']

        self.FuncToAngular          = basisSet['FuncToAngular']
        self.FuncToCenter           = basisSet['FuncToCenter']
        self.FuncToShell            = basisSet['FuncToShell']
        self.IsSpherical            = basisSet['IsSpherical']
        self.Name                   = basisSet['Name']
        self.NumFunctions           = basisSet['NumFunctions']
        self.NumShells              = basisSet['NumShells']
        self.PrimCoeffUnnorm        = basisSet['PrimCoeffUnnorm']
        self.PrimExp                = basisSet['PrimExp']
        self.ShellNumFunctions      = basisSet['ShellNumFunctions']
        self.ShellNumPrimitives     = basisSet['ShellNumPrimitives']
        self.ShellToCenter          = basisSet['ShellToCenter']
        self.ShellTypes             = basisSet['ShellTypes']
        
        self.corePotentialMats     = qdata['corePotentialMats']
        self.coulombMat            = qdata['coulombMat']
        self.densvec               = qdata['densvec']
        self.rho                   = 2 * self.densvec.reshape(self.NumFunctions, self.NumFunctions)
        self.exchangeMat           = qdata['exchangeMat']
        self.hfEnergy              = qdata['hfEnergy']
        self.kineticMat            = qdata['kineticMat']
        self.nucRepEnergy          = qdata['nucRepEnergy']
        self.numElectrons          = qdata['numElectrons']
        self.orbital               = qdata['orbital']
        self.orbitalEnergies       = qdata['orbitalEnergies']
        self.overlapMat            = qdata['overlapMat']
        if (storeTwoElec):
            self.twoElecIntegrals      = qdata['twoElecIntegrals']
        else:
            self.twoElecIntegrals = []
        
        self.cartesian            = self.getCart()




    def onAtom(self, atomNum, orbType):
        """ 
        returns list of basis functions on atom with number atomNum, and
        that are of type orbType = all 2s2p (all 2-level shells [2s, 2px, 2py, 2pz])
                        1s 2s 2p 2px 2py 2pz (explicity reference, 2p means all 3)
                        s p px py pz         (valence orbs)
        Lists are based on FuncToCenter(which has a seperate entry for each orbital)
=======

def adjust_matlab_input(mat, subtract_one = False):
    '''
    The input must be a vector, either by having one dimension or by being
    a 1xN or Nx1 vector. In any case, the returned value has only one
    dimension.
 
    If subtract_one, then 1 is subtracted from all elements, to account for
    atom referencing starting at 0 instead of 1
    '''
    shp = mat.shape
    non_ones = sum(1 for dim in shp if dim > 1)
    if (non_ones != 1):
        raise Exception('adjust_matlab_input should only be passed vectors')

    if (len(shp) == 2):
        res = mat.reshape([-1])
    else:
        res = mat
        
    if subtract_one:
        res = [x-1 for x in res]
        
    return np.asarray(res)

class Geometry(object):
    """
    A molecule with a specific geometry, as represented by a specific level of 
    theory. Holds the geometry of all molecules.
    """

    def __init__(self, zdata, basisSet, hdata, store_two_elec = True ):
        """
        Create Geometry from matlab input

        Required input:
        zdata -- zmatrix from matlab
        basisSet -- basis set from matlab
        hdata -- Hamiltonian data from matlab

        Keyword arguments:
        storeTwoElec -- (True) Can opt to not store two electron integrals in
                        order to save space
        """
        adj = adjust_matlab_input        
        
        self.angles                = adj(zdata['angles']                     )
        self.ang_ref               = adj(zdata['ang_ref']               ,True)
        self.bond_ref              = adj(zdata['bond_ref']              ,True)
        # Not clear what this is or why we need it
        #self.bond_total            = adj(zdata['bond_total']
        self.bonds                 = adj(zdata['bonds']                      )
        self.di_ref                = adj(zdata['di_ref']                ,True)
        self.dihedrals             = adj(zdata['dihedrals']                  )
        self.nums                  = adj(zdata['nums']                  ,True)
        self.z                     = adj(zdata['z'])

        self.FuncToAngular          = adj(basisSet['FuncToAngular']          )
        self.FuncToCenter           = adj(basisSet['FuncToCenter']      ,True)
        self.FuncToShell            = adj(basisSet['FuncToShell']       ,True)
        self.IsSpherical            =     basisSet['IsSpherical'][0][0]
        self.Name                   =     basisSet['Name']             
        self.NumFunctions           =     basisSet['NumFunctions'][0][0]
        self.NumShells              =     basisSet['NumShells'][0][0]
        self.PrimCoeffUnnorm        = adj(basisSet['PrimCoeffUnnorm']        )
        self.PrimExp                = adj(basisSet['PrimExp']                )
        self.ShellNumFunctions      = adj(basisSet['ShellNumFunctions']      )
        self.ShellNumPrimitives     = adj(basisSet['ShellNumPrimitives']     )
        self.ShellToCenter          = adj(basisSet['ShellToCenter']     ,True)
        self.ShellTypes             = adj(basisSet['ShellTypes']             )
        
        self.corePotentialMats     = hdata['corePotentialMats']
        self.coulombMat            = hdata['coulombMat']
        self.densvec               = hdata['densvec']
        self.rho                   = (2 * 
                  self.densvec.reshape(self.NumFunctions, self.NumFunctions))
        self.exchangeMat           = hdata['exchangeMat']
        self.hfEnergy              = hdata['hfEnergy']
        self.kineticMat            = hdata['kineticMat']
        self.nucRepEnergy          = hdata['nucRepEnergy']
        self.numElectrons          = hdata['numElectrons']
        self.orbital               = hdata['orbital']
        self.orbitalEnergies       = hdata['orbitalEnergies']
        self.overlapMat            = hdata['overlapMat']
        if (store_two_elec):
            self.twoElecIntegrals      = hdata['twoElecIntegrals']
        else:
            self.twoElecIntegrals = []
            
        self.rcart                 = self.__initialize_rcart()
        self.connections           = self.__initialize_connections()
    

    def __initialize_rcart(self):
        '''
        Return 3xNatom matrix with cartesian coordinates
        '''
        natoms = len(self.z)
        rc = np.zeros([3,natoms])
        # first atom is at zero, so nothing needed
        # second atom is on x axis
        rc[0,1] = self.bonds[0]
        # third atom is in xy plane
        rbond = self.bonds[1]
        theta = self.angles[0]
        bond3 = self.bond_ref[2]
        if (bond3 == 1):
            theta = 180 - theta
        thetad = np.deg2rad(theta)
        rc[0,2] = rc[1,bond3]+(rbond * np.cos(thetad))
        rc[1,2] = rc[2,bond3]+(rbond * np.sin(thetad))
        for iatom in range(3,natoms): #self.nums[0][3:]:
        # basing on code below, so using that notation
            b = self.bonds[iatom-1]
            a = np.deg2rad(self.angles[iatom-2])
            t = np.deg2rad(self.dihedrals[iatom - 3])
            k1 = self.bond_ref[iatom]
            k2 = self.ang_ref[iatom]
            k3 = self.di_ref[iatom]
            # reference vectors (r# is (x#,y#,z#) below)
            # r1 is unit vector: Atheta - Abond
            r1 = rc[:,k2]-rc[:,k1]
            r1 = r1/np.linalg.norm(r1)
            # r2 is Aphi - Atheta projected into perp plane
            r2 = rc[:,k3]-rc[:,k2]
            r2 = r2 - np.dot(r1.T,r2) * r1
            r2 = r2/np.linalg.norm(r2)
                # r3 is cross product of above
            r3 = np.cross(r1.T,r2.T)
            r3 = r3.T/np.linalg.norm(r3.T)
            temp = rc[:,k1] + b*(np.cos(a)*r1 + np.sin(a)*(np.cos(t)*r2 - np.sin(t)*r3))
            rc[:,iatom] = temp.T
            # Gaussian convention is apparently different, this fixes it
        np.rad2deg(rc)
        return rc[[1,2,0],:]
        
    def __calc_isbonded(self):
        ## from: https://en.wikipedia.org/wiki/Covalent_radius
        covalent_radii = {1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 9:0.57}
        natom = self.natom()
        res = np.zeros([natom, natom])
        for iatom in range(natom):
            for jatom in range(iatom+1, natom):
                rdist = self.dist(iatom,jatom)
                bond_cutoff = ( covalent_radii[self.z[iatom] ]
                              + covalent_radii[self.z[jatom] ]) + 0.35
                if (rdist < bond_cutoff):
                    res[iatom,jatom] = 1
                    res[jatom,iatom] = 1
        return res

    def __initialize_connections(self):
        natom = self.natom()
        res = np.zeros([natom,natom])
        idiag = np.diag_indices_from(res)
        
        bonded = self.__calc_isbonded()
        conn = np.identity(natom)
        for ilevel in range(1,10):
            #extend connections by 1 level
            upConn = np.dot(conn, bonded)
            # set diagonal to zero
            upConn[idiag] = 0
            # set all elements that were already specified to zero
            upConn[res>0] = 0
            res[upConn > 0] = ilevel
            conn = upConn
        return res        
        

    def natom(self):
        '''
        number of atoms in this geometry
        '''
        return len(self.z)
        
    def on_atom(self, atomNum, orbType):
        """        
        Return list of basis functions (atomic orbitals) of a given type on a 
        given atom
        
        Required input:
        atomNum -- atom number with convention that first atom is number 0
        orbType -- string specifying 
                   "1s" "2s" "2px" "2py" "2pz" -- individual orbital
                   "2p"    -- all p orbitals
                   "2s2p"  -- 2s and all 2 p orbitals
                   "all"   -- all orbitals

>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        """
        
        basisList = np.argwhere(self.FuncToCenter == atomNum)
        if basisList.size == 0:    
<<<<<<< HEAD
            print "Error, atom not found"
            return        
            
        basisList = basisList[0:basisList.shape[0], 1]   # argwhere returns a matrix like:
                                                         #   [[0,0],[0,1],[0,2],[0,etc...]]
                                                         #   we only want the [0,1,2,etc...]   
        if basisList.size == 1:
            return basisList
=======
            raise Exception("Geometry.on_atom:  atom not found")
            
        if basisList.size == 1:
            if (not (orbType  in ['1s', 'all']) ):
                raise Exception("Geometry.on_atom: n=2 orb requested on Hydrogen")
            return np.asarray(basisList).reshape(-1)
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        else:
                
            options ={
                "all": basisList,
                "2s2p": basisList[1:5],
<<<<<<< HEAD
                "1s": np.array([basisList[0]]),
                "2s": np.array([basisList[1]]),
                "2p": basisList[2:5],
                "2px": np.array([basisList[2]]), 
                "2py": np.array([basisList[3]]),
                "2pz": np.array([basisList[4]]),
            }
        
            
        return options[orbType]
        
        
    
    def getSubMatrix(self, oper, list1, list2):
        """ 
=======
                "1s": basisList[0],
                "2s": basisList[1],
                "2p": basisList[2:5],
                "2px": basisList[2], 
                "2py": basisList[3],
                "2pz": basisList[4],
            }
        
        t1 = np.asarray(options[orbType])            
        return t1.reshape(-1)
                
    def rbond(self,atom1, atom2):
        '''
        cartesian vector pointing from atom1 to atom2 
        '''
        return (self.rcart[:,atom2] - self.rcart[:,atom1]).reshape(-1)
        
    def dist(self, atom1, atom2):
        '''
        distance between atom2 and atom1
        '''
        return np.linalg.norm(self.rbond(atom1, atom2))
    
    def submatrix(self, oper, orbsLeft, orbsRight):
        """ 
        Get submatrix of a specified opera                
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        for operator oper, get subMatrix (:list1, :list2)
        oper can be: 'KE' 'H1' 'H1nuc' 'J' 'K' 'H1nuc1' 'H1nuc2' ... 
        where H1nuc# means only the interaction with nucleus #
        Call onAtom for each atom in question and pass it's
<<<<<<< HEAD
        list to getSubMatrix
        """  
        if list1.shape[0] > list2.shape[0]:
            temp = list1
            list1 = list2
            list2 = temp
            
        Mfull = None
        
#        if list1 == None or list2 == None:
#            print "Error, at least one list does not exist"
#            return
        
        # loops across all H1nuc matrix slices to see if any were specified in argument
        for i in range(1, self.corePotentialMats.shape[2] + 1):
            if oper == 'H1nuc%d' % i:
=======
        list to getSubMatrixtor

        Required input:
        oper -- string specifying operator type
                "H1" "KE" "J" "K" "S"  -- standard operators
                "MO"  -- molecular orbitals
                "rho"  -- density matrix 
                "F"   -- fock matrix
                "Hnuc" -- sum of all e-nuc operators
                "Hnuc5" -- e-nuc operator for nucleus of atom 5
                           (atom number convention is first atom is 0)
        """  
    
        Mfull = None


        for i in range(1, self.corePotentialMats.shape[2] + 1):
            if oper == 'Hnuc%d' % i:
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
                Mfull = self.corePotentialMats[:,:,i]
                break
        
        operators = {
            "H1": self.corePotentialMats.sum(2) + self.kineticMat,
<<<<<<< HEAD
            "H1nuc": self.corePotentialMats.sum(2),
=======
            "Hnuc": self.corePotentialMats.sum(2),
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
            "KE": self.kineticMat,
            "J": self.coulombMat,
            "K": self.exchangeMat,
            "S": self.overlapMat,
            "c": self.orbital,
            "rho": self.rho
        }
        
        Mfull = operators[oper]
        # Create a new submatrix based on the indices found. 
        #   Will include cases for the 4 dimensional matrix later.
<<<<<<< HEAD
    
        subMatrix = Mfull[list1[0]:list1[-1] + 1, list2[0]:list2[-1] + 1]
        
        return subMatrix
    
    
    

    def getSubList(self, oper, list1):
        '''
        getSubList finds the portion of a list that is on a particular atom and orbital type
        specified in the list. Call OnAtom beforehand to generate list passed to argument
        '''
       
        Lfull = None
        
        if list1 == None:
            print "Error, list does not exist"
            return
            
        options = {
            'e': self.orbitalEnergies,
            "FuncToCenter": self.FuncToCenter,
            "FuncToAngular": self.FuncToAngular,
            "FuncToShell": self.FuncToShell,
        }
        
        Lfull = options(oper)
        
                                  
        # create a new sub list based on the list provided. Special cases for inconsistent data types.
        if oper == 'e':
            subList = Lfull[list1[0]:list1[-1] + 1]
        else:
            subList = Lfull[0, list1[0]:list1[-1] + 1]

        return subList
        
    def printMatrix(self, oper, printType):
        np.set_printoptions(precision = 6, suppress = True)
        
        options = {
            "TEI": self.twoElecIntegrals,
            "H1": self.corePotentialMats.sum(2) + self.kineticMat,
            "H1nuc": self.corePotentialMats.sum(2),
            "KE": self.kineticMat,
            "J": self.coulombMat,
            "K": self.exchangeMat,
            "S": self.overlapMat,
            "c": self.orbital,
            "rho": self.rho,
            
            'angles': self.angles,
            'ang_ref': self.ang_ref,
            'bond_ref': self.bond_ref,
            'bond_total': self.bond_total,
            'bonds': self.bonds,
            'di_ref': self.di_ref,
            'dihedrals': self.dihedrals,
            'nums': self.nums,
            'z': self.z 
        }
        
        matrix = options[oper]
        
        if printType == "file":
            import csv
            if oper == "TEI":
                with open("Data/%s.csv" % oper, "wb") as f:
                    writer = csv.writer(f)
        
                    for slice_2d in matrix:
#                        writer.writerows(oper)
                        writer.writerows(slice_2d)
            else:
     
                with open("Data/%s.csv" % oper, "wb") as f:
                    writer = csv.writer(f)
#                    writer.writerows(oper)
                    writer.writerows(matrix)
        elif printType == "screen":
            print oper, "\n", matrix, "\n"
            
        return
        
    def getCart(self):
        natoms = len(self.nums[0])
        rc = np.zeros([3,natoms])
        # first atom is at zero, so nothing needed
        # second atom is on x axis
        rc[0,1] = self.bonds[0][0]
        # third atom is in xy plane
        rbond = self.bonds[0][1]
        theta = self.angles[0][0]
        bond3 = self.bond_ref[2]
        if (bond3 == 2):
            theta = 180 - theta
        thetad = np.deg2rad(theta)
        rc[0,2] = rc[1,bond3-1]+(rbond * np.cos(thetad))
        rc[1,2] = rc[2,bond3-1]+(rbond * np.sin(thetad))
        for iatom in range(3,len(self.z)): #self.nums[0][3:]:
        # basing on code below, so using that notation
            b = self.bonds[0][iatom-1]
            a = np.deg2rad(self.angles[0][iatom-2])
            t = np.deg2rad(self.dihedrals[0][iatom - 3])
            k1 = self.bond_ref[iatom]-1
            k2 = self.ang_ref[iatom]-1
            k3 = self.di_ref[iatom]-1
            # reference vectors (r# is (x#,y#,z#) below)
            # r1 is unit vector: Atheta - Abond
            r1 = rc[:,k2]-rc[:,k1]
            r1 = r1/np.linalg.norm(r1)
            # r2 is Aphi - Atheta projected into perp plane
            r2 = rc[:,k3]-rc[:,k2]
            r2 = r2 - np.dot(r1.T,r2) * r1
            r2 = r2/np.linalg.norm(r2)
                # r3 is cross product of above
            r3 = np.cross(r1.T,r2.T)
            r3 = r3.T/np.linalg.norm(r3.T)
            temp = rc[:,k1] + b*(np.cos(a)*r1 + np.sin(a)*(np.cos(t)*r2 - np.sin(t)*r3))
            rc[:,iatom] = temp.T
            # Gaussian convention is apparently different, this fixes it
        np.rad2deg(rc)
        return rc[[1,2,0],:]
=======
        return Mfull[np.ix_(orbsLeft,orbsRight)]
        

    def connection_order(self, atom1, atom2):
        '''
        Get minimum number of bonds connecting atom1 to atom2. For example,
        if atom1==atom2, return 0
        if atoms are bonded, return 1
        if atoms are bonded to a common atom, return 2
        The data is precomputed to make evaluation lowcost
        
        In molecular mechanics, nonbonded atoms have a connectivity >=3
        
        Matrix is directly accessible as Geometry.connections        
        
        Require input:
        atom1 -- number of first atom 
        atom2 -- number of second atom
                 (with convention that first atom is 0)
        '''
        return self.connections[atom1,atom2]

    def substituent_level(self, atom_number):
        '''
        Number of heavy (i.e. non-hydrogen) atoms bonded to atom_number
        '''
        return 0
        
    def get_atoms(self, z):
        '''
        Returns list of atom numbers for element type Z
        '''
        return [i for i, x in enumerate(self.z) if x == z]
    
    def get_bonds(self, z1, z2, conn_order_target):
        '''
        Returns list of dubles (atom_number1, atom_number2) where
        atom1 is of type z1, atom2 is of type z2, and the connection order
        agrees with the provided target values. 
        
        If z1 == z2, returned bonds have atom_number1 < atom_number2
        
        Input:
        z1, z2 : element numbers of first and second atom
        conn_order_target: order is the number of bonds along shortest path 
                           between two atoms. Return only those with this
                           value.
        '''
        res = []
        for a1 in self.get_atoms(z1):
            for a2 in self.get_atoms(z2):
                if (not ((z1 == z2) & (a1 >= a2))):                    
                    if self.connection_order(a1,a2) == conn_order_target:
                        res.append( (a1, a2) )
        return res
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
