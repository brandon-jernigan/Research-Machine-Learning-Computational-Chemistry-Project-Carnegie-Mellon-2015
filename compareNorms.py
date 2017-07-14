# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 02:20:34 2015

@author: brandon
"""

import os
import cPickle as pickle
from collections import OrderedDict

import numpy as np

from MoleculeSet import MoleculeSet
from miscFunctions import convertToKcalPerMol
import generateFiles as gen









def generateAllNorms_data(dataPath, dataType, relType, orbType, oper, atoms, molType, overwriteMats = False):
    '''
    generateAllNorms_data will generate a dictionary of norms. If a dictionary of norms is already on file, this function will
    be skipped if overWriteNorms is False. If a list of subMatrices is already on file, this function will use them
    to save time if overwriteMats is False. If no data is found, or if the overwrite bools are True this function 
    will generate the datafiles and return the norms dictionaries for the matrix operation (oper), orbital type, and
    bond of interest.
    '''
#    if atoms != "CC" and atoms != "XC" and atoms != "CH" and atoms != "XH":
#        print "Error, incorrect atoms format, options are: CC, XC, CH, or XH"
#        return
    
    dataPath = os.path.join(dataPath, "dataFiles")
    fullPath = os.path.join(dataPath, "subMatrices") 
    LL = os.path.isfile(dataPath + "/subMatrices/LL_%s_%s_subMatrices" % (atoms, oper))
    QU = os.path.isfile(dataPath + "/subMatrices/QU_%s_%s_subMatrices" % (atoms, oper))
    
    if (LL and QU) and not overwriteMats:
        fullPath = os.path.join(dataPath, "subMatrices", "LL_%s_%s_subMatrices" % (atoms, oper))
        LL_mat = pickle.load(open(fullPath, 'rb'))
        LL_norms = getNorms(LL_mat)
        fullPath = os.path.join(dataPath, "subMatrices", "QU_%s_%s_subMatrices" % (atoms, oper))
        QU_mat = pickle.load(open(fullPath, 'rb'))
        QU_norms = getNorms(QU_mat)
    else:
        LL_mat, QU_mat = generateAllSubMatrix_data(fullPath, relType, orbType, oper, atoms, molType)
        LL_norms = getNorms(LL_mat)
        QU_norms = getNorms(QU_mat)
       
    fullPath = os.path.join(dataPath, "%ss" % dataType) 
    LL_fileName = "LL_%s_%s_%ss" % (atoms, oper, dataType)
    QU_fileName = "QU_%s_%s_%ss" % (atoms, oper, dataType)
    gen.saveDataToFolder_pickle(fullPath, LL_fileName, LL_norms)
    gen.saveDataToFolder_pickle(fullPath, QU_fileName, QU_norms) 

    return LL_norms, QU_norms
    

def generateAllSubMatrix_data(dataPath, relType, orbType, oper, atoms, molType):
    '''  
    generateAllSubMatrix_data is called by generateAllNorms_data and generates files of subMatrices for the for the 
    matrix operation (oper), orbital type, and bond of interest.
    '''
#    if atoms != "CC" and atoms != "XC" and atoms != "CH" and atoms != "XH":
#        print "Error, incorrect atoms format, options are: CC, XC, CH, or XH"
#        return

    
    if relType == "Bonded":
        LL_subMatrices = MoleculeSet("LL", storeTwoElec = False).getAllSubMatricesBonded(atoms, molType, orbType, oper)
        QU_subMatrices = MoleculeSet("QU", storeTwoElec = False).getAllSubMatricesBonded(atoms, molType, orbType, oper)
    elif relType == "Nonbonded":
        LL_subMatrices = MoleculeSet("LL", storeTwoElec = False).getAllSubMatricesNonbonded(atoms, molType, orbType, oper)
        QU_subMatrices = MoleculeSet("QU", storeTwoElec = False).getAllSubMatricesNonbonded(atoms, molType, orbType, oper)
        
    gen.saveDataToFolder_pickle(dataPath, "LL_%s_%s_subMatrices" % (atoms, oper), LL_subMatrices)
    gen.saveDataToFolder_pickle(dataPath, "QU_%s_%s_subMatrices" % (atoms, oper), QU_subMatrices)

    return LL_subMatrices, QU_subMatrices
    
    

def checkFiles(dataPath, dataType, relType, oper, atoms, group, orbital, molType, overwriteNorms = False, overwriteMats = False):
    '''   
    checkFiles is called in getAllNormsBonded_save and checks to see if there is already data on file that can be used
    to save time. If there is no data, or if overwriteNorms is True, it will then call generateAllNorms_data to create it.
    '''
    LL = os.path.isfile(dataPath + "/dataFiles/" + dataType + "s/" + "LL_" + atoms + "_" + group + "_norms")
    QU = os.path.isfile(dataPath + "/dataFiles/" + dataType + "s/" + "QU_" + atoms + "_" + group + "_norms")
    if (LL and QU) and not (overwriteNorms or overwriteMats):
        
        LL_fullPath = os.path.join(dataPath, "dataFiles",  "%ss" % dataType, "LL_%s_%s_norms" % (atoms, group))
        LL_norms = pickle.load(open(LL_fullPath))
        
        QU_fullPath = os.path.join(dataPath, "dataFiles",  "%ss" % dataType, "QU_%s_%s_norms" % (atoms, group))
        QU_norms = pickle.load(open(QU_fullPath))
    else:
         LL_norms, QU_norms = generateAllNorms_data(dataPath, dataType, relType, orbital, oper, atoms, molType, overwriteMats)
         
    return LL_norms, QU_norms



def getAllNorms_save(relType, atoms, molType, overwriteNorms = False, overwriteMats = False): 
    '''
    getAllNormsBonded_save takes all the norms for a particular bond type ("CC", "XC", "CH", or "XH"). It then creates
    a set of folders, files, plots, and a markdown page with a display of the data.
    '''
        # these strings used to generate text in folders, pages, tables, and plots
    if atoms == "CH" and relType == "Bonded":    
        molType = ["C_unsub", "C_sub", "C_disub"]

    dataType = "norm"
    orbital = '2s2p'
    string2 = '## %s(%s) %s orbitals\n\n' % (atoms, relType, orbital)
    xLab, yLab = "LL Submatrix", "QU Submatrix"
    dataPath = os.path.join("Data" , "%sPlots" % dataType , relType , atoms)
    
    for group in ['K', 'H1nuc', 'KE', 'S', 'J']: 
                            
        units = "unitless" 
        title = "LL vs QU: %s(%s) %s " % (atoms, relType, group)
        title += "Norms Across All Sets\n(%s orbitals)" % orbital
        
        # create norms if they don't already exist on file
        LL_norms, QU_norms = checkFiles(dataPath, dataType, relType, group, atoms, group, orbital, molType,
                                        overwriteNorms, overwriteMats)
        
        # all matrices except overlapMat get converted to kcal/mol
        if group != 'S':
            units = "kcal/mol" 
            LL_norms = convertToKcalPerMol(LL_norms)  
            QU_norms = convertToKcalPerMol(QU_norms)  
         
        # generate data table and plot based on norms
        string = gen.composeDataTable(LL_norms, QU_norms, title, xLab, yLab, units)
        
        # add table and plot to markdown page
        string2 += "<p align=\"center\">"
        string2 += "<img src=/Data/normPlots/%s/%s/plots/%s_%s.png /></p>\n\n" % (relType, atoms, atoms, group) + string
        gen.savePlotToFolder(os.path.join(dataPath, "plots"), "%s_%s" % (atoms, group))
        
    
    # save markdown page
    gen.saveTextToFolder(dataPath, "%s" % relType, string2)  
    
    return


def getAllNormsSuperset_save(relType, molType, overwriteNorms = False, overwriteMats = False): 
    '''
    getAllNormsNonbonded_save takes all the norms for a particular pair type. It then creates
    a set of folders, files, plots, and a markdown page with a display of the data.
    '''
        # these strings used to generate text in folders, pages, tables, and plots
    LL_norms = OrderedDict()
    QU_norms = OrderedDict()
    dataType = "norm"
    orbital = '2s2p'
    string2 = '## (%s) %s orbitals\n\n' % (relType, orbital)
    xLab, yLab = "LL Submatrix", "QU Submatrix"
    newDataPath = os.path.join("Data" , "%sPlots" % dataType , relType , relType)
    if relType == "Bonded":
        pairList = ["CH", "XH", "XC", "CC"]

    elif relType == "Nonbonded":
        pairList = ["HH_1", "HH_2", "CH_1", "CH_2", "XH_1", "XH_2", "XC_1", "XX_1", "XX_2"]
    elif relType == "Both":
        pairList1 = ["CH", "XH", "XC", "CC"] 
        pairList2 = ["HH_1", "HH_2", "CH_1", "CH_2", "XH_1", "XH_2", "XC_1", "XX_1", "XX_2"]

        
    for group in ['K', 'H1nuc', 'KE', 'S', 'J']: 
        LL_norms = OrderedDict()
        QU_norms = OrderedDict()                    
        units = "unitless" 
        title = "LL vs QU: (%s) %s " % (relType, group)
        title += "Norms Across All Sets\n(%s orbitals)" % orbital
        
        # create norms if they don't already exist on file
            
        if relType == "Both":     
            for pairs in pairList1:
                dataPath = os.path.join("Data" , "%sPlots" % dataType , "Bonded" , "Bonded", pairs)
                LL1, QU1 = checkFiles(dataPath, dataType, "Bonded", group, pairs, group, orbital, molType,
                                                overwriteNorms, overwriteMats)

                if group != 'S':
                    units = "kcal/mol" 
                    LL1 = convertToKcalPerMol(LL1)  
                    QU1 = convertToKcalPerMol(QU1)  
 
                    
                keys = LL1.keys()
                for key in keys:
                    LL_norms[key] = LL1[key]  
                    
                keys = QU1.keys()
                for key in keys:
                    QU_norms[key] = QU1[key]  
            for pairs in pairList2: 

                dataPath = os.path.join("Data" , "%sPlots" % dataType , "Nonbonded" , "Nonbonded", pairs)                                
                LL2, QU2 = checkFiles(dataPath, dataType, "Nonbonded", group, pairs, group, orbital, molType,
                                      overwriteNorms, overwriteMats)    
                                      
                if group != 'S':
                    units = "kcal/mol" 
                    LL2 = convertToKcalPerMol(LL2)  
                    QU2 = convertToKcalPerMol(QU2) 
                    
                keys = LL2.keys()
                for key in keys:
                    LL_norms[key] = LL2[key]  
                    
                keys = QU2.keys()
                for key in keys:
                    QU_norms[key] = QU2[key]   
        else:
            for pairs in pairList:
                dataPath = os.path.join("Data" , "%sPlots" % dataType , relType , relType, pairs)
                LL1, QU1 = checkFiles(dataPath, dataType, relType, group, pairs, group, orbital, molType,
                                                overwriteNorms, overwriteMats)
                                                
                if group != 'S':
                    units = "kcal/mol" 
                    LL1 = convertToKcalPerMol(LL1)  
                    QU1 = convertToKcalPerMol(QU1)  
                    
                keys = LL1.keys()
                for key in keys:
                    LL_norms[key] = LL1[key]  
                    
                keys = QU1.keys()
                for key in keys:
                    QU_norms[key] = QU1[key]  
            
        # generate data table and plot based on norms
        string = gen.composeDataTable(LL_norms, QU_norms, title, xLab, yLab, units)
        
        # add table and plot to markdown page
        string2 += "<p align=\"center\">"
        string2 += "<img src=/Data/normPlots/%s/%s/plots/%s_%s.png /></p>\n\n" % (relType, relType, relType, group) + string 
        gen.savePlotToFolder(os.path.join(newDataPath, "plots"), "%s_%s" % (relType, group))
        
    
    # save markdown page
    gen.saveTextToFolder(newDataPath, "%s" % relType, string2)  
    
    return
    
def getNorms(matrixList):
    '''
    getNorms takes a list, OrderedDict or list of matrices and converts into thier norms. Can also take a single matrix
    '''
    temp = []

    if type(matrixList) == dict or type(matrixList) == OrderedDict:
        norms = OrderedDict() 
        for group in matrixList.keys():
            for i in range(0, len(matrixList[group])):
            
                if matrixList[group][i].size > 1:
                    temp.append(np.linalg.norm(matrixList[group][i]))
                else:
                    temp.append(matrixList[group][i][0])    # a 1x1 matrix is its own norm
            norms[group] = temp
            temp = []
    elif type(matrixList) == list:
        norms = []
        for i in range(0, len(matrixList)):
        
            if matrixList[i].size > 1:
                temp.append(np.linalg.norm(matrixList[i]))
            else:
                temp.append(matrixList[i][0])    # a 1x1 matrix is its own norm
            norms += [temp]
            temp = []
    elif type(matrixList) == np.ndarray:
        if matrixList.size > 1:
            norms = np.linalg.norm(matrixList)
        else:
            norms = matrixList[0]  # a 1x1 matrix is its own norm

    return norms
    
#%%

if __name__ == "__main__":
#    import sys
#    getAllNormsBonds_save(sys.argv[1], sys.argv[2] = False, sys.argv[3] = False)
#    overwrite = True
#    overwriteNorms = overwrite
#    overwriteMats = overwrite
#    molType = [0, 1, 2, 3]
#    
#    getAllNormsSuperset_save("Bonded", ["all"], overwriteNorms, overwriteMats)
#    getAllNormsSuperset_save("Nonbonded", ["all"])        
#    getAllNormsSuperset_save("Both", ["all"], overwriteNorms, overwriteMats)
    #all of these have to generate the same lists 30 times, can be made faster if only generate once, then make matrices
#    
#    bonds = "CC"
#    getAllNorms_save("Bonded", bonds, molType, overwriteNorms, overwriteMats)
#    bonds = "XC"
#    getAllNorms_save("Bonded", bonds, molType, overwriteNorms, overwriteMats)
#    bonds = "CH"
#    getAllNorms_save("Bonded", bonds, molType, overwriteNorms, overwriteMats)
#    bonds = "XH"
#    getAllNorms_save("Bonded", bonds, molType, overwriteNorms, overwriteMats)
#    pairs = "XX_1"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)
#    pairs = "XX_2"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)
#    pairs = "XC_1"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)
#    pairs = "XH_1"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)
#    pairs = "XH_2"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)
#    pairs = "CH_1"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)
#    pairs = "CH_2"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)
#    pairs = "HH_1"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)
#    pairs = "HH_2"
#    getAllNorms_save("Nonbonded", pairs, molType, overwriteNorms, overwriteMats)

    print "Done!"