# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 18:49:03 2015

@author: brandon
"""


import os
import cPickle as pickle
from collections import OrderedDict

import numpy as np

from MoleculeSet import MoleculeSet
from miscFunctions import convertToKcalPerMol
import generateFiles as gen

def generateAllExpectation_data(dataPath, dataType, relType, orbType, oper, atoms, molType, rhoList, overwriteMats = False):
    '''
    generateAllExpectation_data will generate a dictionary of norms. If a dictionary of norms is already on file, this function will
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
        LL_exp = calcExpectiation(LL_mat, rhoList)
        fullPath = os.path.join(dataPath, "subMatrices", "QU_%s_%s_subMatrices" % (atoms, oper))
        QU_mat = pickle.load(open(fullPath, 'rb'))
        QU_exp = calcExpectiation(QU_mat, rhoList)
    else:
        LL_mat, QU_mat = generateAllSubMatrix_data(fullPath, relType, orbType, oper, atoms, molType)
        LL_exp = calcExpectiation(LL_mat, rhoList)
        QU_exp = calcExpectiation(QU_mat, rhoList)
       
    fullPath = os.path.join(dataPath, "%ss" % dataType) 
    LL_fileName = "LL_%s_%s_%ss" % (atoms, oper, dataType)
    QU_fileName = "QU_%s_%s_%ss" % (atoms, oper, dataType)
    gen.saveDataToFolder_pickle(fullPath, LL_fileName, LL_exp)
    gen.saveDataToFolder_pickle(fullPath, QU_fileName, QU_exp) 

    return LL_exp, QU_exp

    
def generateRho(dataPath, dataType, relType, orbType, atoms, molType, overwriteMats = False):
    r = os.path.isfile(dataPath + "/dataFiles/%ss/%s_rho_subMatrices" % (dataType, atoms))
    
    if (r) and not overwriteMats:
        
        fullPath = os.path.join(dataPath, "/dataFiles/%ss/%s_rho_subMatrices" % (dataType, atoms))
        rhoList = pickle.load(open(fullPath))
    else:
         
        if relType == "Bonded":
            rhoList = MoleculeSet("LL", storeTwoElec = False).getAllSubMatricesBonded(atoms, molType, orbType, "rho")
        
        elif relType == "Nonbonded":
            rhoList = MoleculeSet("LL", storeTwoElec = False).getAllSubMatricesNonbonded(atoms, molType, orbType, "rho")
            
        elif relType == "Self":
            rhoList = MoleculeSet("LL", storeTwoElec = False).getAllSubMatricesSelf(atoms, molType, orbType, "rho")
            
        gen.saveDataToFolder_pickle(dataPath, "%s_rho_subMatrices" % atoms, rhoList)
    
    return rhoList
    
    
    
def generateAllSubMatrix_data(dataPath, relType, orbType, oper, atoms, molType):
    '''  
    generateAllSubMatrix_data is called by generateAllExpectation_data and generates files of subMatrices for the for the 
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
    elif relType == "Self":
        LL_subMatrices = MoleculeSet("LL", storeTwoElec = False).getAllSubMatricesSelf(atoms, molType, orbType, oper)
        QU_subMatrices = MoleculeSet("QU", storeTwoElec = False).getAllSubMatricesSelf(atoms, molType, orbType, oper) 
        
    gen.saveDataToFolder_pickle(dataPath, "LL_%s_%s_subMatrices" % (atoms, oper), LL_subMatrices)
    gen.saveDataToFolder_pickle(dataPath, "QU_%s_%s_subMatrices" % (atoms, oper), QU_subMatrices)

    return LL_subMatrices, QU_subMatrices
    
    

def checkFilesExp(dataPath, dataType, relType, oper, atoms, orbType, molType, 
               rhoList, overwriteExp = False, overwriteMats = False):
    '''   
    checkFiles is called in calcExpectiationBonded_save and checks to see if there is already data on file that can be used
    to save time. If there is no data, or if overwriteNorms is True, it will then call generateAllNorms_data to create it.
    '''
    LL = os.path.isfile(dataPath + "/dataFiles/" + dataType + "s/" + "LL_" + atoms + "_" + oper + "_exp")
    QU = os.path.isfile(dataPath + "/dataFiles/" + dataType + "s/" + "QU_" + atoms + "_" + oper + "_exp")
    if (LL and QU) and not (overwriteExp or overwriteMats):
        
        LL_fullPath = os.path.join(dataPath, "dataFiles",  "%ss" % dataType, "LL_%s_%s_exp" % (atoms, oper))
        LL_exp = pickle.load(open(LL_fullPath))
        
        QU_fullPath = os.path.join(dataPath, "dataFiles",  "%ss" % dataType, "QU_%s_%s_exp" % (atoms, oper))
        QU_exp = pickle.load(open(QU_fullPath))
    else:
         LL_exp, QU_exp = generateAllExpectation_data(dataPath, dataType, relType, orbType, oper, atoms, molType, 
                                                    rhoList, overwriteMats)
         
    return LL_exp, QU_exp


def getExpectiation_save(relType, atoms, molType, overwriteExp = False, overwriteMats = False):
    '''
    calcExpectiation takes all the norms for a particular bond type ("CC", "XC", "CH", or "XH"). It then creates
    a set of folders, files, plots, and a markdown page with a display of the data.
    '''
        # these strings used to generate text in folders, pages, tables, and plots
    if atoms == "CH" and relType == "Bonded":    
        molType = ["C_unsub", "C_sub", "C_disub"]

    dataType = "expectation"
    orbType = '2s2p'
    string2 = '## %s(%s) %s orbitals\n\n' % (atoms, relType, orbType)
    xLab, yLab = "LL Submatrix", "QU Submatrix"
    dataPath = os.path.join("Data" , "%sPlots" % dataType , relType , atoms)
    
    rhoList = generateRho(dataPath, dataType, relType, orbType, atoms, molType, overwriteMats)
        
    for group in ['K', 'H1nuc', 'KE', 'S', 'J']: 
                            
        units = "unitless" 
        title = "LL vs QU: %s(%s) %s " % (atoms, relType, group)
        title += "Expectation Vals Across All Sets\n(%s orbitals)" % orbType
        
        # create norms if they don't already exist on file
        LL_exp, QU_exp = checkFilesExp(dataPath, dataType, relType, group, atoms, orbType, molType, rhoList,
                                        overwriteExp, overwriteMats)
        
        # all matrices except overlapMat get converted to kcal/mol
        if group != 'S':
            units = "kcal/mol" 
            LL_exp = convertToKcalPerMol(LL_exp)  
            QU_exp = convertToKcalPerMol(QU_exp)  
         
        # generate data table and plot based on norms
        string = gen.composeDataTable(LL_exp, QU_exp, title, xLab, yLab, units)
        
        # add table and plot to markdown page
        string2 += "<p align=\"center\">"
        string2 += "<img src=/Data/%sPlots/%s/%s/plots/%s_%s.png /></p>\n\n" % (dataType, relType, atoms, atoms, group) + string
        gen.savePlotToFolder(os.path.join(dataPath, "plots"), "%s_%s" % (atoms, group))
        
    
    # save markdown page
    gen.saveTextToFolder(dataPath, "%s(%s)" % (atoms, relType), string2)  
    
    return
    
def getExpectiationSelf_save(atom, molType, overwriteExp, overwriteMats):
    '''
    calcExpectiation takes all the norms for a particular bond type ("CC", "XC", "CH", or "XH"). It then creates
    a set of folders, files, plots, and a markdown page with a display of the data.
    '''
    # these strings used to generate text in folders, pages, tables, and plots
    relType = "Self"
    units = "kcal/mol" 
    dataType = "expectation"
    orbType = '2s2p'
    string2 = '## %s(%s) %s orbitals\n\n' % (atom, relType, orbType)
    xLab, yLab = "LL Submatrix", "QU Submatrix"
    dataPath = os.path.join("Data" , "%sPlots" % dataType , relType , atom)
    
    rhoList = generateRho(dataPath, dataType, relType, orbType, atom, molType, overwriteMats)
    
    LL_matrices = {
        'KE': [],
        'H1nuc' : [],
        'J' : [],
        'K' : [],
    }
    
    QU_matrices = {
        'KE': [],
        'H1nuc' : [],
        'J' : [],
        'K' : [],
    }
    for group1 in LL_matrices.keys(): 
                            
       
        # create norms if they don't already exist on file
        LL_matrices[group1], QU_matrices[group1] = checkFilesExp(dataPath, dataType, relType, group1, atom, orbType, 
                                                            molType, rhoList, overwriteExp, overwriteMats)
        
        if group1 != 'rho':
            LL_matrices[group1] = convertToKcalPerMol(LL_matrices[group1])  
            QU_matrices[group1] = convertToKcalPerMol(QU_matrices[group1])  
         

         
    LL_matrices["F"] = LL_matrices["KE"]
    QU_matrices["F"] = QU_matrices["KE"]
    
    for i, n in enumerate(LL_matrices["KE"]):
        for j, m in enumerate(n):
            LL_matrices["F"][n][j] = LL_matrices["KE"][n][j] + LL_matrices["H1nuc"][n][j] + 2 * LL_matrices["J"][n][j] - LL_matrices["K"][n][j]
            QU_matrices["F"][n][j] = QU_matrices["KE"][n][j] + QU_matrices["H1nuc"][n][j] + 2 * QU_matrices["J"][n][j] - QU_matrices["K"][n][j]
    title = "LL vs QU: %s(%s) %s " % (atom, relType, "F")
    title += "Expectation Vals Across All Sets\n(%s orbitals)" % orbType
    string = gen.composeDataTable( LL_matrices["F"], QU_matrices["F"], title, xLab, yLab, units)
        
        # add table and plot to markdown page
    string2 += "<p align=\"center\">"
    string2 += "<img src=/Data/%sPlots/%s/%s/plots/%s_%s.png /></p>\n\n" % (dataType, relType, atom, atom, "F") + string
    

    gen.savePlotToFolder(os.path.join(dataPath, "plots"), "%s_%s" % (atom, "F"))
        
    
    # save markdown page
    gen.saveTextToFolder(dataPath, "%s(%s)" % (atom, relType), string2)  
    
    return
    
    
    
def getExpectiationSelfSuperset_save(atom, molType, overwriteExp = False, overwriteMats = False):
    '''
    calcExpectiation takes all the norms for a particular bond type ("CC", "XC", "CH", or "XH"). It then creates
    a set of folders, files, plots, and a markdown page with a display of the data.
    '''
    # these strings used to generate text in folders, pages, tables, and plots
    relType = "Self"
    units = "kcal/mol" 
    dataType = "expectation"
    orbType = '2s2p'
    LL_exp = OrderedDict()        
    QU_exp = OrderedDict()
    for atom in ['H', 'C', 'N', 'O', 'F']:
        string2 = '## %s(%s) %s orbitals\n\n' % (atom, relType, orbType)
        xLab, yLab = "LL Submatrix", "QU Submatrix"
        dataPath = os.path.join("Data" , "%sPlots" % dataType , relType , atom)
        
    
        LL1 = OrderedDict()
        QU1 = OrderedDict()

        for group1 in ["KE", "H1nuc", "J", "K"]:
#            LL1 = {
#            'KE': [],
#            'H1nuc' : [],
#            'J' : [],
#            'K' : [],
#            }
#            
#            QU1 = {
#                'KE': [],
#                'H1nuc' : [],
#                'J' : [],
#                'K' : [],
#            }                    
            dataPath = os.path.join("Data" , "%sPlots" % dataType , relType , relType, atom)
            rhoList = generateRho(dataPath, dataType, relType, orbType, atom, molType, overwriteMats)
            # create norms if they don't already exist on file
            LL1[group1], QU1[group1] = checkFilesExp(dataPath, dataType, relType, group1, atom, orbType, 
                                                                molType, rhoList, overwriteExp, overwriteMats)

            if group1 != 'rho':
                LL1[group1] = convertToKcalPerMol(LL1[group1])  
                QU1[group1] = convertToKcalPerMol(QU1[group1])  
             
    
             
        LL1["F"] = LL1["KE"]
        QU1["F"] = QU1["KE"]
        
        for i, n in enumerate(LL1["KE"]):
            for j, m in enumerate(n):
                LL1["F"][n][j] = LL1["KE"][n][j] + LL1["H1nuc"][n][j] + 2 * LL1["J"][n][j] - LL1["K"][n][j]
                QU1["F"][n][j] = QU1["KE"][n][j] + QU1["H1nuc"][n][j] + 2 * QU1["J"][n][j] - QU1["K"][n][j]
                
        keys = LL1["F"].keys()
        for key in keys:
                LL_exp[key] = LL1["F"][key]  
                    
        keys = QU1["F"].keys()
        for key in keys:
                QU_exp[key] = QU1["F"][key]  
                        
    title = "LL vs QU: (%s) %s " % (relType, "F")
    title += "Expectation Vals Across All Sets\n(%s orbitals)" % orbType
    string = gen.composeDataTable(LL_exp, QU_exp, title, xLab, yLab, units)
        
        # add table and plot to markdown page
    string2 += "<p align=\"center\">"
    string2 += "<img src=/Data/%sPlots/%s/%s/plots/%s_%s.png /></p>\n\n" % (dataType, relType, relType, relType, "F") + string
    
    dataPath = os.path.join("Data" , "%sPlots" % dataType , relType , relType)
    gen.savePlotToFolder(os.path.join(dataPath, "plots"), "%s_%s" % (relType, "F"))
            
    
    # save markdown page
    gen.saveTextToFolder(dataPath, "%s(%s)" % (relType, relType), string2)  
    
    return   
def getAllExpectationSuperset_save(relType, molType, overwriteExp = False, overwriteMats = False): 
    '''
    getAllNormsNonbonded_save takes all the norms for a particular pair type. It then creates
    a set of folders, files, plots, and a markdown page with a display of the data.
    '''
        # these strings used to generate text in folders, pages, tables, and plots
    LL_exp = OrderedDict()
    QU_exp = OrderedDict()
    dataType = "expectation"
    orbType = '2s2p'
    string2 = '## (%s) %s orbitals\n\n' % (relType, orbType)
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
        LL_exp = OrderedDict()
        QU_exp = OrderedDict()                    

        units = "unitless" 
        title = "LL vs QU: (%s) %s " % (relType, group)
        title += "Expectation Vals Across All Sets\n(%s orbitals)" % orbType
        # create norms if they don't already exist on file
            
        if relType == "Both":     
            for pairs in pairList1:
                dataPath = os.path.join("Data" , "%sPlots" % dataType , "Bonded" , "Bonded", pairs)
                rhoList = generateRho(dataPath, dataType, "Bonded", orbType, pairs, molType, overwriteMats)

                LL1, QU1 = checkFilesExp(dataPath, dataType, "Bonded", group, pairs, orbType, molType, rhoList,
                                                overwriteExp, overwriteMats)

                if group != 'S':
                    units = "kcal/mol" 
                    LL1 = convertToKcalPerMol(LL1)  
                    QU1 = convertToKcalPerMol(QU1)  
 
                    
                keys = LL1.keys()
                for key in keys:
                    LL_exp[key] = LL1[key]  
                    
                keys = QU1.keys()
                for key in keys:
                    QU_exp[key] = QU1[key]  
            for pairs in pairList2: 

                dataPath = os.path.join("Data" , "%sPlots" % dataType , "Nonbonded" , "Nonbonded", pairs) 
                rhoList = generateRho(dataPath, dataType, "Nonbonded", orbType, pairs, molType, overwriteMats)                               
                LL2, QU2 = checkFilesExp(dataPath, dataType, "Nonbonded", group, pairs, orbType, molType, rhoList,
                                      overwriteExp, overwriteMats)    
                                      
                if group != 'S':
                    units = "kcal/mol" 
                    LL2 = convertToKcalPerMol(LL2)  
                    QU2 = convertToKcalPerMol(QU2) 
                    
                keys = LL2.keys()
                for key in keys:
                    LL_exp[key] = LL2[key]  
                    
                keys = QU2.keys()
                for key in keys:
                    QU_exp[key] = QU2[key]   
        elif relType == "Bonded" or relType == "Nonbonded":  
            for pairs in pairList:
                dataPath = os.path.join("Data" , "%sPlots" % dataType , relType , relType, pairs)
    
                rhoList = generateRho(dataPath, dataType, relType, orbType, pairs, molType, overwriteMats)
#                units = "unitless" 
#                title = "LL vs QU: %s(%s) %s " % (pairs, relType, group)
#                title += "Expectation Vals Across All Sets\n(%s orbitals)" % orbType
        
                # create norms if they don't already exist on file
                LL1, QU1 = checkFilesExp(dataPath, dataType, relType, group, pairs, orbType, molType, rhoList,
                                        overwriteExp, overwriteMats)
                                                
                if group != 'S':
                    units = "kcal/mol" 
                    LL1 = convertToKcalPerMol(LL1)  
                    QU1 = convertToKcalPerMol(QU1)    
                    
                keys = LL1.keys()
                for key in keys:
                    LL_exp[key] = LL1[key]  
                    
                keys = QU1.keys()
                for key in keys:
                    QU_exp[key] = QU1[key]  

            
        # generate data table and plot based on norms
        string = gen.composeDataTable(LL_exp, QU_exp, title, xLab, yLab, units)
        
        # add table and plot to markdown page
        string2 += "<p align=\"center\">"
        string2 += "<img src=/Data/%sPlots/%s/%s/plots/%s_%s.png /></p>\n\n" % (dataType, relType, relType, relType, group) + string 
        gen.savePlotToFolder(os.path.join(newDataPath, "plots"), "%s_%s" % (relType, group))
        
    
    # save markdown page
    gen.saveTextToFolder(newDataPath, "%s" % relType, string2)  
    
    return
        
def calcExpectiation(matrixList, rhoList):
    expectation = OrderedDict()

    for key in rhoList.keys():
        entryList = []
        for i, n in enumerate(rhoList[key]):
            entry = rhoList[key][i] * matrixList[key][i]
            entry = entry.sum()
            entryList.append(entry)
        expectation[key] = entryList
        
    return expectation
    
if __name__ == "__main__":
#    import sys
#    getAllNormsBonds_save(sys.argv[1], sys.argv[2] = False, sys.argv[3] = False)
    overwrite = True
    overwriteExp = overwrite
    overwriteMats = overwrite
    molType = [0, 1, 2, 3]
#    getAllExpectationSuperset_save("Bonded", ["all"])
#    getAllExpectationSuperset_save("Nonbonded", ["all"])       
#    getAllExpectationSuperset_save("Both", ["all"])
#    getExpectiationSelfSuperset_save("Self", ["all"])
    #all of these have to generate the same lists 30 times, can be made faster if only generate once, then make matrices
#    
#    bonds = "CC"
#    getExpectiation_save("Bonded", bonds, molType, overwriteExp, overwriteMats)
#    bonds = "XC"
#    getExpectiation_save("Bonded", bonds, molType, overwriteExp, overwriteMats)
#    bonds = "CH"
#    getExpectiation_save("Bonded", bonds, molType, overwriteExp, overwriteMats)
#    bonds = "XH"
#    getExpectiation_save("Bonded", bonds, molType, overwriteExp, overwriteMats)
    
#    pairs = "XX_1"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
#    pairs = "XX_2"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
#    pairs = "XC_1"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
#    pairs = "XH_1"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
#    pairs = "XH_2"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
#    pairs = "CH_1"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
#    pairs = "CH_2"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
#    pairs = "HH_1"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
#    pairs = "HH_2"
#    getExpectiation_save("Nonbonded", pairs, molType, overwriteExp, overwriteMats)
    
#    atom = "H"
#    getExpectiationSelf_save(atom, molType, overwriteExp, overwriteMats)
#    atom = "C"
#    getExpectiationSelf_save(atom, molType, overwriteExp, overwriteMats)
#    atom = "N"
#    getExpectiationSelf_save(atom, molType, overwriteExp, overwriteMats)
#    atom = "O"
#    getExpectiationSelf_save(atom, molType, overwriteExp, overwriteMats)
#    atom = "F"
#    getExpectiationSelf_save(atom, molType, overwriteExp, overwriteMats)
    print "Done!"