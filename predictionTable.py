# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 00:22:06 2015

@author: brandon
"""
from collections import OrderedDict
import itertools
import os
<<<<<<< HEAD
#from operator import sub

=======
from operator import sub

import math
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
import numpy as np
import csv
import matplotlib.pyplot as plt

<<<<<<< HEAD
from MoleculeSet import MoleculeSet
from miscFunctions import convertToKcalPerMol, zNumToLetter
import generateFiles as gen
#from rotate import get_rotation_matrix

import scipy.stats as stats
from sklearn import linear_model, grid_search
#from sklearn.metrics import mean_squared_error
from sklearn.kernel_ridge import KernelRidge
from sklearn.dummy import DummyRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale

#np.set_printoptions(suppress = False)


def makePredictionTable(oper, atomZ1, atomZ2, rotate = True):   
 
    pairType = "Bonded"
    orbType = "all"
    molType = ["all"]
    a1 = zNumToLetter(atomZ1)
    A1 = a1 + '1'
    a2 = zNumToLetter(atomZ2)
    A2 = a2 + '2'
=======
from generateFiles import mkdir_p
from MoleculeSet import MoleculeSet
from miscFunctions import convertToKcalPerMol
from miscFunctions import zNumToLetter
import generateFiles as gen

import scipy.stats as stats
from sklearn import linear_model, grid_search
from sklearn.metrics import mean_squared_error
from sklearn.kernel_ridge import KernelRidge


def get_axis_rotation_matrix(axis, theta):
    # http://stackoverflow.com/questions/6721544/circular-rotation-around-an-arbitrary-axis
    ct = math.cos(theta)
    nct = 1 - ct
    st = math.sin(theta)
    r = np.linalg.norm(axis)
    if r == 0.0:
        return np.matrix(np.eye(3))
    ux = axis[0] / r
    uy = axis[1] / r
    uz = axis[2] / r
    rot = np.matrix([
        [ct + ux ** 2 * nct, ux * uy * nct - uz * st, ux * uz * nct + uy * st],
        [uy * ux * nct + uz * st, ct + uy ** 2 * nct, uy * uz * nct - ux * st],
        [uz * ux * nct - uy * st, uz * uy * nct + ux * st, ct + uz ** 2 * nct],
    ])
    return rot
 
 
def get_new_operator(operator, bond_vector):
    e_x = np.array([1., 0., 0.])
    e_bond = bond_vector / np.linalg.norm(bond_vector)
    e_1 = np.cross(e_bond, e_x)
    e_2 = np.cross(e_x, e_1)
    theta_1 = math.atan2(np.dot(e_bond, e_2), np.dot(e_bond, e_x))
 
    R_1 = get_axis_rotation_matrix(e_1, theta_1)
    operator_mid = R_1 * operator * R_1.T
 
    operator_lower = operator_mid[1:, 1:]
    eigen_vals, eigen_vectors = np.linalg.eig(operator_lower)
    # Sort eigenvalues [High, ..., Low]
    idx = eigen_vals.argsort()[::-1]
    eigen_sorted = eigen_vectors[:,idx]
 
    # Lazy way to add 1 to the top of the diagonal
    R_2 = np.matrix(np.eye(3))
    R_2[1:,1:] = eigen_sorted
 
    return R_2 * operator_mid * R_2.T
 
def get_rotation_matrix(bond_vector,operator):
    e_x = np.array([1., 0., 0.])
    e_bond = bond_vector / np.linalg.norm(bond_vector)
    e_1 = np.cross(e_bond, e_x)
    e_2 = np.cross(e_x, e_1)
    theta_1 = math.atan2(np.dot(e_bond, e_2), np.dot(e_bond, e_x))
 
    R_1 = get_axis_rotation_matrix(e_1, theta_1)
    operator_mid = R_1 * operator * R_1.T
 
    operator_lower = operator_mid[1:, 1:]
    eigen_vals, eigen_vectors = np.linalg.eig(operator_lower)
    # Sort eigenvalues [High, ..., Low]
    idx = eigen_vals.argsort()[::-1]
    eigen_sorted = eigen_vectors[:,idx]
 
    # Lazy way to add 1 to the top of the diagonal
    R_2 = np.matrix(np.eye(3))
    R_2[1:,1:] = eigen_sorted
 
    return R_2 * R_1



def makePredictionTable(oper, atomZ1, atomZ2):   
 
    
    a1 = zNumToLetter(atomZ1)
    a2 = zNumToLetter(atomZ2)
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
    bondName = a1 + a2
    
    LL_set = MoleculeSet("LL", storeTwoElec = False)
    QU_set = MoleculeSet("QU", storeTwoElec = False)
    
<<<<<<< HEAD

    if rotate == True:
        LL = LL_set.make_rotated_matrix(pairType, molType, orbType, oper, atomZ1, atomZ2)
        QU = QU_set.make_rotated_matrix(pairType, molType, orbType, oper, atomZ1, atomZ2)
    else:
        LL = LL_set.getPairsAndSelves(oper, atomZ1, atomZ2)
        QU = QU_set.getPairsAndSelves(oper, atomZ1, atomZ2)
        
    table = np.zeros(len(LL[bondName]) * 36).reshape(len(LL[bondName]), 36)
    #print table.shape
    X_names = list(itertools.product(['1s', '2s', '2px', '2py', '2pz'],['1s', '2s', '2px', '2py', '2pz']))
    for i in range(0, len(X_names)):
        X_names[i] =  X_names[i][0] + X_names[i][1]
    #print C_names

    for i in xrange(0, len(LL[bondName])):

        l = 0
        options = OrderedDict()
        
        options["%s_LL_%s_1s1s" % (oper, bondName)] = LL["%s" % bondName][i][0]
        options["%s_LL_%s_2s1s" % (oper, bondName)] = LL["%s" % bondName][i][1]
        options["%s_LL_%s_2px1s" % (oper, bondName)] = LL["%s" % bondName][i][2]
        options["%s_LL_%s_2py1s" % (oper, bondName)] = LL["%s" % bondName][i][3]
        options["%s_LL_%s_2pz1s" % (oper, bondName)] = LL["%s" % bondName][i][4]
        
        if atomZ2 == 1:
            options["%s_LL_%s_1s1s" % (oper, a1)] = LL["%s" % A2][i]
        
        
            
            for j in range(0,5):   
                for k in range(0,5):
                    options["%s_LL_%s_%s" % (oper, A1, X_names[l])] = LL["%s" % A1][i][j][k]
                    l += 1
        elif atomZ1 == 1:
            options["%s_LL_%s_1s1s" % (oper, a2)] = LL["%s" % A1][i]
        
        
            
            for j in range(0,5):   
                for k in range(0,5):
                    options["%s_LL_%s_%s" % (oper, A2, X_names[l])] = LL["%s" % A2][i][j][k]
                    l += 1
        
        options["%s_%s_1s1s" % (oper, bondName)] = QU["%s" % bondName][i][0]
        options["%s_QU_%s_2s1s" % (oper, bondName)] = QU["%s" % bondName][i][1]
        options["%s_QU_%s_2px1s" % (oper, bondName)] = QU["%s" % bondName][i][2]
        options["%s_QU_%s_2py1s" % (oper, bondName)] = QU["%s" % bondName][i][3]
        options["%s_QU_%s_2pz1s" % (oper, bondName)] = QU["%s" % bondName][i][4]
    
        for j, m in enumerate(options):
#            if oper != "rho" and oper != "S":
#                table[i][j] = convertToKcalPerMol(options[m])
#            else: 
#                table[i][j] = options[m]
            table[i][j] = options[m]

=======
    LL = LL_set.getPairsAndSelves(oper, atomZ1, atomZ2)
    QU = QU_set.getPairsAndSelves(oper, atomZ1, atomZ2)
    
    table = np.zeros(len(LL["%s" % bondName]) * 36).reshape(len(LL["CH"]), 36)
    #print table.shape
    C_names = list(itertools.product(['1s', '2s', '2px', '2py', '2pz'],['1s', '2s', '2px', '2py', '2pz']))
    for i in range(0, len(C_names)):
        C_names[i] =  C_names[i][0] + C_names[i][1]
    #print C_names

    for i in xrange(0, len(LL["%s" % bondName])):
        l = 0
        options = OrderedDict()
        
        options["LL_%s_1s1s" % bondName] = LL["%s" % bondName][i][0][0]
        options["LL_%s_2s1s" % bondName] = LL["%s" % bondName][i][1][0]
        options["LL_%s_2px1s" % bondName] = LL["%s" % bondName][i][2][0]
        options["LL_%s_2py1s" % bondName] = LL["%s" % bondName][i][3][0]
        options["LL_%s_2pz1s" % bondName] = LL["%s" % bondName][i][4][0]
        
        options["LL_%s_1s1s"] = LL["%s" % a2][i][0][0]
    
    
        
        for j in range(0,5):   
            for k in range(0,5):
                options["LL_%s_%s" % C_names[l]] = LL["%s" % a1][i][j][k]
                l += 1
                
        options["QU_%s_1s1s" % bondName] = QU["%s" % bondName][i][0][0]
        options["QU_%s_2s1s" % bondName] = QU["%s" % bondName][i][1][0]
        options["QU_%s_2px1s" % bondName] = QU["%s" % bondName][i][2][0]
        options["QU_%s_2py1s" % bondName] = QU["%s" % bondName][i][3][0]
        options["QU_%s_2pz1s" % bondName] = QU["%s" % bondName][i][4][0]
    
        for j, m in enumerate(options):
    
            table[i][j] = convertToKcalPerMol(options[m])
    
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e

    return options.keys(), table


<<<<<<< HEAD
def getPredictionTable(oper, atomZ1, atomZ2, rotate = True, overwrite = False):
=======
def getPredictionTable(oper, atomZ1, atomZ2):
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e

    a1 = zNumToLetter(atomZ1)
    a2 = zNumToLetter(atomZ2)
    bond = a1 + a2
    
    dataPath = os.path.join("Data", "predictionTables", "%s" % bond, "%s" % oper )
    fileName = "%s_%s" % (bond, oper)
    
    
<<<<<<< HEAD
    if os.path.isfile(os.path.join(dataPath, fileName) + ".csv") and overwrite == False:
=======
    if os.path.isfile(os.path.join(dataPath, fileName) + ".csv"):
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        with open(os.path.join(dataPath, fileName) + ".csv", "rb") as f:
            dialect = csv.Sniffer().sniff( f.read( 10*1024 ) )
            f.seek(0)
            reader = csv.reader( f, dialect )
            op_options = reader.next()
            op_table = [row for row in reader]
    else:
<<<<<<< HEAD
        op_options, op_table = makePredictionTable(oper, atomZ1, atomZ2, rotate)
        gen.mkdir_p(dataPath)
=======
        op_options, op_table = makePredictionTable(oper, atomZ1, atomZ2)
        mkdir_p(dataPath)
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        with open(os.path.join(dataPath, fileName) + ".csv", "wb") as f:
            writer = csv.writer(f)
            writer.writerow(op_options)
            writer.writerows(op_table)
            
    op_table = np.array(op_table).astype('float')
            
    
    dataPath = os.path.join("Data", "predictionTables", "%s" % bond, "%s" % "rho" )
    fileName = "%s_%s" % (bond, "rho")
    
    
<<<<<<< HEAD
    if os.path.isfile(os.path.join(dataPath, fileName) + ".csv") and overwrite == False:
=======
    if os.path.isfile(os.path.join(dataPath, fileName) + ".csv"):
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        with open(os.path.join(dataPath, fileName) + ".csv", "rb") as f:
            dialect = csv.Sniffer().sniff( f.read( 10*1024 ) )
            f.seek(0)
            reader = csv.reader( f, dialect )
            rho_options = reader.next()
            rho_table = [row for row in reader]
    else:
<<<<<<< HEAD
        rho_options, rho_table = makePredictionTable("rho", atomZ1, atomZ2, rotate)
        gen.mkdir_p(dataPath)
=======
        rho_options, rho_table = makePredictionTable("rho", atomZ1, atomZ2)
        mkdir_p(dataPath)
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        with open(os.path.join(dataPath, fileName) + ".csv", "wb") as f:
            writer = csv.writer(f)
            writer.writerow(rho_options)
            writer.writerows(rho_table)
    
    rho_table = np.array(rho_table).astype('float')       
    
    options = rho_options[:31] + op_options
    rho_table = rho_table[:,:31]
    table = np.hstack((rho_table, op_table))
    
    return options, table

<<<<<<< HEAD
def runFits(oper, atomZ1, atomZ2, rotate = True, overwrite = False):
    a1 = zNumToLetter(atomZ1)
    a2 = zNumToLetter(atomZ2)
    bond = a1 + a2
    
    options, table = getPredictionTable(oper, atomZ1, atomZ2, rotate, overwrite)
    
    dataPath = os.path.join("Data", "predictionTables", "%s" % bond, "rho_%s" % oper )
    fileName = "%s_%s" % (bond, oper)
    
    gen.mkdir_p(dataPath)
    with open(os.path.join(dataPath, fileName) + ".csv", "wb") as f:
            writer = csv.writer(f)
            writer.writerow(options)
            writer.writerows(table)
            
    distribution = OrderedDict()
    
    for group1 in ["Original", "Dummy", "Linear_Regression", "Kernal_Ridge"]:
        distribution[group1] = OrderedDict()
        for group2 in options[62:67]:
=======
def runFits(oper, atomZ1, atomZ2):
    options, table = getPredictionTable(oper, atomZ1, atomZ2)
    
    distribution = OrderedDict()
    
    for group1 in ["Original", "Linear_Regression", "Ridge_CV", "Kernal_Ridge"]:
        distribution[group1] = OrderedDict()
        for group2 in options[62:67]:
    #    for group2 in options[31:36]:
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
            distribution[group1][group2] = OrderedDict()
            for group3 in ["Means", "Abs_errors"]:
                distribution[group1][group2][group3] = []
        
    #print distribution
    
<<<<<<< HEAD
    for j in xrange(0,2):
        print "iteration %d..." %j
        iterations = j + 1
        np.random.shuffle(table)
        
        #preprocess code:
#        X1 = table[:,:31]
#        X1 = StandardScaler().fit_transform(X1)
#        X2 = table[:,31:62]
#        X2 = StandardScaler().fit_transform(X2)
#        X = np.hstack((X1,X2))  
#        
#        X_train = X[:650,:]
#        X_test =  X[650:,:]

        X_train = table[:650,:62]
        X_test =  table[650:,:62]
        
        y_train = table[:650,62:67]
        y_test =  table[650:,62:67]

        for i, n in enumerate(distribution["Original"]):
            distribution["Original"][n]["Means"].append(y_test[:,i].mean())
    
        lin = linear_model.LinearRegression()
        lin.fit(X_train,y_train)
        
        y_pred_lin = lin.predict(X_test)
=======
    for j in xrange(0,25):
        iterations = j + 1
        np.random.shuffle(table)
        
    #    X_train = table[:700,:31]
    #    X_test =  table[700:,:31]
    #    y_train = table[:700,31:36]
    #    y_test =  table[700:,31:36]    
        X_train = table[:700,:62]
        X_test =  table[700:,:62]
        y_train = table[:700,62:67]
        y_test =  table[700:,62:67]
        
    #    print y_test.mean()
        for i, n in enumerate(distribution["Original"]):
            distribution["Original"][n]["Means"].append(y_test[:,i].mean())
    
        model = linear_model.LinearRegression()
        model.fit(X_train,y_train)
        
        y_pred_lin = model.predict(X_test)
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        
        
        for i, n in enumerate(distribution["Linear_Regression"]):
            distribution["Linear_Regression"][n]["Means"].append(y_pred_lin[:,i].mean())
            distribution["Linear_Regression"][n]["Abs_errors"].append(abs(y_pred_lin[:,i] - y_test[:,i]).mean())
        
<<<<<<< HEAD
        dum = DummyRegressor()
        
        dum.fit(X_train, y_train)
        y_pred_dum = dum.predict(X_test)
        
        for i, n in enumerate(distribution["Dummy"]):
            distribution["Dummy"][n]["Means"].append(y_pred_lin[:,i].mean())
            distribution["Dummy"][n]["Abs_errors"].append(abs(y_pred_dum[:,i] - y_test[:,i]).mean())
            
        if rotate == True:
            parameters = {'alpha':  [1e-3, 1e-2, 1e-1] , 'gamma': [1e-9, 1e-8, 1e-7]}
        else:            
            parameters = {'alpha': [1e-8, 1e-7, 1e-6], 'gamma':[1e-7, 1e-5, 1e-3, 1e-1, 1e0] }             #without prepreprocess
#            parameters = {'alpha': [1e-9, 1e-7, 1e-5, 1e-3,1e-1], 'gamma':[1e-9, 1e-7, 1e-5, 1e-3,1e-1] } #with preprocess

        rdg = KernelRidge(kernel="rbf")
=======
        clf = linear_model.RidgeCV([1e-3, 1e-1, 1e1])
        clf.fit (X_train, y_train) 
            
        y_pred_CV = clf.predict(X_test)
        for i, n in enumerate(distribution["Ridge_CV"]):
            distribution["Ridge_CV"][n]["Means"].append(y_pred_CV[:,i].mean())
            distribution["Ridge_CV"][n]["Abs_errors"].append(abs(y_pred_CV[:,i] - y_test[:,i]).mean())
        
        parameters = {'alpha': [1e-7, 1e-6, 1e-5], 'gamma':[1e-9, 1e-8, 1e-7] }
        rdg = KernelRidge(alpha=0.5, kernel="rbf")
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
        grid = grid_search.GridSearchCV(rdg, parameters)
        grid.fit(X_train, y_train) 
        y_pred_rdg = grid.predict(X_test)
        
        for i, n in enumerate(distribution["Kernal_Ridge"]):
            distribution["Kernal_Ridge"][n]["Means"].append(y_pred_rdg[:,i].mean())
            distribution["Kernal_Ridge"][n]["Abs_errors"].append(abs(y_pred_rdg[:,i] - y_test[:,i]).mean())
            
    iterations = j + 1
    print "Iterations:", iterations
    for group1 in distribution.keys():
        if group1 == "Original":
<<<<<<< HEAD
            continue
        print "\n"
        for group2 in distribution[group1].keys():

=======
            orig = distribution[group1]
            continue
        for group2 in distribution[group1].keys():
            origMeans = orig[group2]["Means"]
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
            
            means = distribution[group1][group2]["Means"]
            
            absErrorMean = np.mean(distribution[group1][group2]["Abs_errors"])
            absErrorStd = np.std(distribution[group1][group2]["Abs_errors"])
<<<<<<< HEAD
            print group1, group2, "Absolute Error: %s   +/- %s" % (convertToKcalPerMol(absErrorMean), convertToKcalPerMol(absErrorStd))
            
            mmean = np.mean(means)
            mstd = np.std(means)
            print group1, group2, "Mean: %s   +/- %s\n" % (convertToKcalPerMol(mmean), convertToKcalPerMol(mstd))
=======
            print group1, group2, "Absolute Error: %s   +/- %s" % (absErrorMean, absErrorStd)
            
            mmean = np.mean(means)
            mstd = np.std(means)
            print group1, group2, "Mean: %s   +/- %s\n" % (mmean, mstd)
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
            means.sort()
            pdf = stats.norm.pdf(means, mmean, mstd)
            plt.plot(means, pdf)
            gen.savePlotToFolder(os.path.join("Data", "predictionTables", "Distributions", group1), group2)
<<<<<<< HEAD
            
            

    
if __name__ == "__main__":  

    runFits("KE", 6, 1, rotate = False, overwrite =  False)
#    theoryType = "LL"
#    atomZ1 = 6
#    atomZ2 = 1
#    a1 = zNumToLetter(atomZ1)
#    A1 = a1 + '1'
#    a2 = zNumToLetter(atomZ2)
#    A2 = a2 + '2'
#    pair = a1 + a2
#    oper = "KE"
#    pairType = "Bonded"
#    molType = ["all"]
#    orbType = "all"
#    list1 = make_rotated_matrix(pairType, molType, orbType, oper, atomZ1, atomZ2)
#    gen.save_mats_to_file(list1[A1])
#    print makePredictionTable("KE", 6, 1)
=======
    
    
runFits("KE", 6, 1)
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
