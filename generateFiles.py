# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 00:58:56 2015

@author: brandon
"""
import collections
import numpy as np
<<<<<<< HEAD
=======
from scipy.stats.stats import pearsonr 
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
import matplotlib.pyplot as plt
import json
import cPickle as pickle
import os
<<<<<<< HEAD
import csv
from sklearn import linear_model
from sklearn.metrics import mean_squared_error


=======
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e

def mkdir_p(path):
    
    try:
        os.makedirs(path)
    except OSError as exc:
        import errno
        if exc.errno == errno.EEXIST and os.path.isdir(path):
           pass 
        else: raise
       
       

def savePlotToFolder(folder, fileName):
    ''' saves plot as image file with specified folder and filename.'''  
    
    mkdir_p(folder)
    plt.savefig(os.path.join(folder, fileName) + '.png', bbox_inches='tight')
    
    plt.gca().cla()




def saveTextToFolder(folder, fileName, string):
    '''creates a markdown page with a specified folder and filename using a provided string'''
    
    mkdir_p(folder)
    with open(os.path.join(folder, fileName) + ".md", 'w') as f:
        f.write(string)




def saveDataToFolder_json(folder, fileName, data):
    '''saves data file to specified folder and filename using provided data'''
    
    mkdir_p(folder)
    with open(os.path.join(folder, fileName), 'wb') as f:
        json.dump(data, f)



def saveDataToFolder_pickle(folder, fileName, data):
    '''saves data file to specified folder and filename using provided data'''
    
    mkdir_p(folder)
    with open(os.path.join(folder, fileName), 'wb') as f:
        pickle.dump(data, f)




def generateDataList(title, stats, units):
    '''This function takes statistics in order to produce a markdown list. Also takes a title, and units.'''

    string = "\n\n\n### {title}:\n\n"
    string += "R^2: ({stats.R^2[0]:.6f}, {stats.R^2[0]:.6f})\n\n\n"
    string += "LL mean: {stats.mean1:.6f} {units}\n"
    string += "QU mean: {stats.mean2:.6f} {units}\n"
    string += "Mean Difference (LL - QU): {stats.MeanDifference:.6f} {units}\n\n\n"
    string += "LL (min, max): ({stats.minMax1[0]:.6f}, {stats.minMax1[1]:.6f}) {units}\n"
    string += "QU (min, max): ({stats.minMax2[0]:.6f}, {stats.minMax2[1]:.6f}) {units}\n"
    string += "Difference (min, max): ({stats.MinMaxDifference[0]:.6f}, {stats.MinMaxDifference[1]:.6f}) {units}\n\n\n"
    string += "LL std: {stats.std1:.6f} {units}\n"
    string += "QU std: {stats.std2:.6f} {units}\n"
    string += "Std Difference (LL - QU): {stats.StdDifference:.6f}"
    
   
    return string




def generateDataTable(title, stats, units):
    '''This function takes statistics in order to produce a markdown table. Also takes a title, units, and a 
    boolean indicating whether it should produce the top of the table of just more rows. 
    '''
        
    template = "<sub>%.6f</sub> | "
    string = "<b><sub>%s</sub></b> | " % title
    for key in stats.keys():
        string += template % stats[key]
    return string + "  \n"

        
    

def composeDataTable(dict1, dict2, title, xLab, yLab, units):
    '''composeDataTable creates a table with a head and tail. Uses generateDataTable for the body. '''
     
     
    stats = compareLists(dict1[dict1.keys()[0]], dict2[dict2.keys()[0]])
    
    string = " "
    for key in stats.keys():
        string += " | <sub>%s</sub>" % key

    string += "\n"    
    
    for i, key in enumerate(stats.keys()):
        if i == 0:
            string += ":---:"
        else:
            string += "|:---:"
        
    string += "|:---:  \n"
    


    for i, group1 in enumerate(dict1.keys()):
        
        N = len(dict1)
        colors = plt.cm.spectral(np.linspace(0, 1, N, endpoint = True))
        stats = compareLists(dict1[group1], dict2[group1])  
        plotLists(dict1[group1], dict2[group1], title, xLab, yLab, units, '.', group1, colors[i])
        string += generateDataTable(group1, stats, units) 

        
    string += "(%s)<br><br><br><br><br>\n\n\n" % units
                             
    return string
    
    

def compareLists(list1, list2):
    '''compareLists takes statistics of two lists'''
    stats = collections.OrderedDict()
<<<<<<< HEAD
    
#   print "R_squared(R^2):", R_squared
    R_squared, RMSE = calcCODandRMSE(list1, list2)
    
=======
    R_squared = pearsonr(list1, list2)[0]**2
#   print "R_squared(R^2):", R_squared


>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
    
    array1 = np.array([list1])
    array2 = np.array([list2])
    
    difference = array1 - array2   
    
    mean1 = np.mean(list1)
    mean2 = np.mean(list2)
    MeanDifference = np.mean(difference)
    
#    print "\nLL mean:", mean1 , units
#    print "QU mean:", mean2 , units
#    print "(LL - QU) mean:", MeanDifference , units
    
#    minMax1 = (min(list1), max(list1))
#    minMax2 = (min(list2), max(list2))
#    MinMaxDifference = (np.min(difference), np.max(difference))
    
#    print "\nLL (min, max):", minMax1 , units
#    print "QU (min, max):", minMax2 , units
#    print "(LL - QU) (min, max):", MinMaxDifference , units
    
    std1 = np.std(list1)
    std2 = np.std(list2)
    StdDifference = np.std(difference)
    
#    print "\nLL std:", std1 , units
#    print "QU std:", std2   , units
#    print "(LL - QU) std:", StdDifference, units
<<<<<<< HEAD
    stats['R^2'], stats['RMSE'] = R_squared, RMSE
=======
    stats['R^2'] = R_squared
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
    stats['LL_mean'], stats['QU_mean'], stats['Mean(LL-QU)'] = mean1, mean2, MeanDifference
#    stats['minMax1'], stats['minMax2'], stats['MinMaxDifference'] = minMax1, minMax2, MinMaxDifference
    stats['LL_std'], stats['QU_std'], stats['Std(LL-QU)'] = std1, std2, StdDifference

    return stats
    
    
    
def plotLists(list1, list2, title, xLabel, yLabel, units, pointType, key, color):
    ''' plotLists produces a scatter plot from two given lists and some parameters'''
    plt.scatter(list1, list2, marker = pointType, lw = 0.5, edgecolor = color, c = color, label = key)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints = 1, fontsize = "small")
    plt.title(title)
    plt.xlabel(xLabel + "(" + units + ")")
    plt.ylabel(yLabel + "(" + units + ")")
    plt.axes().set_aspect('equal', 'box-forced')
    plt.xticks(rotation='vertical')
    plt.grid(True)
    plt.tick_params(axis='both', which='major', labelsize=7)
#    plt.axes().set_aspect('equal', 'box-forced')
#    plt.axis('equal')
    plt.tight_layout()
    return

    
<<<<<<< HEAD
def calcCODandRMSE(list1a, list2):
    if len(list1a) == 0 or len(list2) == 0:
        return float('nan'), float('nan')
    list1b = []
    
    for i, n in enumerate(list1a):
        list1b.append([1,n])
        
    model = linear_model.LinearRegression()
    model.fit(list1b,list2)
    model.score(list1b,list2)
    return model.score(list1b,list2), mean_squared_error(list2,model.predict(list1b))

#def save_mats_to_file(mat, key1 = None, key2 = None):
#    
#    np.set_printoptions(precision = 6, suppress = True)  
#    with open("Data/saved.csv", "wb") as f:
#        writer = csv.writer(f)    
#        for i, item in enumerate(mat):
#            np.array_repr(item, precision = 6, suppress_small = True) 	
#            writer.writerows(item)
#            writer.writerows("\n")
            
            
            

def save_mats_to_file(mat, filename = "saved", key1 = None, key2 = None):
    
    np.set_printoptions(precision = 6, suppress = True)  
    with open("Data/%s.csv" % filename, "wb") as f:   
        for i, item in enumerate(mat):
            if len(item.shape) == 0:
                empty = np.array(float('nan'))
            if len(item.shape) == 1:
                empty = np.empty([1, item.shape[0]])
                empty[:] = float('nan')
            elif len(item.shape) == 2:
                empty = np.empty([1, item.shape[1]])
                empty[:] = float('nan')
            
            item = np.vstack([item, empty])
            np.savetxt(f, item, fmt='%.6f')
  
                
                
=======
>>>>>>> f44231fe3a1d627b7e0fdaf06326ce273003c85e
