# -*- coding: utf-8 -*-

import os
#os.environ["quamboData"] = '/home/yaron/code/quamboML/quambo'
#import numpy as np
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

import scipy.stats as stats
from sklearn import linear_model, grid_search
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn import linear_model
from sklearn.kernel_ridge import KernelRidge
from sklearn.grid_search import GridSearchCV
import sklearn.cross_validation as cv
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
import cPickle as pickle


from Geometry import Geometry
from Molecule import Molecule
from Molecule_set import Molecule_set
#from Molecule import Molecule
#from MoleculeSet import MoleculeSet
#from miscFunctions import convertToKcalPerMol, zNumToLetter, fixMolType

def to_scikit(alist):
    reshaped = map( lambda x: x.reshape(-1), alist)
    return np.vstack(reshaped)

def train_and_evaluate(clf, X_train, y_train):
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_train)
    print "RMS error on full training set:",np.sqrt(mean_squared_error(y_train,y_pred))
    print "Mean abs error:", mean_absolute_error(y_train,y_pred)
    #create a k-fold cross validation iterator of k=5 folds
    kf = cv.KFold(X_train.shape[0], 5, shuffle=True, random_state=33)     
    scores = cv.cross_val_score(clf, X_train, y_train, cv=kf, scoring = 'mean_squared_error')
    scores = np.sqrt(np.abs(scores))
    print "RMS of determination using 5-fold crossvalidation: ", np.mean(scores)," +- ",np.std(scores)
    

    

def get_molecule_set(data_path, theory_type, mol_nums, geom_nums):
    mset = Molecule_set(theory_type)
    for mol_num in mol_nums:
        mol = Molecule(theory_type+str(mol_num))
        basisSet = scipy.io.loadmat(os.path.join(
                                    data_path,'basisSets-%d.mat' % 
                                    (mol_num))) 
        for geom_num in geom_nums:
            # print('loading mol %d geom %d '%(mol_num,geom_num))
            zdata = scipy.io.loadmat(os.path.join(
                                    data_path,'zmats-%d-%d.mat' % 
                                     (mol_num, geom_num)))   
            qdata = scipy.io.loadmat(os.path.join(
                                    data_path,'%ss-%d-%d.mat' % 
                                     (theory_type, mol_num, geom_num)))
            g1 = Geometry(zdata, basisSet, qdata, store_two_elec = store_two_elec)
            mol.add(g1)
        mset.add(mol)
    return mset
        
if (True):
    data_path = '/home/yaron/code/quamboML/quambo';
    # There are 16 molecules available, but 6 is a repeat
    mol_nums = range(1,17)
    mol_nums.remove(6)
    #mol_nums = range(1,3)
    
    # Each molecule has 20 geometries
    geom_nums = range(1,21)
    #geom_nums = range(1,3)
    store_two_elec = False
    
    LL = get_molecule_set(data_path, 'LL', mol_nums, geom_nums)
    QU = get_molecule_set(data_path, 'QU', mol_nums, geom_nums)
    pickle.dump(LL,open("LL.p","wb"))
    pickle.dump(QU,open("QU.p","wb"))
else:
    LL = pickle.load(open("LL.p","rb"))
    QU = pickle.load(open("QU.p","rb"))

if (False):
    g = LL[1][1]
    
    Cs = g.get_atoms(6)
    Hs = g.get_atoms(1)
    CH = g.get_bonds(6,1,1)
    HC = g.get_bonds(1,6,1)
    CC = g.get_bonds(6,6,1)

allC = LL.atoms(6)
allCH = LL.bonds(6,1,1)

allCKE = LL.submatrices(allC, 'KE', '2s','2s')

Xlist = []
Xlist.append(to_scikit( LL.submatrices(allCH,'KE','all','1s') ))
Xlist.append(to_scikit( LL.submatrices(allCH,'KE','all','1s', bondAtom=0) ))
Xlist.append(to_scikit( LL.submatrices(allCH,'KE','all','1s', bondAtom=1) ))
Xlist.append(to_scikit( LL.submatrices(allCH,'rho','all','1s') ))
Xlist.append(to_scikit( LL.submatrices(allCH,'rho','all','1s', bondAtom=0) ))
Xlist.append(to_scikit( LL.submatrices(allCH,'rho','all','1s', bondAtom=1) ))


ylist = []
ylist.append(to_scikit( QU.submatrices(allCH,'KE','2s2p','1s') ))

X = np.hstack(Xlist) * 627.7
y = np.hstack(ylist) * 627.7
#scalerX = StandardScaler().fit(X_train)
#scalery = StandardScaler().fit(y_train)
#X_train = scalerX.transform(X_train)


#clf_sgd = linear_model.SGDRegressor(loss='squared_loss', 
 #   penalty='l2',  random_state=42)
    
#print 'LINEAR REGRESSION'
#train_and_evaluate(clf_sgd, X,y[:,1])

Xtrain, Xtest, ytrain, ytest = cv.train_test_split(X, y, test_size = 0.25, random_state = 20)

#clf = GridSearchCV(KernelRidge(kernel="rbf"), 
#                   {'alpha': [1e-7, 1e-6, 1e-5],
#                    'gamma': [1e-9, 1e-8, 1e-7]})

clf = GridSearchCV(SVR(kernel="rbf", epsilon= 0.01), 
                   {'C': [1e4, 1e6, 1e8],
                   'gamma': [1e-8, 1e-7, 1e-6]})

iorb = 0
y1 = ytrain[:,iorb]
clf.fit(Xtrain,y1)
y1_pred = clf.predict(Xtrain)
print 'mean absolute error on train = ', mean_absolute_error(y1,y1_pred)

y2 = ytest[:,iorb]
y2_pred = clf.predict(Xtest)
print 'mean absolute error on test = ', mean_absolute_error(y2,y2_pred)




