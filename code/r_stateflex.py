#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Usage: 
    r_stateflex.py <k> <roiname> <width> <sublist>
    
Arguments:
    <k> k for k clustering
    <roiname> ROI set name
    <width> Width
    <sublist> File name of subject list
    
"""

#Imports necessary packages.
import os
from docopt import docopt
import numpy as np
import pandas as pd
import multiprocess as mp
from sklearn_extra.cluster import KMedoids

#Gets the window vectors for the subject.
def makedFCvec(subject):
    
    #Gets the vectorized dFC matrices.
    subvec = pd.read_csv('../outputs/r_vectorized/'+roiname+'/'+width+
                             '/'+subject+'/vector.csv',header=None)
    
    #Converts to numpy array for future functions.
    subvec = subvec.to_numpy()
    
    #Return window vectors for the subject.
    return subvec

#Averages k-clustered windows to get states for the cluster.
def statefind(i):
    
    #Subsets windows into the cluster.
    subset = totaldFC[(totaldFC[:,0] == i)] 
    
    #Averages the window values across a column within a cluster.
    subset = np.mean(subset[:,2:],axis=0) 
    
    #Labels window vector with cluster name.
    subset = np.insert(subset,0,i)
    
    #Return labelled vector.
    return subset

#Finds the total number of window changes between clusters for the subject.
def stateflex(clusteredsub):
    
    #Finds subject ID from first row second column. Gets rid of decimal for
    #legibility.
    subid = int(clusteredsub[0,1])
    
    #Finds total number of changes of windows between clusters.
    flexscore = len(np.where(np.diff(clusteredsub[:,0]))[0])
    
    #Returns both together.
    return [subid,flexscore]

#Required for parallel processing function.
if __name__ == '__main__':
    __spec__ = None
    
    #Catches command-line arguments.
    args = docopt(__doc__)
    
    #Imports integer k selected for calculations.
    k = int(args['<k>'])
    
    #Imports ROI name.
    roiname = args['<roiname>']
    
    #Imports window width.
    width = args['<width>']
    
    #Reads in chosen list of subject IDs.
    sublist = args['<sublist>']
    with open(sublist) as f:
        subjects = [subject.rstrip() for subject in f]
    
    #Combine subject window vectors. Remove unused variables to save memory.
    #Convert list to a numpy array to save memory.
    dFCvec = []
    for onesub in subjects:
        onesubvec = makedFCvec(onesub) #Produce subject windows.
        dFCvec.append(onesubvec)
    del onesubvec
    dFCvec = np.array(dFCvec)
    
    #Set output folder path. If output folder doesn't exist, creates it.
    basepath = '../outputs/r_stateflex/'+roiname+'/'+width+'/'+str(k)+'/'
    if not os.path.exists(basepath):   
        os.makedirs(basepath)
            
    #Get rid of subject level. Remove unused variables to save memory.
    #Convert list to a numpy array to save memory.
    totaldFC = []
    for subj in dFCvec:
        for win in subj:
            totaldFC.append(win)
    del dFCvec
    del subj
    del win
    totaldFC = np.array(totaldFC)
    
    #Conduct KMedoids clustering on the dFC vector data in the array using 
    #the Manhattan metric - slices out subject labels prior.
    kresult = KMedoids(n_clusters=k,metric='cityblock').fit(totaldFC[:,1:])
    
    #Insert the labels found from clustering into the before the subject labels
    #in the dFC vector array. Remove unused variables to save memory.
    totaldFC = np.insert(totaldFC,0,kresult.labels_,axis=1)
    del kresult
    
    #Add subject level back into the clustered dFC vector array - finds where
    #the subject label changes and splits where those occur.
    splitdFC = np.split(totaldFC,np.where(np.diff(totaldFC[:,1]))[0]+1)
    
    #Set up clusters to be processed.
    clusters = list(range(k))
    
    #Set up parallel loop.
    pool = mp.Pool()
    
    #Find states for each cluster in parallel.
    states = pool.map(statefind,clusters)
    
    #Save states in a matrix by stacking list elements. Removes unused 
    #variables to save memory.
    pd.DataFrame(states).to_csv(basepath+'states.csv',header=None,index=None)
    del totaldFC
    del states
    
    #Find state flexibility for each subject in parallel.
    flexresult = pool.map(stateflex,splitdFC)
    
    #Close parallel threads.
    pool.close()
    pool.join()
    
    #Save flexibility measures in a matrix by stacking list elements.
    pd.DataFrame(flexresult).to_csv(basepath+'stateflex.csv',header=None,
                                    index=None)
    