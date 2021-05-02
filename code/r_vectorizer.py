#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Usage: 
    r_vectorizer.py <subject> <width> <rois> 
    
Arguments:
    <subject> Subject ID
    <width> Width
    <rois> ROIs file (name of set, ROI index, ROI index2,...)
    
"""

#Get all required packages.
import os
from docopt import docopt
import pandas as pd
import multiprocess as mp

#Make the vectorized windows for each subject.
def makedFCvec(singcols):
    
    #Read in submatrix.
    submatrix = pd.read_csv('../outputs/r_slid_dFC/'+width+'/'+subject+
                            '/r_dFC100.csv',usecols=singcols,header=None)
    
    #Add rois as column names.
    submatrix.columns = rois
        
    #Label vector with subject then sequentially add connections using 
    #the row and column indices of the connections.
    subvec = [subject]
    for indices in wantindices:
        subvec.append(submatrix.loc[indices[0],indices[1]])
    
    #Return the vectorized windows for the subject.
    return subvec

#Required for parallel processing function.
if __name__ == '__main__':
    __spec__ = None
    
    #Catches command-line arguments.
    args = docopt(__doc__)
    
    #Reads in subject ID.
    subject = args['<subject>']
    
    #Imports window width, then finds number of windows based on it.
    width = args['<width>']
    nwin = (1200 - int(width) + 1)*2
    
    #Imports the list of ROIs and finds their number.
    roilist = args['<rois>']
    roilist = 'roi_whole.txt'
    with open(roilist) as f:
        labrois = [subject.rstrip() for subject in f] #Reads into separate cells.
    roiname = labrois.pop(0) #Remove the roi label from the first line.
    rois = list(map(int,labrois)) #Convert to integers to allow calculation.
    nrois = len(rois) #Finds number of rois.
    
    #Finds the indices of each unique connection between 360 regions.
    wantindices = []
    empty = pd.DataFrame(index=list(range(360)),columns=rois)
    for row,allconn in empty.iterrows():
        for col,conn in allconn.iteritems():
            
            #For all possible combinations of indices, if it is not a 
            #repeat connection and it is not an auto-connection, then
            #append it to a list.
            if (row != col) and ([col,row] not in wantindices):
                wantindices.append([row,col])
    
    #First, add the the roi indices of the first window to the list.
    firstlist = rois[:]
    
    #Then, for all windows beyond the first window which is subtracted out, 
    #and for each roi, process them. Remove unused variables to save memory.
    wantcols = [firstlist]
    for win in range(nwin-1):
        subcols = []
        for roival in rois:
                
            #Each index for a roi is 360 indices away from the next index.
            #If you multiply the current window number (e.g. 1 for first 
            #window, by adding 1 to first value of 0 in the range object) 
            #by 360, you will get the next index by addition. Each row index 
            #is saved in the list.
            subcols.append(roival+(360*(win+1)))
        wantcols.append(subcols)
    del subcols
    
    #Set up parallel loop for vectorizing connections for each window. Remove
    #unused variables to save memory.
    pool = mp.Pool()
    dFCvec = pool.map(makedFCvec,wantcols) #Vectorize windows for each subject.
    del wantcols
    
    #Closes parallel threads.
    pool.close()
    pool.join()
    
    #If output folder doesn't exist, creates it.
    basepath = '../outputs/r_vectorized/'+roiname+'/'+width+'/'+subject+'/'
    if not os.path.exists(basepath):
        os.makedirs(basepath)
    
    #Saves the vectorized windows as a csv.
    pd.DataFrame(dFCvec).to_csv(basepath+'vector.csv',header=None,index=None)
