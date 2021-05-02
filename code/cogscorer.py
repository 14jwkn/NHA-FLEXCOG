#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Finds cognitive scores.

"""

#Imports necessary packages.
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
def main():
    
    #Get all columns of interest (subject, demograpics, cognitive test) from 
    #the HCP1200 data dictionary, and the list of all 100 subjects of interest.
    inpath = '../inputs/data/HCP1200_Data_Dictionary.csv'
    cogscores = pd.read_csv(inpath,header=0,usecols=[
                            'Subject',
                            'Gender',
                            'Age',
                            'CogFluidComp_Unadj',
                            'CogCrystalComp_Unadj',
                            'CardSort_Unadj',
                            'Flanker_Unadj',
                            'ListSort_Unadj',
                            'PicSeq_Unadj',
                            'PicVocab_Unadj',
                            'ProcSpeed_Unadj',
                            'ReadEng_Unadj',
                            'PMAT24_A_CR',
                            'IWRD_TOT',
                            'VSPLOT_TC'
                            ])
    with open('r_sub100.txt') as f:
            subjects = [subject.rstrip() for subject in f]
    
    #Convert subjects of interest gotten as strings into integers in order to
    #find matching subjects in the dataframe which contains integers.
    numsubint = list(map(int,subjects))
    cogint = cogscores[cogscores.Subject.isin(numsubint)]
    
    #Set the specific cogntiive test names needed for the g-factor.
    cogtests = ['CardSort_Unadj',
                'Flanker_Unadj',
                'ListSort_Unadj',
                'PicSeq_Unadj',
                'PicVocab_Unadj',
                'ProcSpeed_Unadj',
                'ReadEng_Unadj',
                'PMAT24_A_CR',
                'IWRD_TOT',
                'VSPLOT_TC']
    
    #Slice out corresponding columns and convert it to numpy for processing.
    pgtable = cogint.loc[:,cogtests]
    ngtable = pgtable.to_numpy()
    
    #Z-score standardization of scores by demeaning and dividing by sample
    #standard deviation.
    mgtable = ngtable - ngtable.mean(axis=0,keepdims=True)
    zgtable = mgtable / ngtable.std(axis=0,ddof=1,keepdims=True)
    
    #Conduct PCA and set the first principal component score as g-factor score.
    pca = PCA()
    gscore = pca.fit_transform(zgtable)[:,0]
    
    #Add the score to the list and isolate the scores of interest. Based on 
    #negative loadings of all cognitive tests on first principal component,
    #convert to positive for better interpretability (i.e. as g-factor score
    #increases, other cognitive scores should increase.
    outcogtable = cogint.assign(gPCA=-gscore)
    outcogtable = outcogtable.loc[:,[
                        'Subject',
                        'Gender',
                        'Age',
                        'gPCA',
                        'CogFluidComp_Unadj',
                        'PMAT24_A_CR',
                        'CogCrystalComp_Unadj',
                        'PicVocab_Unadj'
                        ]]
    
    #Verify PCA validity, first sets the names of the components.
    components = ['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']
    
    #Compute factor loadings and creates data frame from it.
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    loading_mat = pd.DataFrame(loadings,columns=components,index=cogtests)
    
    #Add row for explained variance of the component.
    loading_mat.loc['Variance Explained'] = pca.explained_variance_ratio_
    
    #Save both cognitive scores and factor loading/explained variance for 
    #PCA components to outpath.
    outpath = '../outputs/cognition/'
    outcogtable.to_csv(outpath +'cogscores.csv',index=False)
    loading_mat.to_csv(outpath+'cogPCAvalidity.csv',header=True,index=True) 
main()
        