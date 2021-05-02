#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Finds subjects of interest.

"""

#Get all required packages.
import pandas as pd
import numpy as np
import os

def main():

    #Get all participants with all MSMAll-processed 3T REST files.
    allrest = []
    for subject in os.listdir("../outputs/meants/"):
        L1path = ("../outputs/meants/" + subject + 
                  "/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv")
        R1path = ("../outputs/meants/" + subject + 
                  "/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv")
        L2path = ("../outputs/meants/" + subject + 
                  "/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv")
        R2path = ("../outputs/meants/" + subject + 
                  "/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv")
        if (os.path.exists(L1path) and os.path.exists(R1path) and os.path.exists(L2path)
            and os.path.exists(R2path)):
            allrest.append(subject)
    
    #Sort subjects into order for readability.
    allrest.sort()
    
    #Get all participants with all cognitive test scores of interest. 
    HCPdict = pd.read_csv("../inputs/data/HCP1200_Data_Dictionary.csv", index_col=0)
    HCPabrv = HCPdict[["CogFluidComp_Unadj",
                       'CogCrystalComp_Unadj',
                       "CardSort_Unadj",
                       "Flanker_Unadj",
                       "ListSort_Unadj",
                       "PicSeq_Unadj",
                       "PicVocab_Unadj",
                       "ProcSpeed_Unadj",
                       "ReadEng_Unadj",
                       "PMAT24_A_CR",
                       "IWRD_TOT",
                       "VSPLOT_TC"]]
    
    #Drop all participants with NAs in cogntiive score.
    cleanHCPabrv = HCPabrv.dropna(axis=0)
    
    #Get the subject id and convert it to a string for later matching.
    allcog = cleanHCPabrv.index
    allcog = list(map(str,list(allcog)))
    
    #Gets all participants who have had the same pre-processing done.
    sameprocess = HCPdict[(HCPdict.Release != 'Q1')& 
                          (HCPdict.Release != 'Q2')&
                          (HCPdict.Release != 'Q3')]
    
    #Get the subject id and convert it to a string for later matching.
    allproc = sameprocess.index
    allproc = list(map(str,list(allproc)))
    
    #Produces a sorted list of the intersection between the participants.
    complist = sorted(list(set.intersection(*map(set,[allcog,allrest,allproc]))))
    
    #Saves the list.
    np.savetxt('r_sub.txt', complist, delimiter="\n", fmt="%s")
main()
