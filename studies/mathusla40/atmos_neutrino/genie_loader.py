import collections
import os

import ROOT
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import numpy as np
from tqdm import tqdm


def extract_array(StdHep, n=4):
    if n>1:
        return np.reshape([StdHep[i] for i in range(len(StdHep))], (len(StdHep)//n, n))
    else:
        return np.array([StdHep[i] for i in range(len(StdHep))])
        
        
def get_final(Tree):
    N = Tree.StdHepN # Number of entries in this event
    fd,ld = extract_array(Tree.StdHepFd, n=1), extract_array(Tree.StdHepLd, n=1)
    inds = np.flatnonzero((fd==-1)&(ld==-1)) # find index of outgoing particle
    
    finals_p4=[]
    finals_pdg=[]
    for i in inds:
        finals_p4.append([Tree.StdHepP4[i*4+0],Tree.StdHepP4[i*4+1],Tree.StdHepP4[i*4+2],Tree.StdHepP4[i*4+3]] )
        finals_pdg.append(Tree.StdHepPdg[i])
    return finals_p4, finals_pdg