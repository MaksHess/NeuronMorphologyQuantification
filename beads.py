import utility
import matplotlib.pyplot as plt
import numpy as np
from classify import traceBranch
import sys


def countBeads(mainbranch, annotations=None, suppression_window=5, local_window=32, scale=(0.22, 0.22, 0.3)):
    mean_thickness = 0
    std_thickness = 0.54 #empirically determined for some neurons
    mainbranch_beads = mainbranch.copy()
    mainbranch_beads[:,1] = 0
    mean_distance = np.mean(utility.calculateDistancesTree(mainbranch, scale=scale, return_sum=False))
    
    endpoints = utility.findEndpoints(mainbranch)
    localThickness = np.zeros(len(mainbranch))
    localThickness2 = localThickness.copy()

    n_indices = int(round(local_window/mean_distance/2))
    for i, node in enumerate(mainbranch):
        window = mainbranch[max(0, (i-n_indices)) : min(len(mainbranch), i+n_indices)]
        localThickness[i] = np.mean(window[:, 5])
   
    realThickness = mainbranch[:,5] - localThickness
    #realThickness2 = mainbranch[:, 5] / localThickness
    debugThickness = realThickness.copy()
    
    mask_existing_structures = (annotations[:, 5]==3) | (annotations[:, 1]!=0)
    it = int(np.floor(suppression_window/2/mean_distance))
    mask_existing_structures = utility.dilate_array(mask_existing_structures, it)
    realThickness = realThickness * np.invert(mask_existing_structures)
    debugThickness2 = realThickness.copy()
    i = 0
    while realThickness.max() > (mean_thickness + 3*std_thickness):
        #print(i)
        mainbranch_beads[realThickness.argmax(), 1] = 1
        mainbranch_beads[realThickness.argmax(), 5] *= 2
        lower = int(realThickness.argmax()-np.floor(suppression_window/2/mean_distance))
        upper = int(realThickness.argmax()+np.floor(suppression_window/2/mean_distance)) + 1
        if lower < 0:
            realThickness[0:upper] = 0
        elif upper > len(realThickness):
            realThickness[lower:] = 0
        else:
            realThickness[lower:upper] = 0
        i += 1
    
    bead_count = sum(mainbranch_beads[:,1])
    return mainbranch_beads, bead_count

if __name__ == '__main__':
    mainbranch = utility.readSWC('mainbranch.swc')
    ann = utility.readSWC('ann.swc')
    m, c = countBeads(mainbranch, ann)
    print(c)