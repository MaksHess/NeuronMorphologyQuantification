#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright (c) 2017-2018 Max Hess, Alvaro Gomariz, ETH Zurich
'''


import utility
import numpy as np

def kinkPositions(fully_annotated_mainbranch, scale=(0.22, 0.22, 0.3)):
    n_kinks = np.sum(fully_annotated_mainbranch[:,5]>0.5)
    n_outgrowths = np.sum(fully_annotated_mainbranch[:,1]==2)
    kink_positions = np.zeros(n_kinks)
    outgrowth_positions = np.zeros(n_outgrowths)
    #fully_annoated_mainbranch has nodes of raius 0.5, nodes of radius >0.5 are kinks.
    endpoints = utility.findEndpoints(fully_annotated_mainbranch, return_node=True)
    current_node = endpoints[0]
    next_node = utility.nextNode(current_node, fully_annotated_mainbranch)
    distance = 0
    i = 0
    j = 0
    
    while isinstance(next_node, np.ndarray):
        distance += utility.dist3D(current_node, next_node, scale=scale)
        if next_node[5] > 0.5: #if next_node is a kink
            kink_positions[i] = distance
            i += 1
        if next_node[1] == 2: #if next_node is a branching point
            outgrowth_positions[j] = distance
            j += 1
        current_node = next_node
        next_node = utility.nextNode(current_node, fully_annotated_mainbranch)
    return kink_positions, outgrowth_positions
    