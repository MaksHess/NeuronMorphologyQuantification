#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright (c) 2017-2018 Max Hess, Alvaro Gomariz, ETH Zurich
Load a APP2 traced .swc file and clean it up, save a cleaned .swc file.
'''

import utility
import numpy as np
import matplotlib.pyplot as plt
from classify import traceBranch
np.set_printoptions(precision=2, suppress=True)


def interpolateNodes(start, end, idx, radius=None):
    """
    Interpolate nodes between `start` and `end`.

    Parameters
    ----------
    start : np.ndarray
        Start node.
    end : np.ndarray
        End node.
    idx : int
        Lower bound to start index count for interpolated nodes.
    radius : int or None
        If int, constant radius for interpolated nodes, else radius is interpolated linearly.

    Returns
    -------
    nodes : np.ndarray
        Interpolated nodes.
    """

    distance_in_pixels = utility.dist3D(start, end)
    nodes = np.zeros((int(distance_in_pixels), 7))
    nodes[:,0] = np.arange(idx+1, idx+1+len(nodes))
    nodes[:,6] = nodes[:,0] - 1
    if end[2]-start[2] != 0:
        nodes[:,2] = np.linspace(start[2], end[2], len(nodes))
    else:
        nodes[:,2] = start[2]
    if end[3]-start[3] != 0:
        nodes[:,3] = np.linspace(start[3], end[3], len(nodes))
    else:
        nodes[:,3] = start[3]
    if end[4]-start[4] != 0:
        nodes[:,4] = np.linspace(start[4], end[4], len(nodes))
    else:
        nodes[:,4] = start[4]
    nodes[:,1] = start[1]
    if radius:
        nodes[:,5] = radius
    else:
        if end[5]-start[5] != 0:
            nodes[:,5] = np.linspace(start[5], end[5], len(nodes))
        else:
            nodes[:,5] = start[5]
    return nodes

def cleanup(tree,
            neurontype='PLM',
            scale=(0.223, 0.223, 0.3),
            visualize=False):
    """
    Cleanup tracing errors where outgrowth events move parallel along the mainbranch.

    Parameters
    ----------
    tree : np.ndarray
        Tree on wich to perfom cleanup.
    neurontype : str {'ALM', 'PLM'}
        Process ALM or PLM neurons.
    scale :  tuple of floats
        x, y and z scales of the images underlying the analysis.
    visualize : bool
        Wheter to visualize the results.

    Returns
    -------
    full_clean_tree:
        Clean tree.
    """
    
    endpoints = utility.findEndpoints(tree)
    
    # For ALM neurons detect the soma_nodes i.e. all nodes connected to the root that have radius above a threshold
    if neurontype=='ALM':
        soma_nodes = utility.findSomaNodes(tree, scale=scale)
    else:
        soma_nodes = []
    
    # Trace from every endpoint to a root and save the corresponding branches, select the longest as mainbranch
    branches = []
    lengths = np.zeros(len(endpoints))
    for i in range(len(endpoints)):
        branch, length = traceBranch(endpoints[i], tree, soma_nodes=soma_nodes, scale=scale)
        branches.append(branch)
        lengths[i] = length
    mainbranch = branches[lengths.argmax()]
    mainbranch_length = lengths.max()-mainbranch[-1][5]*scale[0] #the last node is part of the soma and its radius gets subtracted from the final length
    
    
    # Trace from every endpoint to a node on the mainbranch to find sidebranches
    side_branches = []
    side_lengths = np.zeros(len(endpoints))
    for i in range(len(endpoints)):
        branch, length = traceBranch(endpoints[i], tree, main_nodes=mainbranch, soma_nodes=soma_nodes, scale=scale)
        side_branches.append(np.flip(branch, axis=0))
        side_lengths[i] = length-branch[-1][5]*scale[0]
    
    
    # check if sidebranches are close and parallel to mainbranch
    if visualize:
        fig, axes = plt.subplots(3, 1, sharex='col')
    all_distances = []
    all_slopes = []
    clean_side_branches = []
    windows = []
    for side_branch in side_branches:
        root = utility.findRoots(side_branch, return_node=True)[0]
        if root.tolist() in soma_nodes:
            window = [root] #set the searching window to the root node in case of alm soma_outgrowth side_branch.
        else:
            window = utility.findWindow(root, mainbranch, window_size=40, scale=scale)
        windows.append(window)
        min_distance_from_mainbranch = []
        min_distance_from_mainbranch2 = []
        for node in side_branch:
            distances = []
            distances2 = []
            for main_node in window:
                distances.append(utility.dist3D(node, main_node, scale=scale))
                distances2.append(utility.dist3DWithRadius(node, main_node, scale=scale))
            min_distance_from_mainbranch.append(min(distances))
            min_distance_from_mainbranch2.append(min(distances2))
            if min(distances)>5:
                break
        #min_distance_from_mainbranch = min_distance_from_mainbranch[5:]
        #min_distance_from_mainbranch2 = min_distance_from_mainbranch2[5:]
        if visualize:
            axes[0].plot(min_distance_from_mainbranch)
            axes[1].plot(min_distance_from_mainbranch2)
            axes[2].plot(np.cumsum(min_distance_from_mainbranch2)/(np.arange(len(min_distance_from_mainbranch2))+1))
        all_distances.append(min_distance_from_mainbranch2)

    if visualize:    
        plt.show()

    start_node_index = np.zeros(len(all_distances))
    d_th = 0.25
    for i, distance in enumerate(all_distances):
        for j, dist in enumerate(distance):
            if dist < d_th:
                start_node_index[i] = j

    #start_node_index[start_node_index<5]=0
    idx = np.max(tree[:, 0])+1
    for i in range(len(start_node_index)):
        if start_node_index[i] == 0:
            clean_side_branches.append(side_branches[i])
        else:
            new_side_branch = side_branches[i]
            radius = np.mean(new_side_branch[:int(start_node_index[i])+1, 5])
            new_side_branch = new_side_branch[int(start_node_index[i]):]
            
            
            window = windows[i]
            distances = np.zeros(len(window))
            for i in range(len(window)):
                distances[i] = utility.dist3D(new_side_branch[0], window[i])
                
                
            connection_node = window[distances.argmin()]
            
            nodes = interpolateNodes(connection_node, new_side_branch[0], idx, radius)
            idx = np.max(nodes[:, 0])+1
            new_side_branch[0][6] = nodes[-1][0]
            nodes[0][6] = connection_node[0]
            real_side_branch = np.concatenate((nodes, new_side_branch))
            #tree = np.concatenate((tree, nodes))
            clean_side_branches.append(real_side_branch)
    
    
    #connect everything again and save clean .swc file
    full_clean_tree = []
    for node in mainbranch:
        full_clean_tree.append(node)
    for node in soma_nodes:
        full_clean_tree.append(node)
    for clean_side_branch in clean_side_branches:
        for node in clean_side_branch:
            full_clean_tree.append(node)
    
    full_clean_tree = np.array(full_clean_tree)
    return utility.removeDoubleNodes(full_clean_tree)
    
if __name__ == '__main__':
    cleanup()