#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright (c) 2017-2018 Max Hess, Alvaro Gomariz, ETH Zurich
Classifie all the nodes in a neuron into mainbranch, neuriteoutgrowth, somanodes, somaoutgrowth, pmv and blobs
'''

import utility
import numpy as np
import itertools
#import matplotlib.pyplot as plt
#from scipy.stats import linregress #as linregress
#import os

def traceBranch(endpoint, tree, main_nodes=[], soma_nodes=[], scale=(1,1,1)):
    '''Trace from an endpoint to a root node or any other node specified in main_branch  or soma_nodes.'''
    if isinstance(endpoint, (float, int)):
        endpoint_index = endpoint
    else:
        endpoint_index = endpoint[0]
    if isinstance(tree, list):
        tree = np.array(tree)
    if isinstance(main_nodes, np.ndarray):
        main_nodes = main_nodes.tolist()
    if isinstance(soma_nodes, np.ndarray):
        soma_nodes = soma_nodes.tolist()
 
    
    #any tracing eventually stops in a root_node
    root_nodes_array = utility.findRoots(tree, return_node=True)
    root_nodes = root_nodes_array.tolist()

    
    branch = []
    length = 0
    current_node = utility.thisNode(endpoint_index, tree, as_list=False)
    branch.append(current_node)
    count = 0
    while current_node.tolist() not in root_nodes and current_node.tolist() not in main_nodes and current_node.tolist() not in soma_nodes:
        count += 1
        if count > 10000:
            break

        next_node = utility.nextNode(current_node, tree)
        dist = utility.dist3D(current_node, next_node, scale=scale)
        try:
            length += dist
            current_node = next_node
            branch.append(current_node)
        except TypeError:
            break
    
    branch_array = np.array(branch)
    return branch_array, length


def classify(tree,
             neurontype='PLM',
             length_threshold=3,
             plm_outgrowth_suppression=5,
             scale=(0.223, 0.223, 0.3),
             debug=False):
    '''
    '''

    tree = utility.removeDoubleNodes(tree)
    endpoints = utility.findEndpoints(tree)
    
    #For ALM neurons detect the soma_nodes i.e. all nodes connected to the root that are above a threshold
    if neurontype=='ALM':
        soma_nodes = utility.findSomaNodes(tree, scale=scale).tolist()
    else:
        soma_nodes = []
     
    #Trace from every endpoint to a root and save the corresponding branches, select the longest as mainbranch
    branches = []
    lengths = np.zeros(len(endpoints))
    for i in range(len(endpoints)):
        branch, length = traceBranch(endpoints[i], tree, soma_nodes=soma_nodes, scale=scale)
        branches.append(branch)
        lengths[i] = length
    if not np.any(branches[lengths.argmax()][:, 1]==4):
        mainbranch = branches[lengths.argmax()]
        mainbranch_length = lengths.max()-mainbranch[-1][5]*scale[0] #the last node is part of the soma and its radius gets subtracted from the final length
    else:
        br = [branch for branch in branches if not np.any(branch[:, 1]==4)]
        le = np.array([lengths[i] for i, branch in enumerate(branches) if not np.any(branch[:, 1]==4)])
        mainbranch = br[le.argmax()]
        mainbranch_length = le.max()-mainbranch[-1][5]*scale[0]
    
    #To counteract false neurite detection at the end of PLM neurons due to neuron tracing errors at the edges the last couple of nodes
    #in PLM neurons are silenced for outgrowth detection (i.e. out growths detected there are dubbed SilencedOutgrowths and excluded from analysis).
    if neurontype=='PLM':
        silenced_nodes = utility.findSilencedNodes(mainbranch, scale=scale, max_distance=plm_outgrowth_suppression).tolist()
    else:
        silenced_nodes = []    
    
    #Trace from every endpoint to a node on the mainbranch to find sidebranches
    side_branches = []
    side_lengths = np.zeros(len(endpoints))
    for i in range(len(endpoints)):
        branch, length = traceBranch(endpoints[i], tree, main_nodes=mainbranch, soma_nodes=soma_nodes, scale=scale)
        side_branches.append(np.flip(branch, axis=0))
        side_lengths[i] = length-branch[-1][5]*scale[0]
    
    
    #Concatenate side branches that end in the same final_node to side_trees
    
    final_nodes = []
    for side_branch in side_branches:
        final_nodes.append(side_branch[0].tolist())
    final_nodes.sort()
    final_nodes = list(final_nodes for final_nodes, _ in itertools.groupby(final_nodes))
    
    
    side_trees = []
    for final_node in final_nodes:
        side_tree = []
        for side_branch in side_branches:
            if side_branch[0].tolist()==final_node:
                side_tree.append(side_branch)
        side_tree_noduplicates = [node.tolist() for branch in side_tree for node in branch]
        side_tree_noduplicates.sort()
        side_tree_noduplicates = list(node for node, _ in itertools.groupby(side_tree_noduplicates))
        side_trees.append(side_tree_noduplicates)
    
    side_trees_array = [np.array(side_tree) for side_tree in side_trees]
        
            
    
    
    
    #Find PMV branch in PLM neurons
    if neurontype=='PLM':
        pmv_nodes = []
        for i in range(len(side_trees_array)):
            maximum = side_trees_array[i][:,5].max() 
            mean = side_trees_array[i][:,5].mean()
            length = utility.calculateDistancesTree(side_trees_array[i], scale=scale)
            string = '{0:.3f}   {1:.3f}   {2:.3f}'.format(maximum, mean, length)
            if (maximum > 5 and length > 12) or np.any(side_trees_array[i][:, 1]==2):
                pmv_node = utility.findRoots(side_trees_array[i], return_node=True)[0]
        try:
            pmv_nodes.extend(utility.findWindow(pmv_node, mainbranch, window_size=10, scale=scale).tolist())
        except NameError:
            pass
    else:
        pmv_nodes = []
    
    #classify the side_trees according to their final node (pmv_nodes -> pmv-branch, soma_nodes -> soma-outgrowth, main_nodes -> neurite-outgrowth)
    annotated_mainbranch = mainbranch.copy()
    annotated_mainbranch[:,1] = 0
    annotated_mainbranch[:,5] = 1
    main_nodes = mainbranch.tolist()
    side_category = np.zeros(len(side_trees_array))
    last_nodes = []
    tree_lengths = []
    tree_classes = []
    tree_mean_radii = []
    tree_max_radii = []
    tree_orders = []
    side_trees_clean = []
    
    for i in range(len(side_trees_array)):
        tree_length = utility.calculateDistancesTree(side_trees[i], scale=scale, return_sum=True)
        tree_mean_radius = side_trees_array[i][:,5].mean()
        tree_max_radius = side_trees_array[i][:,5].max()
        
        if tree_length > length_threshold:
            root = utility.findRoots(side_trees[i], return_node=True)[0].tolist()
            last_nodes.append(root)
            tree_lengths.append(tree_length)
            tree_mean_radii.append(tree_mean_radius)
            tree_max_radii.append(tree_max_radius)
            
            if root in pmv_nodes:
                side_trees_array[i][:,1] = 6
                side_category[i] = 6
                tree_classes.append('PMV')
                annotated_mainbranch[np.argwhere(annotated_mainbranch[:,0]==root[0])[0][0], 1] = 6
                    
            elif root in soma_nodes:
                side_trees_array[i][:,1] = 3
                side_category[i] = 3
                tree_classes.append('SomaOutgrowth')
            
            elif root in silenced_nodes:
                side_trees_array[i][:,1] = 10
                side_category[i] = 10
                tree_classes.append('SilencedOutgrowth')
                annotated_mainbranch[np.argwhere(annotated_mainbranch[:,0]==root[0])[0][0], 1] = 10
                
            elif root in  main_nodes and (np.any(side_trees_array[i][:, 1]==1) or np.any(side_trees_array[i][:, 1]==3)):
                side_trees_array[i][:,1] = 7
                side_category[i] = 7
                tree_classes.append('PMV_Connection')
                annotated_mainbranch[np.argwhere(annotated_mainbranch[:,0]==root[0])[0][0], 1] = 7
            
            elif root in main_nodes and ((side_trees_array[i][:,5].mean() > 2 and tree_length < 5) or side_trees_array[i][:,5].mean() > 3) or (np.any(side_trees_array[i][:, 1]==5)):
                side_trees_array[i][:,1] = 5
                side_category[i] = 5
                tree_classes.append('Blob')
                annotated_mainbranch[np.argwhere(annotated_mainbranch[:,0]==root[0])[0][0], 1] = 5
        
            elif root in main_nodes:
                side_trees_array[i][:,1] = 2
                side_category[i] = 2
                tree_classes.append('NeuriteOutgrowth')
                annotated_mainbranch[np.argwhere(annotated_mainbranch[:,0]==root[0])[0][0], 1] = 2
            
            else:
                side_trees_array[i][:,1] = 9
                side_category[i] = 9
                tree_classes.append('Unknown')                
                
            side_trees_clean.append(side_trees_array[i])
    
    
    
    
    
    #save the final classified tree
    full_classified_tree = []
    mainbranch[:,1] = 1
    main_nodes = mainbranch.tolist()
    soma_nodes = np.array(soma_nodes)
    try:
        soma_nodes[:,1] = 0
    except IndexError:
        pass
    soma_nodes = soma_nodes.tolist()
    for node in main_nodes:
        full_classified_tree.append(node)
    for tree in side_trees_clean:
        for node in tree:
            full_classified_tree.append(node.tolist())
    for node in soma_nodes:
        full_classified_tree.append(node)
    full_classified_tree_array = np.array(full_classified_tree)
    
    tree_lengths.append(mainbranch_length)
    tree_classes.append('MainBranch')
    tree_mean_radii.append(mainbranch[:,5].mean())
    tree_max_radii.append(mainbranch[:,5].max())
    if debug:
        for i in range(len(tree_lengths)):
            
            print('Class: {0:16}    Length: {1:>6.2f}    Max_r: {2:>4.2f}    Mean_r: {3:>4.2f}'.format(tree_classes[i], tree_lengths[i], tree_max_radii[i], tree_mean_radii[i]))
            
    return tree_lengths, tree_classes, tree_mean_radii, tree_max_radii, mainbranch, annotated_mainbranch, full_classified_tree_array


if __name__ == '__main__':
    tree = utility.readSWC(r'E:\debug_data\ALM\traces_manually\COL10_D21_S4_ALM08_gs0-7.swc')
    tree_lengths, tree_classes, tree_mean_radii, tree_max_radii, mainbranch, annotated_mainbranch, full_classified_tree_array = classify(tree, neurontype='ALM', debug=True)
    a = 3