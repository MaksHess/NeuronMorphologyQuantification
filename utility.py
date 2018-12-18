#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright (c) 2017-2018 Max Hess, Alvaro Gomariz, ETH Zurich
Some much used utility functions.
'''

import os
import numpy as np

def dilate_array(array, n):
    '''Perform n iteratios of a morphological dilation on a 1-D boolean array.'''
    if not np.all((array==1)|(array==0)):
        raise ValueError('Only booolean arrays ar valid')
    
    out = np.zeros(array.shape, dtype=bool)
    for idx in np.where(array==1)[0]:
        out[max(0, idx-n):min(len(out), idx+n+1)] = 1
    return out


def readSWC(filename, as_list=False):
    '''Read an .swc file into a np.array or a list of lists.'''
    with open(filename, 'r') as f:
        trace = []
        for line in f:
            if line[0] != '#':
                strings = line.split()
                nums = [float(c) for c in strings]
                trace.append(nums)
    trace_array = np.array(trace)
    if as_list:
        return trace
    else:
        return trace_array


def findStartPLM(stack):
    '''Take an input stack and return a pixel at the edge of the stack with maximum value that serves as a root for PLM neuron tracing.'''
    max_values = np.zeros(stack.shape[0])
    max_coordinates = np.zeros((stack.shape[0],2))
    for i in range(stack.shape[0]):
        image = stack[i]
        image[1:image.shape[0]-1,1:image.shape[1]-1] = 0
        max_values[i] = image.max()
        coords = np.unravel_index(image.argmax(), image.shape)
        max_coordinates[i,0]=coords[0]; max_coordinates[i,1]=coords[1]
    z = max_values.argmax()
    x, y = max_coordinates[z]
    manual_start = (int(z),int(x),int(y))
    return manual_start



def findStartPLMwithBorder(stack):
    max_values = np.zeros(stack.shape[0])
    max_coordinates = np.zeros((stack.shape[0],2))
    #for i in range(stack.shape[0]):
        #image = stack[i]
        #image[1:image.shape[0]-1,1:image.shape[1]-1] = 0
        #max_values[i] = image.max()
        #coords = np.unravel_index(image.argmax(), image.shape)
        #max_coordinates[i,0]=coords[0]; max_coordinates[i,1]=coords[1]
    for i in range(stack.shape[0]):
        im = stack[i]
        im[3:im.shape[0]-3,3:im.shape[1]-3] = 0
        im2 = im[2:-2, 2:-2]
        im3 = np.pad(im2, 2, 'constant')
        max_values[i] = im3.max()
        coords = np.unravel_index(im3.argmax(), im3.shape)
        max_coordinates[i,0]=coords[0]; max_coordinates[i,1]=coords[1]
    z = max_values.argmax()
    x, y = max_coordinates[z]
    manual_start = (int(z)+1,int(x)+1,int(y)+1)
    return manual_start

def saveSWC(filename, branch, set_roots_to_minus_one=True):
    '''Create an empty file and write branch in .scw format, set all the root nodes parents nodes to -1.'''
    if set_roots_to_minus_one:
        branch = np.array(branch)
        indices = branch[:,0]
        parents = branch[:,6]
        roots = list(set(parents)-set(indices))
        for element in roots:
            parents[np.where(parents==element)] = -1
        branch[:,6] = parents
        #branch[np.where(np.isin(branch[:,6], list(roots)))][:,6] = -1
    if isinstance(branch, np.ndarray):
        branch = branch.tolist()
    with open(filename, 'w') as f:
        f.write('')
    with open(filename, 'a') as f:
        for node in branch:
            node = [int(element)  if element%1 == 0.0 else element for element in node]
            line = [str(element) for element in node]
            line = ' '.join(line)
            f.write(line + '\n')
          

def thisNode(index, tree, as_list=True):
    '''Return the node at index from tree as list or np.array.'''
    if isinstance(tree, list):
        tree = np.array(tree)  
    current_node = tree[np.where(tree[:,0]==index)]
    current_node = current_node.reshape(current_node.shape[1])
    if as_list:
        return current_node.tolist()
    else:
        return current_node
    
def findEndpoints(tree, return_node=False):
    '''Return the indices (or nodes) of endpoints.'''
    if isinstance(tree, list):
        tree = np.array(tree)
    indices = tree[:,0]
    parents = tree[:,6]
    endpoints = list(set(indices) - set(parents))
    if return_node:
        endpoint_nodes = []
        for index in endpoints:
            node = thisNode(index, tree)
            endpoint_nodes.append(node)
        return np.array(endpoint_nodes)
    else:
        return np.array(endpoints)


def findRoots(tree, return_node=False):
    '''Return the indices (or nodes) of roots.'''
    if isinstance(tree, list):
        tree = np.array(tree)
    indices = tree[:,0]
    parents = tree[:,6]
    roots = list(set(parents) - set(indices))
    root_nodes = []
    for index in roots:
        nodes = tree[np.where(tree[:,6]==index)].tolist()
        root_nodes.extend(nodes)
    root_nodes_array = np.array(root_nodes)
    if return_node:
        return root_nodes_array
    else:
        return root_nodes_array[:,0]

def findSomaNodes(tree, scale=(1,1,1), return_list=False):
    soma_diameter_threshold = 22*scale[0]
    soma_nodes = []
    soma_node = findRoots(tree, return_node=True)
    current_nodes = [node.tolist() for node in soma_node]
    i = 0
    while True:
        next_nodes = []
        for node in current_nodes:
            if node[5] > soma_diameter_threshold:
                soma_nodes.append(node)
                next_nodes.append(prevNode(node, tree))
        current_nodes = [node.tolist() for children in next_nodes for node in children]
        if current_nodes==[]:
            break

    if return_list:
        return soma_nodes
    else:
        return np.array(soma_nodes)
    
def findSilencedNodes(mainbranch, scale=(1,1,1), max_distance=5, return_list=False):
    silenced_nodes= []
    current_node = findRoots(mainbranch, return_node=True)[0]
    silenced_nodes.append(current_node)
    distance = 0
    while distance < max_distance:
        next_node = prevNode(current_node, mainbranch)[0]
        distance += dist3D(current_node, next_node, scale=scale)
        current_node = next_node
        silenced_nodes.append(current_node)
    if return_list:
        return silenced_nodes
    else:
        return np.array(silenced_nodes)
        
def nextNode(node_or_index, tree):
    '''Return the parent node if the input is a root node return its parent node value (generally -1).'''
    if not isinstance(tree, np.ndarray):
        tree = np.array(tree)    
    if isinstance(node_or_index, (int, float)):
        current_index = node_or_index
    else:
        current_index = node_or_index[0]
    
    
    current_node = thisNode(current_index, tree, as_list=True)
    try:
        next_node = tree[np.where(tree[:,0]==current_node[6])]
        next_node = next_node.reshape(next_node.shape[1])
        return next_node
    except ValueError:
        return current_node[6]

def removeDoubleNodes(tree):
    '''Remove any nodes that appear twice based on their index.'''
    if isinstance(tree, list):
        tree = np.array(tree)
    _, indices = np.unique(tree[:,0], return_index=True)
    return tree[indices]

def prevNode(node_or_index, tree):
    '''Return an array of child nodes.'''
    if not isinstance(tree, np.ndarray):
        tree = np.array(tree)
    if isinstance(node_or_index, (int, float)):
        current_index = node_or_index
    else:
        current_index = node_or_index[0]
    children = tree[np.where(tree[:,6]==current_index)]
    return children

def calculateDistancesTree(tree, scale=(1,1,1), return_sum=True):
    '''Clean up the tree first so that there are no nodes occuring twice!'''
    distances = []
    for node in tree:
        parent = nextNode(node, tree)
        try:
            dist = dist3D(node, parent, scale=scale)
            if not dist==None:
                distances.append(dist)
        except IndexError:
            pass
    
    if return_sum:
        d = sum(distances)
        d = d - (tree[-1][5] * scale[0])
        return d
    else:
        return distances

def dist3D(start, end, scale=(1,1,1)):
    '''Calculate the distance between 2 points (x,y,z) or .swc nodes multiplied by a scaling factor and return the distance'''
    start = np.array(start)
    end = np.array(end)
    if start.size>3:
        start = start[2:5]
    if end.size>3:
        end = end[2:5]
    try:
        dist = np.sqrt(((start[0]-end[0])*scale[0])**2+((start[1]-end[1])*scale[1])**2+((start[2]-end[2])*scale[2])**2)
        return dist
    except IndexError:
        pass
        #print('Index Error in utility.dist3D!!!')
        
def dist3DWithRadius(start, end, scale=(1,1,1)):
    '''Calculate the distance between 2 .swc nodes multiplied by a scaling factor and return the distance minus the radii of the two nodes.'''
    start = np.array(start)
    end = np.array(end)
    dist = np.sqrt(((start[2]-end[2])*scale[0])**2+((start[3]-end[3])*scale[1])**2+((start[4]-end[4])*scale[2])**2)-(start[5]+end[5])*scale[0]
    return dist


def findWindow(node, branch, window_size=4, scale=(1,1,1)):
    if isinstance(node, np.ndarray):
        node = node.tolist()
        #node[1] = 2
    parent_distance = 0
    child_distance = 0
    parent_nodes = []
    child_nodes = []
    # select nodes in distance window_size/2 around node
    current_parent_node = node
    while parent_distance < window_size/2:
        parent_nodes.append(current_parent_node)
        next_parent_node = nextNode(current_parent_node, branch)
        try:
            next_parent_node = next_parent_node.tolist()
        except AttributeError:
            #print('not enough parent nodes!')
            break
        if not isinstance(next_parent_node, list):
            #print('not enough parent nodes!')
            break
        parent_distance += dist3D(current_parent_node, next_parent_node, scale=scale)
        current_parent_node = next_parent_node
    current_child_node = node
    while child_distance < window_size/2:
        child_nodes.append(current_child_node)
        try:
            next_child_node = prevNode(current_child_node, branch).tolist()[0]
        except IndexError:
            #print('not enough child nodes!')
            break            
        if not isinstance(next_child_node, list):
            #print('not enough child nodes!')
            break              
        child_distance += dist3D(current_child_node, next_child_node, scale=scale)
        current_child_node = next_child_node
    child_nodes = child_nodes[1:]
    try:
        window_nodes = np.concatenate((np.array(parent_nodes), np.array(child_nodes)), axis=0)
    except ValueError:
        if parent_nodes==[]:
            window_nodes = child_nodes
        elif child_nodes==[]:
            window_nodes = parent_nodes
        else:
            raise NameError('NoWindowFound')
    return np.array(window_nodes)

def saveTree(filename, tree):
    '''Save a .scw file from a tree (list of lists)'''
    with open(filename, 'w') as f:
        f.write('')
    with open(filename, 'a') as f:
        for branch in tree:
            for node in branch:
                node = [int(element)  if element%1 == 0.0 else element for element in node]
                line = [str(element) for element in node]
                line = ' '.join(line)
                f.write(line + '\n')
    

def vectorAngle3D(a_in, b_in, scale=(1,1,1), angle_in_radians=False):
    '''Calculate and return the angle between two vectors (a [x1, y1, z1,] and b [x2, y2, z2]).'''
    scale = np.array(scale)
    a = a_in*scale
    b = b_in*scale
    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)
    cosine_angle = np.dot(a, b) / (a_norm * b_norm)
    if cosine_angle > 1 or cosine_angle < -1: #Catch a numerical error that leads to cosine_angle > 1 and subsequent miscalculation of angle
        if abs(round(cosine_angle)-cosine_angle) < 1e-15:
            cosine_angle = round(cosine_angle)
    angle = np.arccos(cosine_angle)
    if angle_in_radians:
        return angle
    else:
        return np.degrees(angle)