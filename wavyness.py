#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright (c) 2017-2018 Max Hess, Alvaro Gomariz, ETH Zurich
'''

import numpy as np
import utility
import matplotlib.pyplot as plt
from scipy import interpolate


def findWindow(node, branch, window_size=4, scale=(0.22, 0.22, 0.3)):
    """
    Given a `node` on a `branch`, return a window of size `window_size` [um].

    Parameters
    ----------
    node : np.ndarray
        Node around which to return a window.
    branch : np.ndarray
        Parent branch of `node`.
    window_size : float [um]
        Maximum distance of the returned window along the branch.
    scale : tuple of floats
        x, y and z scales of the images underlying the analysis.

    Returns
    -------
    window_nodes : np.ndarray
        Array of nodes belonging to the window.
    """
    if isinstance(node, np.ndarray):
        node = node.tolist()
        node[1] = 2
    parent_distance = 0
    child_distance = 0
    parent_nodes = []
    child_nodes = []
    # select nodes in distance window_size/2 around node
    current_parent_node = node
    while parent_distance < window_size/2:
        parent_nodes.append(current_parent_node)
        next_parent_node = utility.nextNode(current_parent_node, branch).tolist()
        if not isinstance(next_parent_node, list):
            #print('not enough parent nodes!')
            break
        parent_distance += utility.dist3D(current_parent_node, next_parent_node, scale=scale)
        current_parent_node = next_parent_node
    current_child_node = node
    while child_distance < window_size/2:
        child_nodes.append(current_child_node)
        next_child_node = utility.prevNode(current_child_node, branch).tolist()
        if not next_child_node == []:
            next_child_node = next_child_node[0]
        if not isinstance(next_child_node, list):
            #print('not enough child nodes!')
            break
        child_distance += utility.dist3D(current_child_node, next_child_node, scale=scale)
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


def calculateAnglesWithLinearRegression(node, branch, window_size=3.2, scale=(0.22, 0.22, 0.3), visualize=True,
                                        fixed_node=False):
    """

    Parameters
    ----------
    node : np.ndarray
        Node for which to calculate the angle.
    branch : np.ndarray
        Parent branch of `node`.
    window_size : float [um]
        Size of the window on which to do the linear regression.
    scale : tuple of floats
        x, y and z scales of the images underlying the analysis.
    visualize : bool
        Visualization of linear regressoin (for debugging).
    fixed_node : bool
        Whether the linear regression needs to go through `node`.

    Returns
    -------
    angle : float
        Local angle of `node`.
    """
    if isinstance(node, np.ndarray):
        node = node.tolist()
    parent_distance = 0
    child_distance = 0
    parent_nodes = []
    child_nodes = []
    
    #select parent nodes in given window
    current_parent_node = node
    while parent_distance < window_size/2:
        parent_nodes.append(current_parent_node)
        try:
            next_parent_node = utility.nextNode(current_parent_node, branch).tolist()
        except:
            return 180
        if not isinstance(next_parent_node, list):
            #print('not enough parent nodes!')
            return 180
        parent_distance += utility.dist3D(current_parent_node, next_parent_node, scale=scale)
        current_parent_node = next_parent_node
    
    #select child nodes in given window    
    current_child_node = node
    while child_distance < window_size/2:
        child_nodes.append(current_child_node)
        next_child_node = utility.prevNode(current_child_node, branch).tolist()
        if not next_child_node ==[]:
            next_child_node = next_child_node[0]
        
        if not isinstance(next_child_node, list):
            #print('not enough child nodes!')
            return 180
        try:
            child_distance += utility.dist3D(current_child_node, next_child_node, scale=scale)
            current_child_node = next_child_node
        except:
            return 180
    
    #take the coordinates from the nodes
    parent_nodes = np.array(parent_nodes)
    child_nodes = np.array(child_nodes)
    parent_points = parent_nodes[:, 2:5]
    child_points = child_nodes[:, 2:5]
    
    #calculate the mean of the points
    parent_mean = parent_points.mean(axis=0)
    child_mean = child_points.mean(axis=0)
    
    #calculate svd's
    parent_uu, parent_dd, parent_vv = np.linalg.svd(parent_points - parent_mean)
    child_uu, child_dd, child_vv = np.linalg.svd(child_points - child_mean)
    
    parent_uu_fixednode, parent_dd_fixednode, parent_vv_fixednode = np.linalg.svd(parent_points - parent_points[0])
    child_uu_fixednode, child_dd_fixednode, child_vv_fixednode = np. linalg.svd(child_points - child_points[0])
    
    #calculate vectors and angle
    parent_vector = parent_vv[0]
    if utility.dist3D(parent_points[0]+parent_vector, parent_points[-1]) > utility.dist3D(parent_points[0]-parent_vector, parent_points[-1]):
        parent_vector *= -1
    child_vector = child_vv[0]
    if utility.dist3D(child_points[0]+child_vector, child_points[-1]) > utility.dist3D(child_points[0]-child_vector, child_points[-1]):
        child_vector *= -1
    angle = utility.vectorAngle3D(parent_vector, child_vector)
    
    parent_vector_fixednode = parent_vv_fixednode[0]
    if utility.dist3D(parent_points[0]+parent_vector_fixednode, parent_points[-1]) > utility.dist3D(parent_points[0]-parent_vector_fixednode, parent_points[-1]):
        parent_vector_fixednode *= -1
    child_vector_fixednode = child_vv_fixednode[0]
    if utility.dist3D(child_points[0]+child_vector_fixednode, child_points[-1]) > utility.dist3D(child_points[0]-child_vector_fixednode, child_points[-1]):
        child_vector_fixednode *= -1    
    angle_fixednode = utility.vectorAngle3D(parent_vector_fixednode, child_vector_fixednode)
    
    
    
    #visualization
    if visualize:
        linspace = np.reshape(np.linspace(-10,10,2), (2,1))
        parent_line = parent_vector * linspace
        child_line = child_vector * linspace
        parent_line += parent_mean
        child_line += child_mean
        
        linspace_fixednode = np.reshape(np.linspace(-20,0,2), (2,1))
        parent_line_fixednode = parent_vector_fixednode * linspace_fixednode
        child_line_fixednode = child_vector_fixednode * linspace_fixednode
        parent_line_fixednode += parent_points[0]
        child_line_fixednode += child_points[0]

        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d as m3d
        lins = np.reshape(np.linspace(0,1,2), (2,1))
        a = parent_points - parent_mean
        a_line = parent_vv[0] * lins
        b = child_points - child_mean
        b_line = child_vv[0] * lins
        c = parent_points - parent_points[0]
        c_line = parent_vv_fixednode[0] * lins
        d = child_points - child_points[0]
        d_line = child_vv_fixednode[0] * lins
        
        ax = m3d.Axes3D(plt.figure())
        ax.scatter3D(*parent_points.T, color='red')
        ax.quiver(parent_points[0][0], parent_points[0][1], parent_points[0][2], parent_vector[0], parent_vector[1], parent_vector[2], color='red')
        ax.quiver(parent_points[0][0], parent_points[0][1], parent_points[0][2], parent_vv_fixednode[0][0], parent_vv_fixednode[0][1], parent_vv_fixednode[0][2], color='orangered')
        ax.scatter3D(*child_points.T, color='blue')
        ax.quiver(child_points[0][0], child_points[0][1], child_points[0][2], child_vector[0], child_vector[1], child_vector[2], color='blue')
        ax.quiver(child_points[0][0], child_points[0][1], child_points[0][2], child_vv_fixednode[0][0], child_vv_fixednode[0][1], child_vv_fixednode[0][2], color='cyan')
        string = 'angles: ' + str(angle) + '/' + str(angle_fixednode)
        ax.text(child_points[0][0]+1, child_points[0][1], child_points[0][2], s=string)
        plt.show()
    if fixed_node:
        return angle_fixednode
    else:
        return angle
        

def wavyness(mainbranch,
             angle_threshold=145,
             window_size_linear_regression=4.0,
             window_size_maximum_supression=4.0,
             n_colors=10,
             scale=(0.223,0.223,0.3),
             fix_node=False,
             plot_cdf=False):
    """
    Detects sharp bends along a branch.

    Parameters
    ----------
    mainbranch : np.ndarray
        Mainbranch along with to count the sharp bends (kinks).
    angle_threshold : float or int
        Threshold value for a node to be considered a sharp bend.
    window_size_linear_regression : float [um]
        Size of the window on which to perform linear regression for angle calculation.
    window_size_maximum_supression : float [um]
        Size of the window for non-maximum supression (preventing double counting of sharp bends).
    n_colors : int
        Number of colors for visualization of angles.
    scale : tuple of floats
        x, y and z scales of the images underlying the analysis.
    fix_node : bool
        If True, linear regressions need to go through the location of the node on which the angle is calculated.
    plot_cdf : bool
        Whether to plot the cdf of resulting angles (for debugging).

    Returns
    -------
    kinks_count : int
        Number of sharp bends (kinks) detected.
    annotated_mainbranch : np.ndarray
        Mainbranch with kinks and outgrowth events annotated.
    kinks : np.ndarray
        List of nodes that were detected as sharp bends.
    tree : np.ndarray
        Mainbranch with classification according to local angles (for visualization).
    """
    
    if plot_cdf:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.grid(True)
        ax.set_title('Cumulative distribution')
        ax.set_xlabel('Angle (deg)')
        ax.set_ylabel('Percentage')
        #ax.axis([50,cutoff_angle,0,1])    
    
    tree = np.array(mainbranch)
    tree[:,5] = 0.5
    annotated_mainbranch=tree.copy()
    angles = np.zeros(len(tree))
    for i in range(len(tree)):
        angles[i] = calculateAnglesWithLinearRegression(tree[i], tree, window_size=window_size_linear_regression, visualize=False, fixed_node=fix_node)
    angles = np.reshape(angles, (len(angles), 1))
    angles[np.where(np.isnan(angles))] = 180
    angles[:20]=180   # set the first and last 20 nodes to 180 as they generally don't correspond to real kinks
    angles[-20:]=180
    data = np.concatenate((tree, angles), axis=1)
    sample_numbers = data[:,0]
    
    kinks_count = 0
    kinks = []
    while min(data[:,7]) < angle_threshold:
        annotated_mainbranch[np.argwhere(annotated_mainbranch[:,0]==data[data[:,7].argmin()][0])[0][0],5]=3
        kinks.append(data[data[:,7].argmin()].tolist())
        tree[data[:,7].argmin(), 5] = 3
        kinks_count += 1
        w = findWindow(data[data[:,7].argmin()][:7], tree, window_size=window_size_maximum_supression, scale=scale)
        indices = np.argwhere(np.isin(data[:,0], w[:,0])).reshape(len(w))
        for index in indices:
            data[index,7] = 180
    
    kinks = np.array(kinks)
    try:
        kinks[:,5] = 3
    except IndexError:
        pass
    
    if plot_cdf:
        n, bins, patches = ax.hist(kinks[:,7], bins=10000, normed=1, histtype='step',
                                   cumulative=True, label='neuronname')
        patches[0].set_xy(patches[0].get_xy()[:-1])
        ax.legend(loc='center left')
        

    #print(angles.min())    
    m = interpolate.interp1d([0, 180], [1, n_colors])
    normalized_angles = np.round(m(angles)).reshape(len(angles))
    tree[:, 1] = normalized_angles
    
    
    if plot_cdf:
        plt.show()
    return kinks_count, annotated_mainbranch, kinks, tree
        
if __name__ == '__main__':
    wavyness()