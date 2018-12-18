#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright (c) 2017-2018 Max Hess, Alvaro Gomariz, ETH Zurich
'''


import numpy as np
import utility
import SimpleITK as sitk
import skimage.io
from skimage import color, io
from skimage.viewer import CollectionViewer
import matplotlib.pyplot as plt
import os
import csv
from PIL import Image
from tqdm import trange

def getSomaPosition(filename_swc_or_tree):
    if isinstance(filename_swc_or_tree, str):
        tree = utility.readSWC(filename_swc)
    else:
        tree = filename_swc_or_tree
    roots = tree[np.where(tree[:,1]==0)]
    seeds = roots[:,2:5].astype(int).tolist()
    return seeds

def getSomaPositionSingleSeed(filename_swc_or_tree):
    if isinstance(filename_swc_or_tree, str):
        tree = utility.readSWC(filename_swc)
    else:
        tree = filename_swc_or_tree
    roots = tree[np.where(tree[:,6]==-1)]
    seeds = roots[:,2:5].astype(int).tolist()
    return seeds

def calculateSomaVolume(filename_tif, 
                        seeds, 
                        scale=(1,1,1), 
                        visualize=False, 
                        outfolder=None, 
                        flip_image=True,
                        numberOfIterations=1,
                        multiplier=2.5,
                        initialNeighborhoodRadius=1,
                        replaceValue=1,
                        vectorRadius=(3,3,3)):
    #Load Image and run region growing segmentation
    img_array = skimage.io.imread(filename_tif)
    if flip_image:
        img_array = np.flip(img_array, axis=1) #flip the image along the x axis as the APP2 output is mirrored along the x axis
    img = sitk.GetImageFromArray(img_array)
    seg = sitk.Image(img.GetSize(), sitk.sitkUInt8)
    seg.CopyInformation(img)
    for element in seeds:
        seg[element] = 1
    seg_conf = sitk.ConfidenceConnected(img, seedList=seeds,
                                        numberOfIterations=numberOfIterations,
                                        multiplier=multiplier,
                                        initialNeighborhoodRadius=initialNeighborhoodRadius,
                                        replaceValue=replaceValue)
    
    #Cleanup
    vectorRadius = vectorRadius
    kernel = sitk.sitkBall
    seg_clean = sitk.BinaryMorphologicalOpening(seg_conf,
                                                vectorRadius,
                                                kernel)
    seg_array = sitk.GetArrayFromImage(seg_clean)
    
    #Calculate Soma volume
    v_voxel = scale[0] * scale[1] * scale[2]
    n_soma_pixels = seg_array.sum()
    soma_volume = n_soma_pixels * v_voxel
    
    #Visualization
    if visualize:
        outimage = []
        for i in trange(len(img_array)):
            img = img_array[i]
            mask = seg_array[i]
            rows, cols = img.shape
            if i == 13:
                a=3
            
            color_mask = np.zeros((rows, cols, 3))
            color_mask[:,:,0] = mask
            #color_mask[170:270, 40:120] = [0, 1, 0] # Green block
            #color_mask[200:350, 200:350] = [0, 0, 1] # Blue block
            
            img_color = np.dstack((img, img, img))
            img_hsv = color.rgb2hsv(img_color)
            color_mask_hsv = color.rgb2hsv(color_mask)
            
            img_hsv[..., 0] = color_mask_hsv[..., 0]
            img_hsv[..., 1] = color_mask_hsv[..., 1] * 0.6
            img_masked = color.hsv2rgb(img_hsv)*255
            img_final = np.array(img_masked, dtype='uint8')
            outimage.append(Image.fromarray(img_final, mode='RGB'))
            #print(i)
        
        outimage_name = filename_tif.split('/')[-1][:-4] + '_soma.gif'
        outimage[0].save(outfolder + outimage_name, format="GIF", save_all=True, append_images=outimage[1:])
        #print('done!')
    return soma_volume


def somavolume(filename_swc_or_tree='data/trees/alm.swc',
               filename_tif='data/trees/alm.tif',
               scale=(0.223,0.223,0.3),
               visualize=True,
               outfolder_gif='data/trees/',
               soma_volume_singleseed=True,
               numberOfIterations=2,
               multiplier=2.5,
               initialNeighborhoodRadius=3,
               replaceValue=1,
               vectorRadius=(3,3,3)):
    
    if soma_volume_singleseed:
        seeds = getSomaPositionSingleSeed(filename_swc_or_tree)
    else:
        seeds = getSomaPosition(filename_swc_or_tree)
    #print(seeds)
    #seed = (33,77,13)
    volume = calculateSomaVolume(filename_tif, 
                                 seeds, 
                                 scale=scale, 
                                 visualize=visualize, 
                                 outfolder=outfolder_gif,
                                 numberOfIterations=numberOfIterations,
                                 multiplier=multiplier,
                                 initialNeighborhoodRadius=initialNeighborhoodRadius,
                                 replaceValue=replaceValue,
                                 vectorRadius=vectorRadius)
    return volume
    
    
if __name__ == '__main__':
    volume=somavolume()
    print(volume)