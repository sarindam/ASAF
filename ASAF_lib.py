#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 14:48:18 2022

@author: arindam
"""

# import h5py
import healpy as hp
import numpy as np
# import matplotlib.pyplot as plt
# import time



def projscatter2(idx, nside, c='red', s=10, marker='o', label=None):
    theta, phi = hp.pix2ang(nside, idx)
    hp.projscatter(theta, phi, c=c, s=s, marker=marker, label=label)
    return

def get_o1neighbours(nside):
    
    npix = hp.nside2npix(nside)
    neighbours_arr = []
    for i in range(npix):
        #Helpy fn. Returns 8 neighbours of the pixel i
        neighbours = hp.get_all_neighbours(nside, i) 
        
        #The next step is necessary as sometimes a pixel will have less than 8 neighbours.
        #Then the missing pixel is labelled as -1. I want to remove that -1 from the list.
        del_arr = []
        for j in range(neighbours.size):
            if neighbours[j] == -1:
                del_arr.append(j)
                
        neighbours = np.delete(neighbours, del_arr)
        
        neighbours_arr.append(neighbours)
        
    return neighbours_arr


def get_o2neighbours(nside):
    
    npix = hp.nside2npix(nside)
    
    o2neighbours = []
    
    for i in range(npix):
        neighbours1 = hp.get_all_neighbours(nside, i)
        
        del_arr = []
        for j in range(neighbours1.size):
            if (neighbours1[j]== -1):
                del_arr.append(j)
                
        neighbours1 = np.delete(neighbours1, del_arr)
        
        
        
        neighbours2 = -np.ones((neighbours1.size,8))
        
        for j in range(neighbours1.size):
            
            neighbours2[j,:] = hp.get_all_neighbours(nside, neighbours1[j])
        
        neighbours2 = np.array(neighbours2, dtype='int')
    
        neighbours2 = np.unique(neighbours2)
        
        
        del_arr = []
        for j in range(neighbours2.size):
            if neighbours2[j] == -1:
                del_arr.append(j)
                
        neighbours2 = np.delete(neighbours2, del_arr)
        
        '''Remove the pixel i from the list of neighbours since it is the pixel 
        whose neighbours are being computed'''
        for k in range(neighbours2.size):
            if neighbours2[k] == i:
                neighbours2 = np.delete(neighbours2, k)
                break
                
        
    o2neighbours.append(neighbours2)

        # peak_idx = np.array(peak_idx, dtype='int')
    return  o2neighbours


def get_o3neighbours(nside):
    
    npix = hp.nside2npix(nside)
    
    o3neighbours = []
    
    for i in range(npix):
        neighbours1 = hp.get_all_neighbours(nside, i)
        
        del_arr = []
        for j in range(neighbours1.size):
            if (neighbours1[j]== -1):
                del_arr.append(j)
                
        neighbours1 = np.delete(neighbours1, del_arr)
        
        
        
        #O2 Neighbours
        neighbours2 = -np.ones((neighbours1.size,8))
        
        for j in range(neighbours1.size):
            
            neighbours2[j,:] = hp.get_all_neighbours(nside, neighbours1[j])
        
        neighbours2 = np.array(neighbours2, dtype='int')
    
        neighbours2 = np.unique(neighbours2)
        
        
        del_arr = []
        for j in range(neighbours2.size):
            if neighbours2[j] == -1:
                del_arr.append(j)
                
        neighbours2 = np.delete(neighbours2, del_arr)
        
        
        
        #O3 neighbours 
        neighbours3 = -np.ones((neighbours2.size, 8))
        
        for j in range(neighbours2.size):
            neighbours3[j,:] = hp.get_all_neighbours(nside, neighbours2[j])
        
        #Convert neighbours3 to a numpy array
        neighbours3 = np.array(neighbours3, dtype='int')
        
        #Delete redundant pixels
        neighbours3 = np.unique(neighbours3)
        
        #Delete entries with -1. These occur when a pixel has less than 8 neighbours
        del_arr = []
        for j in range(neighbours3.size):
            if neighbours3[j] == -1:
                del_arr.append(j)
                
        neighbours3 = np.delete(neighbours3, del_arr)
        
        #Delete centre pixel whose neighbours you are considering 
        for k in range(neighbours3.size):
            if neighbours3[k] == i:
                neighbours3 = np.delete(neighbours3, k)
                break
            
        o3neighbours.append(neighbours3)
        
    return  o3neighbours


def get_o4neighbours(nside):
    
    npix = hp.nside2npix(nside)
    
    o4neighbours = []
    
    for i in range(npix):
        neighbours1 = hp.get_all_neighbours(nside, i)
        
        del_arr = []
        for j in range(neighbours1.size):
            if (neighbours1[j]== -1):
                del_arr.append(j)
                
        neighbours1 = np.delete(neighbours1, del_arr)
        
        
        
        
        neighbours2 = -np.ones((neighbours1.size,8))
        
        for j in range(neighbours1.size):
            
            neighbours2[j,:] = hp.get_all_neighbours(nside, neighbours1[j])
        
        neighbours2 = np.array(neighbours2, dtype='int')
    
        neighbours2 = np.unique(neighbours2)
        del_arr = []
        for j in range(neighbours2.size):
            if neighbours2[j] == -1:
                del_arr.append(j)
                
        neighbours2 = np.delete(neighbours2, del_arr)
        
        
        
        
        neighbours3 = -np.ones((neighbours2.size, 8))
        
        for j in range(neighbours2.size):
            neighbours3[j,:] = hp.get_all_neighbours(nside, neighbours2[j])
        
        #Convert neighbours3 to a numpy array
        neighbours3 = np.array(neighbours3, dtype='int')
        
        #Delete redundant pixels
        neighbours3 = np.unique(neighbours3)
        
        #Delete entries with -1. These occur when a pixel has less than 8 neighbours
        del_arr = []
        for j in range(neighbours3.size):
            if neighbours3[j] == -1:
                del_arr.append(j)
                
        neighbours3 = np.delete(neighbours3, del_arr)
        
        
        
        
        
        neighbours4 = -np.ones((neighbours3.size, 8))
        
        for j in range(neighbours3.size):
            neighbours4[j,:] = hp.get_all_neighbours(nside, neighbours3[j])
        
        #Convert neighbours4 to a numpy array
        neighbours4 = np.array(neighbours4, dtype='int')
        
        #Delete redundant pixels
        neighbours4 = np.unique(neighbours4)
        
        #Delete entries with -1. These occur when a pixel has less than 8 neighbours
        del_arr = []
        for j in range(neighbours4.size):
            if neighbours4[j] == -1:
                del_arr.append(j)
                
        neighbours4 = np.delete(neighbours4, del_arr)
        
        #Delete centre pixel whose neighbours you are considering
        for k in range(neighbours4.size):
            if neighbours4[k] == i:
                neighbours4 = np.delete(neighbours4, k)
                break
                
        o4neighbours.append(neighbours4)

    return  o4neighbours




def get_peaks_new(SNRmap, neighbours_arr):
    '''Returns the indices of peaks of ONE map
    For every pixel, this code checks its neighbours 
    and counts it as a peak if its SNR is greater than ALL its neighbours.
    
    First calculate all the neighbours according to the req. nside and order
    and pass that array as input to this function
    
    SNRmap: Shape = npix X 1
    neighbour_arr: A list of arrays. ith array in the list gives all the neighbours of
                   the ith pixel
    '''
    
    npix = SNRmap.size
    # nside = hp.npix2nside(npix)
    
    if len(neighbours_arr) != npix:
        raise Exception("neighbours_arr not of the right size")
        
    peak_idx = []
    for i in range(npix):

        neighbours = neighbours_arr[i]
        neighbour_SNR = SNRmap[neighbours]
    
        if (SNRmap[i] >= neighbour_SNR).all():
            peak_idx = np.append(peak_idx, i)
     
    peak_idx = peak_idx.astype('int')
    return peak_idx



    
    
    
