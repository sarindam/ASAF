#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 14:48:18 2022

@author: arindam
"""

# import h5py
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
# from scipy.optimize import curve_fit
# import time



def projscatter2(idx, nside, **kwargs):
    theta, phi = hp.pix2ang(nside, idx)
    hp.projscatter(theta, phi, **kwargs)
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
    
    First find the healpix indices of all neighbours of each pixel
    and pass that array as input to this function
    
    SNRmap: Shape = npix
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




def pvalSNRthres_old(SNRmap, bins, pval, plot=False, intg_left=True):
    counts, bin_edges = np.histogram(SNRmap, density=True, bins=bins)
    

    if plot==True:
        fig, ax  = plt.subplots()
        ax.hist(SNRmap, bins=bins, density=True)
        
    
    if intg_left:
        intg = 0
        i = 0
        while intg <= (1-pval):
            x = bin_edges[i+1] - bin_edges[i]
            fx = counts[i]
            
            intg = intg + fx * x
            
            i = i+1
            
        out = bin_edges[i]   
        
        # print('p-val m1=', intg)    
    
    else:
        intg = 0
        
        for j in range(bin_edges.size-1, 0, -1):
            x = bin_edges[j] - bin_edges[j-1]
            fx = counts[j-1]
            intg = intg + fx*x
            # print(bin_edges[j], bin_edges[j-1], intg)
    
            if intg > pval:
                # print('pval 0.05 at SNR', bin_edges[j-1])
                out = bin_edges[j-1]
                break
    
        # print('p-val m1=', intg)
            
    return round(out,3)

def pvalSNRthres_new(SNRmap, bins, pval, plot=False, intg_left=True):
    counts, bin_edges = np.histogram(SNRmap, density=True, bins=bins)
    
    x = (bin_edges[1:] + bin_edges[:-1])/2

    spline = interp1d(x, counts, kind='cubic')
    # print('Done')
    
    x_arr = np.linspace(x[0], x[-1], 1000)
    y = spline(x_arr)
    
    # tot_intg = quad(lambda x: spline(x), x[0], x[-1])[0]
    # print('integral check =', tot_intg)
    if plot==True:
        fig, ax  = plt.subplots()
        ax.hist(SNRmap, bins=bins, density=True)
        ax.scatter(x, counts)
        ax.plot(x_arr, y, c='r')
        # plt.scatter(x_arr, y, s=10, c='r')
    
    if intg_left:
        
        intg = 0
        i = 0
        while intg<= 1 - pval:
            intg = intg + quad(lambda x: spline(x), x_arr[i], x_arr[i+1])[0]
            # print(x_arr[i], intg)
            i = i + 1
        
        # print('i = ', i)
        out = x_arr[i]       
        # print('integral =', intg)
        # print(out)
        
        # print('p-val m2=', intg)
        
    else:
        
        intg = 0
        i=-1
        while intg <= pval:
            intg = intg + quad(lambda x: spline(x), x_arr[i-1], x_arr[i])[0]
            i = i - 1
        
        out = x_arr[i-1]
        
        # print('p-val m2=', intg)
        
    return round(out,3)



def inject_map2(fisher, fisher_diag, pix_inj, snr_inj, nmaps, seed=10, ret_noise=False):
    '''
    More efficient code to inject sources into multiple maps at the same time. 
    The dirty map is given by X = \Gamma . \hat{P} + n
    n is the random noise generated from a multivariate Gaussian distribution with pixel to 
    pixel correlation given by the Fisher matrix \Gamma.

    \hat{P} is zero for all pixels except the pixels where we are injecting sources.
    These pixels have SNRs given by snr_inj.
    Parameters
    ----------
    fisher : The full Fisher matrix for the frequency we want to make injected maps.
             We use the same Fisher matrix for all realisations of noise. 
             Shape = npix X npix
    fisher_diag : Diagonalised fisher matrix. 
                We use the same Fisher diagonal for all realisations of noise. Shape = npix
    pix_inj : If we want to do only one injection and the same one in each map then
              Shape = 1, else Shape = Ninj X Nmaps. MUST be a numpy array
    snr_inj : If we want to do only one injection and the same one in each map then
        Shape = 1, else Shape = Ninj X Nmaps. MUST be a numpy array
    nmaps : no. of maps we want to generate
    ret_noise: Whether to return noise maps or not

    Raises
    ------
    Exception
        pix_inj and snr_inj should have the shape
        

    Returns
    -------
    dirty_SNR_mat : Shape = npix X nmaps
    noise_SNR_mat: Shape = npix X nmaps

    '''
    
    npix = fisher.shape[0]

    
    if pix_inj.shape != snr_inj.shape:
        raise Exception('Shapes of pix_inj and snr_inj do not match')
        
    sqrt_fisher_diag = np.sqrt(fisher_diag)


    
    u,s,v=np.linalg.svd(np.real(fisher))
    rng = np.random.default_rng(seed)
    x = rng.normal(0,1,size=(npix, nmaps))
    
    #Regularisation
    s_module_reg=np.zeros_like(s)
    s_module_reg[:]=s[:]
    s_module_reg_inv=np.zeros_like(s)
    s_module_reg_inv[:]=s[:]
    neg = (np.sum(u.T * v, axis=1) < 0) & (s > 0)
    if np.any(neg):  
        s_module_reg[neg] = 0.
        s_module_reg_inv[neg] = np.inf
    
    
    noise = np.dot(u,np.dot(np.diag(s_module_reg**0.5),x))
    
    
    Phat = np.zeros((npix, nmaps))
    

    if (pix_inj.ndim == 1):
        ninj = pix_inj.size
        
        pix_inj = pix_inj.reshape(ninj,1)
        pix_inj_mat = np.dot(pix_inj, np.ones((1,nmaps)))
        pix_inj_mat = pix_inj_mat.astype('int')
        
        snr_inj = snr_inj.reshape(ninj,1)
        snr_inj_mat = np.dot(snr_inj, np.ones((1,nmaps)))
        # if ninj==1:
        #     pix_inj_mat = pix_inj * np.ones((ninj, nmaps), dtype=int)
        #     snr_inj_mat = snr_inj * np.ones((ninj, nmaps))
        # elif ninf>1:
            
    elif (pix_inj.ndim > 1):
        # ninj = pix_inj.shape[0]
        pix_inj_mat = pix_inj
        snr_inj_mat = snr_inj
    else:
        raise Exception('Idx Injection matrix of wrong shape')
        
    for j in range(nmaps):
        for i,p in enumerate(pix_inj_mat[:,j]):
            Phat[p,j] = snr_inj_mat[i,j] 
            
        Phat[:,j] = Phat[:,j]/sqrt_fisher_diag
    
    X_mat = np.dot(fisher, Phat) + noise
    X_mat = X_mat.real
    
    dirty_SNR_mat = np.zeros((npix, nmaps))
    noise_SNR_mat = np.zeros_like(dirty_SNR_mat)
    for j in range(nmaps):
            dirty_SNR_mat[:,j] = X_mat[:,j]/sqrt_fisher_diag
            
            noise_SNR_mat[:,j] = noise[:,j]/sqrt_fisher_diag
            
    
    if ret_noise==False:
        return dirty_SNR_mat
    else:
        return dirty_SNR_mat, noise_SNR_mat





def get_noise_maps(fisher, fisher_diag, nmaps, nside):

    # sqrt_fisher_diag = np.sqrt(fisher_diag[0])
    npix = hp.nside2npix(nside)
    
    
    u,s,v = np.linalg.svd(np.real(fisher))
    rng = np.random.default_rng(10)
    x = rng.normal(0,1, size=(npix, nmaps))

    #Regularisation
    s_module_reg=np.zeros_like(s)
    s_module_reg[:]=s[:]
    s_module_reg_inv=np.zeros_like(s)
    s_module_reg_inv[:]=s[:]
    neg = (np.sum(u.T * v, axis=1) < 0) & (s > 0)
    if np.any(neg):  
        s_module_reg[neg] = 0.
        s_module_reg_inv[neg] = np.inf
    
    
    noise = np.dot(u,np.dot(np.diag(s_module_reg**0.5),x))
    
    noise_SNR = np.zeros((npix, nmaps))
    for j in range(nmaps):
        noise_SNR[:,j] = noise[:,j]/np.sqrt(fisher_diag)
        
    return noise_SNR

