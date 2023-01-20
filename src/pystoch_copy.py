# pystoch - a Python based code for SGWB mapping from GW interferometer data
# This is the principle code of the PyStoch Package.
# authors: Anirban Ain (anirban.ain@ligo.org), Jishnu Suresh (jishnu.suresh@ligo.org),
#    Sudhagar Suyamprakasam (sudhagar.suyamprakasam@ligo.org) and Sanjit Mitra (sanjit.mitra@ligo.org)
#-------------------------------------------------------------------------------

import healpy as hp
import numpy as np
from scipy import interpolate

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from multiprocessing import Pool, current_process, cpu_count
import re
import sys

if sys.version_info >= (3, 0):
    import configparser
else:
    import ConfigParser as configparser

import time
import h5py
import glob

from detectors import *
#-------------------------------------------------------------------------------
print ('\nPyStoch - a Python based code for SGWB mapping.')

class Pystoch_param:
    #Class to pack all the PyStoch parameters.
    #Initialized by location of the parameter file.
    def __init__(self,param_file):
        param = configparser.ConfigParser()
        param.read(param_file)

        #HealPix nside parameter (power of 2)
        self.nside = param.getint('parameters','nside')
        #dpi of the output map images
        self.map_res = param.getint('parameters','output_map_dpi')

        #minimum and maximum of the frequency range for which maps will be calculated
        self.f_max = param.getfloat('parameters','f_max')
        self.f_min = param.getfloat('parameters','f_min')

        # Spectral index, either from a file or a simple index
        self.Hf_file = param.get('parameters','Hf_file')
        self.alpha = param.getfloat('parameters','alpha')
        # to enable notching option and load the list from a text file
        self.notch_list = param.get('parameters','notch_list')
        self.notching = param.getboolean('parameters','notching')

        #wheather to use pycbc or not.
        self.use_pycbc = param.getboolean('parameters','pycbc')
        #wheather multiprocessing will be used or not
        self.multi_core = param.getboolean('parameters','multiprocessing')
        #if multiprocessing used, then how many threads
        self.n_cpu = param.getint('parameters','multi_threads')
        #wheather injection will be tested or not
        self.injection = param.getboolean('parameters','injection')
        #wheather image of maps (png) will be made or not
        self.draw_maps = param.getboolean('parameters','draw_maps')
        #wheather image of maps (png) for all frequencies will be made or not
        self.draw_all_narrowband_maps = param.getboolean('parameters','draw_all_narrowband_maps')
        #wheather image of maps (png) for given frequencies will be made or not
        self.draw_narrowband_maps = param.get('parameters','draw_narrowband_maps')
        #the frequencies for which maps will be made
        self.draw_maps_f = [float(i) for i in re.findall(r"[-+]?\d*\.\d+|\d+",self.draw_narrowband_maps)]

        #location where maps (in png and hdf5 format) will be saved
        self.output_map_path = param.get('parameters','output_map_location')
        #location where orf seeds and pre-loaded frames will be saved in a fast loadable format
        self.tmp_results_path = param.get('parameters','tmp_results_location')
class Frameset_param:
    #Class to pack all relevant parameters for a set of framesets.
    #Initialized by location of framesets file and name of frameset folder
    def __init__(self,framesets_file,set_name):
        framesets = configparser.ConfigParser()
        framesets.read(framesets_file)

        #location of the frame files
        self.path =  framesets.get(set_name,'path')
        #total number of gwf files
        self.total_frames = framesets.getint(set_name,'total_frames')
        #wheather to process (i.e. use these frames to make maps) or not
        self.process =  framesets.getboolean(set_name,'process')
        #interferometer 1
        self.ifo1 =  framesets.get(set_name,'ifo1')
        #interferometer 2
        self.ifo2 =  framesets.get(set_name,'ifo2')
        #wheather to do overlap correction or not
        self.overlap =  framesets.getboolean(set_name,'overlap')
        #size of frequency bins in the frames
        self.deltaF =  framesets.getfloat(set_name,'deltaF')
        #top and bottom of the frequency range in the frames
        self.fhigh =  framesets.getfloat(set_name,'fhigh')
        self.flow =  framesets.getfloat(set_name,'flow')
        #length of data (seconds) in each frames
        self.segDuration = framesets.getint(set_name,'segDuration')
        #begening and end of data (in GPStime), this is virtual if using folded data
        self.GPSStart =  framesets.getint(set_name,'GPSStart')
        self.GPSEnd =  framesets.getint(set_name,'GPSEnd')
        #self.files =  glob.glob(self.path+'*.gwf')
        self.winFactor = framesets.getfloat(set_name,'winFactor')
        self.w1w2bar = framesets.getfloat(set_name,'w1w2bar')
        self.bias = framesets.getfloat(set_name,'bias')

#pystoch parameters loaded from file
framesets = configparser.ConfigParser()
framesets_file = "../parameters/framesets.ini"
framesets.read(framesets_file)
datasets_full = framesets.sections()

param_default = "../parameters/parameters.ini"
param_file = param_default

datasets = []
datasets_input = []

passed_arg = sys.argv[1:]
if len(passed_arg) > 0:
    for item in passed_arg:
        if item.endswith('.ini'):
            param_file = item
        else:
            datasets_input.append(item)
            if item in datasets_full:
                datasets.append(item)
            else:
                print('Dataset not found:',item)

try:
    parameters = Pystoch_param(param_file)
    print('using parameters from',param_file)
except:
    print('file not found',param_file)
    print('using default file',param_default)
    parameters = Pystoch_param(param_default)

if len(datasets) == 0:
    if len(datasets_input) > 0:
        sys.exit('None of given dataset found. Exiting program.')

    for dataset in datasets_full:
        if framesets.getboolean(dataset,'process'):
            datasets.append(dataset)

if len(datasets) == 0:
    sys.exit("No frameset found")
else:
    print ('\nFollowing framesets will be processed. \n %s \n' % datasets)
    
# import pycbc module only if specified otherwise use pystoch module
if  parameters.use_pycbc:
    import pycbc.detector

# Spectral index H(f)
# Keep the H(f) data locaded once and for all
#Hf_data = np.loadtxt(parameters.Hf_file)
# Call this function only once in the code, otherwise it will read the file everytime
def spectral_index (freq,Hf_data):
    # Linear interpolation
    H_f = interpolate.interp1d(np.log(Hf_data[:,0]), np.log(Hf_data[:,1]))
    H = np.exp(H_f(np.log(freq)))
    return H


def calculate_t_delay_antenna_response_maps_pycbc(*args):
    # calculating the seed ORF matrices using pycbc functions

    # Condition for multiprocessing analysis
    if len(args) == 1:
        ifo1, ifo2, nside, gpst = args[0]
    # Condition for single process analysis
    if len(args) == 4:
        ifo1, ifo2, nside, gpst = args

    #print ('Calculating ORF seed maps for time {:12.2f}. '.format(float(gpst)),current_process().name)
    print ('Calculating ORF seed maps for time {:12.2f}. '.format(float(gpst)))
    sys.stdout.write("\033[F")

    # creating a 'HealPix sized' mesh of angles
    (theta, phi) = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    dec = (np.pi/2) - theta
    ra = phi
    #ra = rotate_map(phi,0,180)

    # using pycbc methods to calculate antenna response
    ifo1_plus, ifo1_cross = ifo1.antenna_pattern(ra, dec, 0, gpst)
    ifo2_plus, ifo2_cross = ifo2.antenna_pattern(ra, dec, 0, gpst)

    # combine the antenna responses
    combined_antenna_resp = (ifo1_plus * ifo2_plus) + (ifo1_cross * ifo2_cross)

    # signal travel time map using pycbc methods
    t_delay = np.vectorize(ifo1.time_delay_from_detector)(ifo2,ra,dec,gpst)

    #these two matrices are the seed matrices for the given GPStime
    return t_delay, combined_antenna_resp


def calculate_t_delay_antenna_response_maps(*args):
    # calculating the seed ORF matrices

    # Condition for multiprocessing analysis
    if len(args) == 1:
        ifo1, ifo2, nside, gpst = args[0]
    # Condition for single process analysis
    if len(args) == 4:
        ifo1, ifo2, nside, gpst = args

    #print ('Calculating ORF seed maps for time {:12.2f}. '.format(float(gpst)),current_process().name)
    print ('Calculating ORF seed maps for time {:12.2f}. '.format(float(gpst)))
    sys.stdout.write("\033[F")

    #pixtheta, pixphi = nside_pix2ang(nside)
    # up or down
    pixtheta, pixphi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    #pixphi = rotate_map(pixphi,0,180)

    #
    combined_antenna_resp, t_delay = combined_antenna_response_t_delay(ifo1, ifo2, np.array([[gpst]]), pixphi, pixtheta)

    #these two matrices are the seed matrices for the given GPStime
    return t_delay, combined_antenna_resp

def load_frame_data(parameters,frame_param,dataset):
    # function to load csd and psd (inv sigma2) and frequency and time Information

    # time and frequency mesh (series) of the csd and psd
    #GPStime = np.arange(frame_param.GPSStart,frame_param.GPSEnd,frame_param.segDuration)
    f_data = np.arange(frame_param.flow, frame_param.fhigh+frame_param.deltaF/2.0, frame_param.deltaF)

    #this string if the firt part of a channel name (H1L1, V1L1 etc.)
    baseline = frame_param.ifo1 + frame_param.ifo2

    # name of the file where the data should be stored (by convert_frames.py)
    frame_data_name = frame_param.path
    #frame_data_name = frame_param.path + '/' + baseline + '_compressed'  +'.hdf5'
    # check if the tmp file exists if yes load csd and psd from there
    try:
        print ('Loading Frame data for {} Hz to {} Hz from following file.'.format(parameters.f_min,parameters.f_max))
        print (frame_data_name)
        frames_preloaded = h5py.File(frame_data_name, "r")
        # if the previous action does not fail then the file exists, loading csd and psd
        csd_data = frames_preloaded['csd'][:]
        sigma_sq_inv_data = frames_preloaded['sigma_sq_inv'][:]
        GPSmid = frames_preloaded['gps_times_mid'][:]
        frames_preloaded.close()
        #print ('frames loaded from preloaded file.')
    except:
        print ('No data found for frameset {}'.format(dataset))

    # using available range unless specified
    if parameters.f_max == 0:
        f_max_used = frame_param.fhigh
    else:
        f_max_used = parameters.f_max

    if parameters.f_min == 0:
        f_min_used = frame_param.flow
    else:
        f_min_used = parameters.f_min

    # creating an array for frequencies which will be processed.
    # also trimming the csd and psd accordingly
    f_all = []
    csd = []
    sigma_sq_inv = []
    for ii,f in enumerate(f_data):
        if f >= parameters.f_min and f <= parameters.f_max:
            f_all.append(f)
            csd.append(list(csd_data[:,ii]))
            sigma_sq_inv.append(list(sigma_sq_inv_data[:,ii]))

    print ('CSD and PSD for {} frequecny bins for range {} Hz to {} Hz loaded.'.format(np.size(f_all),f_all[0],f_all[-1]))
    return GPSmid,f_all,csd,sigma_sq_inv

def seed_matrices(GPSmid,parameters,frame_param,dataset):
    #function to load or calculate the seed seed_matrices
    #time series for which ORF seeds are to be calculated
    #GPStime = np.arange(frame_param.GPSStart,frame_param.GPSEnd,frame_param.segDuration)
    # File where combined antenna response, time delay maps is/will be saved.
    data_file_name = parameters.tmp_results_path+'seed_maps_'+dataset+'_'+np.str(np.int(GPSmid[0]))+'-'+\
    np.str(np.int(GPSmid[-1]))+'-'+np.str(np.int(frame_param.segDuration))+'_nside-'+np.str(parameters.nside)+'.hdf5'

    #check if the orf seeds exist, load them if they do.
    try:
        maps_precalced = h5py.File(data_file_name, "r")
        #if the previous action does not fail then the file exists, loading ORF seeds
        combined_antenna_response = maps_precalced['combined_antenna_response'][:]
        t_delay = maps_precalced['t_delay'][:]
        maps_precalced.close()
        print ('ORF seeds loaded from following file.')
        print (data_file_name)
    #if it does not, calculate combined antenna response, time delay maps
    except:
        start = time.time()
        print ('Calculating ORF seeds for {} time segments (from {:12.2f} to {:12.2f}).'.format(GPSmid.size,float(GPSmid[0]),float(GPSmid[-1])))

        if parameters.use_pycbc:
            # using pycbc objects
            ifo1_pycbc = pycbc.detector.Detector(frame_param.ifo1)
            ifo2_pycbc = pycbc.detector.Detector(frame_param.ifo2)
        if parameters.multi_core:
            # preparing a list argument for multiprocessing
            rr = []
            for tt in GPSmid:
                if parameters.use_pycbc:
                    rr.append([ifo1_pycbc,ifo2_pycbc, parameters.nside, tt])
                else:
                    rr.append([frame_param.ifo1, frame_param.ifo2, parameters.nside, tt])

            # number of cores to be used.
            if parameters.n_cpu == 0:
                # the program will use all available cores (use this only in personal computers)
                pool = Pool(processes=cpu_count())
            else:
                pool = Pool(processes=parameters.n_cpu)
            if parameters.use_pycbc:
                combined_tdelay_and_antenna = np.array(pool.map(calculate_t_delay_antenna_response_maps_pycbc,rr))
            else:
                combined_tdelay_and_antenna = np.array(pool.map(calculate_t_delay_antenna_response_maps,rr))
            pool.close()
            pool.join()

            # separating the outputs
            t_delay = combined_tdelay_and_antenna[:,0,:]
            combined_antenna_response = combined_tdelay_and_antenna[:,1,:]
            # attempt to clean the screen a bit
            print ('\n' * min([x for x in (parameters.n_cpu, cpu_count()) if x != 0]))
        else:
            # using single process to calculate ORF seeds
            t_delay = []
            combined_antenna_response = []
            # simple loop over all GPS times
            for tt in GPSmid:
                if parameters.use_pycbc:
                    t_delay_tt, combined_antenna_response_tt = calculate_t_delay_antenna_response_maps_pycbc(ifo1_pycbc,ifo2_pycbc, parameters.nside, tt)
                else:
                    t_delay_tt, combined_antenna_response_tt = calculate_t_delay_antenna_response_maps(frame_param.ifo1, frame_param.ifo2, parameters.nside, tt)
                t_delay.append(t_delay_tt)
                combined_antenna_response.append(combined_antenna_response_tt)
            t_delay = np.array(t_delay)
            combined_antenna_response = np.array(combined_antenna_response)

        #open a tmp file and save the ORF seeds for future use
        maps_precalced = h5py.File(data_file_name, "w")
        maps_precalced.create_dataset('combined_antenna_response', data = combined_antenna_response)
        maps_precalced.create_dataset('t_delay', data = t_delay)
        maps_precalced.close()
        #this is an estimate of time you will save in future runes
        print ('ORF seeds calculated and saved in {} sec.     '.format( int(time.time() - start)))

    return t_delay, combined_antenna_response

#def calculate_maps(f, segDuration, csd_f, sigma_sq_inv_f,combined_antenna_response,t_delay):
def calculate_maps(*args):
    # calculating the narrowband maps for a frequencybin

    # Condition for multi process analysis
    if len(args) == 1:
        f, segDuration, csd_f, sigma_sq_inv_f,combined_antenna_response,t_delay = args[0]
    # Condition for single process analysis
    if len(args) == 6:
        f, segDuration, csd_f, sigma_sq_inv_f,combined_antenna_response,t_delay = args

    #print ('Calculating maps for frequency {} Hz.       '.format(f),current_process().name)
    print ('Calculating maps for frequency {} Hz.       '.format(f))
    sys.stdout.write("\033[F")

    # the Overlap Reduction Function (ORF) is calculated using the ORF seed matrices
    #gamma = combined_antenna_response * np.exp(2*np.pi*1j*f*t_delay)
    #gamma_star = np.conjugate(gamma)
    gamma_star = combined_antenna_response * np.exp(-1j*(2*np.pi*f)*t_delay)
    map_dirty_f = np.dot(np.array(csd_f)*np.array(sigma_sq_inv_f),gamma_star)*segDuration
    fisher_diag_f = np.dot(sigma_sq_inv_f,np.square(combined_antenna_response))*segDuration**2

    # calculation of maps if there are injections
    #csd_inj = np.dot(gamma,map_inj)*segDuration
    #csd_inj = csd_inj*sigma_sq_inv[:,ii] + csd[:,ii]
    #map_dirty_f_inj = np.dot(csd_inj,gamma_star)*segDuration

    return map_dirty_f,fisher_diag_f

#an empty skymap and an snr map
map_dirty_total = np.zeros(hp.nside2npix(parameters.nside))
fisher_diag_total = np.zeros(hp.nside2npix(parameters.nside))
map_snr_total = np.zeros(hp.nside2npix(parameters.nside))
#Spectrum cue loaded from file
Hf_data = np.loadtxt(parameters.Hf_file)

#Notching cue loaded from file if notching is expected
def make_notch_array(f_all,notching,notch_list):
    f_min = f_all[0]
    f_max = f_all[-1]
    deltaf = f_all[1]-f_all[0]
    numFreqs = np.size(f_all)
    notching_array = np.ones(np.size(f_all))
    if notching:
        notch = np.loadtxt(notch_list)
        freqsToRemove = notch[:,0]
        nBinsToRemove = notch[:,1]
    #notching begins
        for ii in range(0,len(freqsToRemove)):
            if nBinsToRemove[ii]>0:
                index = (freqsToRemove[ii]-f_min)/deltaf  + 1
                if nBinsToRemove[ii]%2 ==0:
                    index_low  = index - nBinsToRemove[ii]/2 + 1
                    index_high = index + nBinsToRemove[ii]/2
                else:
                    index_low  = index - (nBinsToRemove[ii]-1)/2
                    index_high = index + (nBinsToRemove[ii]-1)/2
                if index_low < 1:
                    index_low=1

                if index_high > np.int(numFreqs):
                    index_high = numFreqs

                if index_high >= index_low:
                    notching_array[np.int(index_low)-1:np.int(index_high)]=np.zeros(1)

                if (freqsToRemove[ii] < f_min or freqsToRemove[ii] > f_max):
                    print ('WARNING: Some of the requested frequencies to notch are ouside the search parameter space')

    return notching_array

for dataset in datasets:
    #loop over all framesets

    start = time.time()
    print (" *** processing {} ***".format(dataset))
    #loading parameters of that frameset from input file
    frame_param = Frameset_param(framesets_file,dataset)

    #printing some info to display which frameset and what parameters
    print ("Baseline is {}{}, Total {} frames of {} seconds are available\n"\
    .format(frame_param.ifo1, frame_param.ifo2,frame_param.total_frames,frame_param.segDuration))

    #fetching frame data and ORF seeds
    GPSmid,f_all,csd,sigma_sq_inv = load_frame_data(parameters,frame_param,dataset)
    t_delay, combined_antenna_response = seed_matrices(GPSmid,parameters,frame_param,dataset)


    # SGWB spectrum H(f)
    H = spectral_index(f_all,Hf_data)
    #H = np.power((np.array(f_all) * (1/100.0)),(parameters.alpha-3))

    # use the parallel processing codes if specified
    if parameters.multi_core:
        # an array of arguments to feed the parallel threads
        qq = []
        for ll,f in enumerate(f_all):
    	    qq.append([f,frame_param.segDuration,csd[ll:ll+1],sigma_sq_inv[ll:ll+1],combined_antenna_response,t_delay])

        if parameters.n_cpu == 0:
            #the program will use all available cores (use this only in personal computers)
            pool = Pool(processes=cpu_count())
            chunksize = int(np.ceil(float(np.size(f_all))/cpu_count()))
        else:
            pool = Pool(processes=parameters.n_cpu)
            chunksize = int(np.ceil(float(np.size(f_all))/parameters.n_cpu))

        # starting threads for calculations of dirty maps
        dirty_maps_diag_fisher = np.array(pool.map(calculate_maps,qq,chunksize=chunksize))
        #dirty_maps_diag_fisher = np.array(pool.starmap(calculate_maps,\
        #        [(f,frame_param.segDuration,csd[ll:ll+1],sigma_sq_inv[ll:ll+1],combined_antenna_response,t_delay)\
        #        for ll,f in enumerate(f_all)], chunksize=chunksize))
        pool.close()
        pool.join()

        # extracting the results from the output of the parallel threads
        map_dirty_mat = dirty_maps_diag_fisher[:,0,:]
        fisher_diag_mat = dirty_maps_diag_fisher[:,1,:]

        # free up some memory
        dirty_maps_diag_fisher = None

        # cleaning up the screen
        print ('\n' * min([x for x in (parameters.n_cpu, cpu_count()) if x != 0]))

    # in case of no multiprocessing
    else:
        map_dirty_mat = []
        fisher_diag_mat = []
        # a loop over f, all maps for the frequencies will be calculated one by one
        for ll,f in enumerate(f_all):
            # calculating the dirty maps
            map_dirty_f, fisher_diag_f = calculate_maps(f,frame_param.segDuration, csd[ll:ll+1], sigma_sq_inv[ll:ll+1],\
            combined_antenna_response, t_delay)
            map_dirty_mat.append(map_dirty_f)
            fisher_diag_mat.append(fisher_diag_f)

    notch_array = make_notch_array(f_all,parameters.notching,parameters.notch_list)
    map_dirty_mat = np.squeeze(map_dirty_mat)
    fisher_diag_mat = np.squeeze(fisher_diag_mat)
 
    # summing up all the maps for the broadband results
     #all sky all frequency maps
#     map_dirty_h = np.expand_dims(H,axis=1) * np.squeeze(map_dirty_mat)
#     fisher_diag_h = np.expand_dims(np.square(H),axis=1) * np.squeeze(fisher_diag_mat)
     
    #unnormalised broadband maps
#     map_dirty = np.sum( map_dirty_h ,axis=0)
#     fisher_diag = np.sum(fisher_diag_h, axis=0)
    #ptEst = 2*(frame_param.segDuration)*np.real(np.divide(map_dirty,fisher_diag))
    #sig = np.sqrt(frame_param.segDuration/((frame_param.deltaF)*(frame_param.winFactor))) *frame_param.bias* np.divide(1.0,np.sqrt(np.real(fisher_diag)))
    #map_snr = np.divide(ptEst,sig)

    #map_dirty_inj = np.sum(map_dirty_mat_inj ,axis=0)

    #file to save the outputs in a hdf5 file (for an asaf run)
    maps_mat_output = h5py.File(parameters.output_map_path + 'Map_dirty_'+dataset+'_'+str(parameters.nside)+'_preproc.hdf5', "w")
    maps_mat_output.create_dataset('map_dirty', data = map_dirty_mat) 
    maps_mat_output.create_dataset('fisher_diag', data = fisher_diag_mat)
    maps_mat_output.create_dataset('notch_array', data = notch_array)
    #maps_mat_output.create_dataset('map_snr', data = map_snr)
    maps_mat_output.create_dataset('f_all', data = f_all)
    maps_mat_output.close()

    # adding the map to the main map (will be repeated over datasets)
#     map_dirty_total = map_dirty_total + map_dirty
#     fisher_diag_total = fisher_diag_total + fisher_diag
    

# ptEst_total = 2*(frame_param.segDuration)*np.real(np.divide(map_dirty_total,fisher_diag_total))
# sig_total = np.sqrt(frame_param.segDuration/((frame_param.deltaF)*(frame_param.winFactor))) *frame_param.bias* np.divide(1.0,np.sqrt(np.real(fisher_diag_total)))
# map_snr_total = np.divide(ptEst_total,sig_total)
    


# saving the combined result from all baselines.
# maps_output = h5py.File(parameters.output_map_path + 'Map_dirty'+str(parameters.nside)+'.hdf5', "w")
# maps_output.create_dataset('ptEst_total', data = ptEst_total)
# maps_output.create_dataset('sig_total', data = sig_total)
# maps_output.create_dataset('map_snr_total', data = map_snr_total)
# maps_output.close()

# create png images of the maps
if parameters.draw_maps:
    print ('Map drawing started.')
    # maps are in a mollview projection
    hp.mollview(np.real(np.squeeze(map_dirty_total)),rot=(180,0,0),flip='astro', title="Dirty Map",nest=False)
    plt.savefig(parameters.output_map_path+'Map_dirty-'+str(parameters.nside)+'.png',dpi = parameters.map_res)
    hp.mollview(np.real(np.squeeze(map_snr_total)),rot=(180,0,0),flip='astro', title="SNR Map",nest=False)
    plt.savefig(parameters.output_map_path+'Map_SNR-'+str(parameters.nside)+'.png',dpi = parameters.map_res)

print ('Complete.')
