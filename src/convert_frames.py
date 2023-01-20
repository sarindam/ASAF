# convert_frames - a Python based code to convert gwf files into hdf5 files
# This code is a part of the PyStoch Package.
# authors: Anirban Ain (anirban.ain@ligo.org), Jishnu Suresh (jishnu.suresh@ligo.org),
#    Sudhagar Suyamprakasam (sudhagar.suyamprakasam@ligo.org) and Sanjit Mitra (sanjit.mitra@ligo.org)
#-------------------------------------------------------------------------------

import os
import glob
import h5py
import numpy as np

import time
import sys

if sys.version_info >= (3, 0):
    import configparser
else:
    import ConfigParser as configparser

param = configparser.ConfigParser()
param_file = "../parameters/parameters.ini"
param.read(param_file)

use_pycbc = param.getboolean('parameters','pycbc')

if  use_pycbc:
    import pycbc.frame
else:
    from gwpy.timeseries import TimeSeriesDict

#frame parameters loaded from file
framesets = configparser.ConfigParser()
framesets_file = "../parameters/framesets.ini"
framesets.read(framesets_file)

#full list of available framesets
datasets_full = framesets.sections()

for dataset in datasets_full:
    if framesets.getboolean(dataset,'process'):
        frames_path = framesets.get(dataset,'path')
        files = sorted(glob.glob(frames_path+'/*.gwf'))
        file_names = [ff for ff in os.listdir(frames_path) if os.path.isfile(os.path.join(frames_path, ff))]
        #print files
        baseline = file_names[0][0]+'1'+file_names[0][1]+'1'
        csd_cnl_rl = baseline + ':' + param.get('parameters','csd_channel_real')
        csd_cnl_im = baseline + ':' + param.get('parameters','csd_channel_imaginary')
        single_psd_cnl = param.getboolean('parameters','psd_combined')
        time_cnl = baseline + ':GPStime'
        if single_psd_cnl:
            psd_cnl = baseline + ':' + param.get('parameters','psd_channel')
        else:
            psd_cnl1 = file_names[0][0]+'1'+':' + param.get('parameters','psd_channel')
            psd_cnl2 = file_names[0][1]+'1'+':' + param.get('parameters','psd_channel')

        #start = time.time()

        csd = []
        sigma_sq_inv = []
        gps_times = []
        #print csd_cnl_rl, csd_cnl_im, psd_cnl
        print ("\nDataset: {}, total {} frames".format(dataset, framesets.getint(dataset,'total_frames')))
        frame_data_name = frames_path + '/' + baseline + '_compressed'  +'.hdf5'

        #try:
        #    frames_preloaded = h5py.File(frame_data_name, "r")
        #    print ('Loading Frame data from saved hdf5 file.'
        #    csd = frames_preloaded['csd'][:]
        #    sigma_sq_inv = frames_preloaded['sigma_sq_inv'][:]
        #    frames_preloaded.close()
        #except:
        print ('Reading and converting frames.')
        start = time.time()
        if  use_pycbc:
            for ii,frame in enumerate(files):
                temp1, temp2, temp4 = pycbc.frame.read_frame(frame,[csd_cnl_rl,csd_cnl_im,time_cnl])
                csd.append(np.array(temp1+(1j*temp2)))
                gps_times.append(temp4)
                if single_psd_cnl:
                    temp3 = pycbc.frame.read_frame(frame,[psd_cnl])
                    sigma_sq_inv.append(np.array(temp3))
                else:
                    temp3, temp4 = pycbc.frame.read_frame(frame,[psd_cnl1,psd_cnl2])
                    sigma_sq_inv.append(np.array(np.reciprocal(np.multiply(temp3,temp4))))

                print (ii,'frames loaded.')
                sys.stdout.write("\033[F")
        else:
            for ii,frame in enumerate(files):
                if single_psd_cnl:
                    data = TimeSeriesDict.read(frame,[csd_cnl_rl,csd_cnl_im ,psd_cnl,time_cnl])
                    csd.append(data[csd_cnl_rl].value+(1j*data[csd_cnl_im].value))
                    sigma_sq_inv.append(data[psd_cnl].value)
                    gps_times.append(data[time_cnl].value)
                else:
                    data = TimeSeriesDict.read(frame,[csd_cnl_rl,csd_cnl_im ,psd_cnl1,psd_cnl2,time_cnl])
                    csd.append(data[csd_cnl_rl].value+(1j*data[csd_cnl_im].value))
                    sigma_sq_inv.append(np.reciprocal(np.multiply(data[psd_cnl1].value, data[psd_cnl2].value)))
                    gps_times.append(data[time_cnl].value)

                print (ii,'frames loaded.')
                sys.stdout.write("\033[F")

        # shift GPS time to the midpoint of the segment
        half_seg = float(framesets.get(dataset,'segDuration'))/2.0
        gps_times_mid = [x + half_seg for x in gps_times]

        frames_preloaded = h5py.File(frame_data_name, "w")
        frames_preloaded.create_dataset('csd', data = csd)
        frames_preloaded.create_dataset('sigma_sq_inv', data = sigma_sq_inv)
        frames_preloaded.create_dataset('gps_times_mid', data = gps_times_mid)
        frames_preloaded.close()
        #display how long it took to load the frames.
        print ('frames loaded and saved in {} seconds.'.format( int(time.time() - start)))
