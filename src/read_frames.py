# read_frames - a code to read GW interferometer data and generate a parameter fil
# This code is a part of the PyStoch Package.
# authors: Anirban Ain (anirban.ain@ligo.org), Jishnu Suresh (jishnu.suresh@ligo.org),
#    Sudhagar Suyamprakasam (sudhagar.suyamprakasam@ligo.org) and Sanjit Mitra (sanjit.mitra@ligo.org)
#-------------------------------------------------------------------------------

import numpy as np
import glob
import os
import datetime
import sys

if sys.version_info >= (3, 0):
    import configparser
else:
    import ConfigParser as configparser


print ("\nA code to check the available frames and prepare parameter list for PyStoch\n")

# Reading the location where gwf framefiles are stored from the main PyStoch parameter file.
param = configparser.ConfigParser()
param.read("../parameters/parameters.ini")
frames_path = param.get('parameters','frames_location')
use_pycbc = param.getboolean('parameters','pycbc')

# if pycbc is not installed then this is useful to use gwpy modules
if  use_pycbc:
    import pycbc.frame
else:
    from gwpy.timeseries import TimeSeriesDict

# frameset: subdirectories in the frames directory. there should be frames from a single run.
# the framesets.ini file should have parameters for all framesets
parameter_file = open("../parameters/framesets.ini","w")

# this is a list for all framesets (sub-diectories) in the frames folder.
directory_list = list()

def read_channels(frame, chnl_list):
    # a wrapper function that uses gwpy to get data from frames
    data = TimeSeriesDict.read(frame, chnl_list)
    output = []
    for chnl in chnl_list:
        output.append(data[chnl].value)
    return output

# this loop is over all elements inside the frames directory
for root, dirs, files in os.walk(frames_path, topdown=False):
    # this loop is over all subdirectories
    for name in dirs:
        # full path to a subdirectories
        dir_full = os.path.join(root, name)
        # writing the path to the file
        #parameter_file.write("[{}]".format(name))
        # adding the directory in the directry list
        directory_list.append(dir_full)
        print ("Frame directory :", dir_full)
        # this is a list of all files in the frameset
        files = [ff for ff in os.listdir(dir_full) if os.path.isfile(os.path.join(dir_full, ff))]
        # this is a list of full paths to all files in the frameset
        files_full_path = sorted(glob.glob(dir_full+'/*.gwf'))
        if sum(np.shape(files_full_path)) == 0:
            print ('There are no frames in {}\n'.format(dir_full))
            continue
        # number of frames (i.e. files) in that frameset
        nn = np.shape(files_full_path)[0]
        print (nn,'Frames available for processing in this directory.')
        # the interferometers are identified from the begening of frame names.
        # it is very iportant to maintain this convention correctly.
        ifo1_name =  files[0][0]+'1'
        ifo2_name =  files[0][1]+'1'
        # baseline is just ifo1 to ifo2
        baseline = ifo1_name+ifo2_name

        try:
            # These are the parameters required to process a frameset
            chnl_list = [baseline+':deltaF',baseline+':fhigh', baseline+':flow',baseline+':foldedSegmentDuration',baseline+':GPStime',baseline+':winFactor',baseline+':w1w2bar',baseline+':bias' ]
            if use_pycbc:
                # useing pycbc if specified
                deltaF,fhigh,flow,segDuration,GPSStart,winFactor,w1w2bar,bias = np.array(pycbc.frame.read_frame(files_full_path[0],chnl_list))
            else:
                # using gwpy if pycbc is not available
                deltaF,fhigh,flow,segDuration,GPSStart,winFactor,w1w2bar,bias =  read_channels(files_full_path[0],chnl_list)

            GPSEnd = GPSStart + segDuration*nn
            #print (deltaF,fhigh,flow,segDuration,GPSStart,GPSEnd
            print ("Frequency Information: deltaF {}, flow {}, fhigh {}.".format(float(deltaF), float(flow), float(fhigh)))
            print ('Time Information: segDuration {}, GPSStart {}, GPSEnd {}.'.format(int(segDuration), int(GPSStart), int(GPSEnd)))

            #these are the frameset parameters which will be written in file for all framesets
            parameter_file.write("[{}]".format(name))
            # path of the frameset
            parameter_file.write("\npath:  {}".format(dir_full))
            # number of frames in that frameset
            parameter_file.write("\ntotal_frames: {}".format(nn))
            # wheather to include this frameset in processing
            parameter_file.write("\nprocess:  True")
            # first interferometer of the baseline
            parameter_file.write("\nifo1:  {}".format(ifo1_name))
            # second interferometer of the baseline
            parameter_file.write("\nifo2:  {}".format(ifo2_name))
            # whether overlap correction has to be done or not
            parameter_file.write("\noverlap:  False")
            # Frequency Information
            parameter_file.write("\ndeltaF:  {}".format(float(deltaF)))
            parameter_file.write("\nfhigh:  {}".format(float(fhigh)))
            parameter_file.write("\nflow:  {}".format(float(flow)))
            #normalization information
            parameter_file.write("\nwinFactor: {}".format(float(winFactor)))
            parameter_file.write("\nw1w2bar: {}".format(float(w1w2bar)))
            parameter_file.write("\nbias: {}".format(float(bias)))
            # Time Information (For FSID this is virtual, may not match the run)
            parameter_file.write("\nsegDuration:  {}".format(int(segDuration)))
            parameter_file.write("\nGPSStart:  {}".format(int(GPSStart)))
            parameter_file.write("\nGPSEnd:  {}\n".format(int(GPSEnd)))
        except:
            # if an expected channel is not found then the above try will fail.
            print (chnl_list)
            print (files_full_path[0])
            print ('One or more expected channel(s) not found.\n')

        print ("Baseline: {}\n".format(baseline))

parameter_file.close()

print ("Frame list prepared. Check frame details, edit parameter file then run PyStoch.\n")
