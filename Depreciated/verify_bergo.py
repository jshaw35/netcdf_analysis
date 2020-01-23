#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 14112019

@author: jonahks
Validate that modifications to micro_mg2_0, micro_mg_cam, and hetfrz_classnuc_cam do not
change model output when inp and wbf tags are set to zero.
"""

import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import mpl_toolkits
import matplotlib
#matplotlib.use('Agg')                   # I am unsure of this
import matplotlib.pyplot as plt
import os
import sys
from datetime import datetime
import pandas as pd

from cplot import cplot
from cplot import CDFextract


def main():
    '''
    User inputs:
    '''
    
    bool_online = True
    bool_save = True
    
#    dir_data = '~/Documents/rsyncd/' #'~/Documents/cesm1data' #'/Users/jonahshaw/Desktop/Fulbright/Bergen files'  # Where the data files should be accessed from
    dir_data = '/home/jonahks/Documents/rsyncd'
    dir_out = dir_data                                            # Where output folders should be created
#    file_ext = 'cam.h0.2004'                                      # str to separate files of interest within the working directory
    noresmraw_ext = '1437'
    noresmmod_ext = 'slf1_inp1'
    
    ind_vars = ['lat', 'lon', 'lev'] # Values shared by all output files
    dep_vars = ['TS', 'MEANSLF_ISOTM', 'CLDTOT_ISOTM'] # dependent variables (unique for each output)
    
    """
    END USER INPUTS
    """
    
    print(os.getcwd())
#    os.chdir('~/')

    # Move to directory where data is stored. Quit if the directory does not exist.
    try:
        os.chdir(dir_data)
        
    except OSError:
        print('Directory ' + dir_data + ' does not exist. Exiting...')
        sys.exit(1)
    
    exit

    # Create folder to store output figures
    day = datetime.now()
    daystamp = day.strftime("/%m%d%Y")
    tstamp = day.strftime("%H%M%S")
    output_dir = dir_out + daystamp
    
    # Create output directory if it does not already exist.
    # Update daily to avoid clutter/stay organized
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    filepaths00 = os.listdir(dir_data)                       # All objects in the data directory
    filepaths_raw = [x for x in filepaths00 if noresmraw_ext in x]   # Objects within the root directory containing with the right extension
    filepaths_raw = sorted(filepaths_raw)                          # Sort to make sure order is correct
    
    filepaths_mod = [x for x in filepaths00 if noresmmod_ext in x]   # Objects within the root directory containing with the right extension
    filepaths_mod = sorted(filepaths_mod)                          # Sort to make sure order is correct
    
#    print(filepaths00)
#    print(filepaths_raw, filepaths_mod)

    var = ['MNUCCDOhet','BERGO', 'CLOUD']
    
    _rawtest = xr.open_dataset(filepaths_raw[1])
    _modtest = xr.open_dataset(filepaths_mod[1])
    
    berg_raw = _rawtest[var[2]][0,:,:,:]
    berg_mod = _modtest[var[2]][0,:,:,:]
    
    _comp1 = (berg_raw == berg_mod)  # Array of booleans indicating equivalence

    print(np.all(_comp1).values)    # This checks whether everything is True or not.
#    print(_comp1)
#    print(np.nanmean(_rawtest['BERGO']))
    
#    print(np.max(_rawtest['BERGO'][0,:,:,:]))
#    print(np.max(_modtest['BERGO'][0,:,:,:]))

#    print(np.median(_rawtest['BERGO'][0,:,:,:]))
#    print(np.median(_modtest['BERGO'][0,:,:,:]))

# average columns
    
#    ivar_dict = CDFextract([filepaths0[0]], ind_vars)  # Create dictionary for independent variables
#    dvar_dict = CDFextract(filepaths0, dep_vars)       # Create dictionary for dependent variables
#    
#    lat_cam = ivar_dict['lat'][0,:] # Pull latitude array
#    lon_cam = ivar_dict['lon'][0,:] # Pull longitude array
#    
#    # Create SLF array
#    slf = np.nanmean((dvar_dict['MEANSLF_ISOTM']/dvar_dict['CLDTOT_ISOTM']), axis = 0)[0,:,:,:]
#    arr_temp = [-40, -35, -30, -25, -20, -15, -10, -5, 0]
#    
#    print(np.shape(slf))
#    
#    weights0 = np.cos(np.deg2rad(lat_cam))  # Create list of weights by latitude. But this is only a 1-d array, and the data is...
#    weights_cam = weights0 / np.nansum(weights0)  # Normalizes weights, so that final value will make sense
#    
#    # Average by longitude
#    slf_lon = np.nanmean(slf, axis = 2)
#    
#    numbins = 2
#    lat_bins = create_bins(50, 20, numbins)
##    print(lat_cam)
#    lat_binned = np.zeros((2, 9)) * np.nan  # first index is the bin, the second is divided between date and the accumulated weight
#    
#    # For each isoterm, bin by latitude
#    for i, bins in enumerate(lat_bins):                     # For each pre-determined latitude bin (indexed by i) 
#        for j, isos in enumerate(slf_lon):                  # For each isotherm (indexed by j)
#            temp_avg = 0                                    # 
#            temp_weight = 0
#            for k, lats in enumerate(isos):                 # For each latitude indexed by latitude
#                ind_bin = find_bin(lat_cam[k], lat_bins)
##                print(ind_bin); print(lat_cam[k])
#                if ind_bin == i:                            # If that latitude goes in the jth bin, put it their weighted by latitude and keep track of the weights
#                    temp_avg += weights0[k] * lats
#                    temp_weight += weights0[k]
#            print(i); print(j); print(temp_avg); print(temp_weight);
#            lat_binned[i,j] = temp_avg / temp_weight        # After going through each value, divide out the accumulated weighting and stored
#            
#    print(lat_binned)
#    
#    # Alternate bining
##
##    headers = ['lat', 'slf']
##    new_bins = [-90.1,50,70,90.1]
##    for j,isos in enumerate(slf_lon):
##        temp = np.transpose([lat_cam, isos])
##        data_no_headers = pd.read_csv("hubble_data_no_headers.csv", names = headers)
#    
#def create_bins(lower_bound, width, quantity):
#    """ create_bins returns an equal-width (distance) partitioning. 
#        It returns an ascending list of tuples, representing the intervals.
#        A tuple bins[i], i.e. (bins[i][0], bins[i][1])  with i > 0 
#        and i < quantity, satisfies the following conditions:
#            (1) bins[i][0] + width == bins[i][1]
#            (2) bins[i-1][0] + width == bins[i][0] and
#                bins[i-1][1] + width == bins[i][1]
#    """
#    
#    bins = []
#    for low in range(lower_bound, 
#                     lower_bound + quantity*width - 1, width): # I changed this because it made too many bins
#        bins.append((low, low+width))
#    return bins
#
#def find_bin(value, bins):
#    """ bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
#        binning returns the smallest index i of bins so that
#        bin[i][0] <= value < bin[i][1]
#    """
#    
#    for i in range(0, len(bins)):
#        if bins[i][0] < value <= bins[i][1]: # inclusive on the high side
#            return i
#    return -1            

if __name__ == '__main__':
    main()