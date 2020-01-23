#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 10:17:38 2019

@author: jonahshaw
Test cplot.py functions cplot and CDFextract
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

from cplot import cplot
from cplot import CDFextract

'''
User inputs:
'''

bool_online = True
bool_save = True

dir_data = '/uio/kant/geo-gjest-u1/jonahks/Documents/cesm1data' #'~/Documents/cesm1data' #'/Users/jonahshaw/Desktop/Fulbright/Bergen files'  # Where the data files should be accessed from
#dir_data = '/false/false'
dir_out = dir_data                                            # Where output folders should be created
file_ext = 'cam.h0.2004'                                      # str to separate files of interest within the working directory

ind_vars = ['lat', 'lon', 'lev'] # Values shared by all output files
dep_vars = ['TGCLDCWP', 'TGCLDLWP', 'TGCLDIWP', 'LWC', 'IWC', 'TS', 'MEANSLF_ISOTM', 'CLDTOT_ISOTM'] # dependent variables (unique for each output)

"""
END USER INPUTS
"""
tosave = []
# Move to directory where data is stored. Quit if the directory does not exist.
try:
    os.chdir(dir_data)
    
except OSError:
    print('Directory ' + dir_data + ' does not exist. Exiting...')
    sys.exit(1)

# Create folder to store output figures
day = datetime.now()
daystamp = day.strftime("/%m%d%Y/")
tstamp = day.strftime("%H%M%S")
output_dir = dir_out + daystamp

# Create output directory if it does not already exist.
# Update daily to avoid clutter/stay organized
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
filepaths00 = os.listdir(dir_data)                       # All objects in the data directory
filepaths0 = [x for x in filepaths00 if file_ext in x]   # Objects within the root directory containing with the right extension
filepaths0 = sorted(filepaths0)                          # Sort to make sure order is correct

ivar_dict = CDFextract([filepaths0[0]], ind_vars)  # Create dictionary for independent variables
dvar_dict = CDFextract(filepaths0, dep_vars)       # Create dictionary for dependent variables

lat_cam = ivar_dict['lat'][0,:] # Pull latitude array
lon_cam = ivar_dict['lon'][0,:] # Pull longitude array

varplot, mapp = cplot(dvar_dict['TS'][:,0,:,:], lat_cam, lon_cam)  # Plot monthly surface temperature (in one step!)

varplot.subplots_adjust(hspace = 0.25)
tsax = varplot.get_axes()

# Add colorbar and adjust
cbar = varplot.colorbar(mapp, ax=tsax, orientation="horizontal")
cbar.set_label("TS")
cbar.ax.tick_params(labelsize=8)

# Set the title of each subplot
for i, aax in enumerate(tsax):
    aax.set_title(filepaths0[i][-10:-3])
    if bool_online: aax.coastlines()
    
tosave.append([varplot, 'ts'])
# Manually move the colorbar down to fix positioning
#cbarxy = np.array(tsax[-1].get_position())
#left = cbarxy[0][0]; bottom = cbarxy[0][1]; width = cbarxy[1][0] - cbarxy[0][0]; height = cbarxy[1][1] - cbarxy[0][1]
#newcbar = [left, bottom - 0.15, width, height]
#tsax[-1].set_position(newcbar)

figs.append([varplot, "surface temps"])

# Now have some fun with SLFs
slf = np.nanmean((dvar_dict['MEANSLF_ISOTM']/dvar_dict['CLDTOT_ISOTM']), axis = 0)[0,:,:,:]
arr_temp = [-40, -35, -30, -25, -20, -15, -10, -5, 0]

newfig, mapp = cplot(slf, lat_cam, lon_cam)
newfig.subplots_adjust(hspace = 0.50)

tsax = newfig.get_axes()
for i, aax in enumerate(tsax):  # Ignore the last axis, which is the colorbar
    aax.set_title(str(arr_temp[i]) + '$^\circ$C')

    if bool_online: aax.coastlines()
    
# Add colorbar and adjust
cbar = varplot.colorbar(mapp, ax=tsax, orientation="horizontal")
cbar.set_label("SLF")
cbar.ax.tick_params(labelsize=8)

tosave.append([newfig, 'slf'])

if bool_save:
    for i in tosave:
        temp_fig = i[0]; temp_name = i[1];
        filename = "/" + temp_name + tstamp 
        temp_fig.savefig(output_dir + filename  + '.pdf')
        temp_fig.clf()

# And we're done, yay!
else: plt.show()
