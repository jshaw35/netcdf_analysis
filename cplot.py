#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 18:27:32 2019

@author: jonahshaw
"""

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import xarray as xr
#import os

# assume [dates/alt, lat, lon] format of var_arr
def cplot(var_arr, arr_lat, arr_lon, lev = 41):

    # Calculate the range over all data and create a linear stepset for the colorbar
    cmin_p = np.nanmin(var_arr)
    cmax_p = np.nanmax(var_arr)

    cmap_p = 'bwr'
    nlevels = lev
    cmap2 = plt.get_cmap(cmap_p)
    
    if cmin_p == cmax_p:
       cmax_p = cmax_p + 0.00001

    levels = np.linspace(cmin_p,cmax_p,nlevels)
    
#    datap[datap < cmin_p] = cmin_p  # Selects the data within set bounds, likely removing nans
#    datap[datap > cmax_p] = cmax_p
    
    # Calculate the correct subplot dimensions
    timesteps = np.shape(var_arr)[0]
    col = np.int(timesteps**(0.5)); rem = timesteps % (timesteps**(0.5));
    if rem > 0: col += 1

    fig = plt.figure() # Create the figure object
    
    for i, data in enumerate(var_arr): # Iterate over the first array dimension, time.
        i += 1
        sp = fig.add_subplot(col, col, i, projection=ccrs.PlateCarree()) # add subplot with correct projection
        mpbl = sp.contourf(arr_lon, arr_lat, data, levels, cmap=cmap2, transform=ccrs.PlateCarree()) # plot the monthly data

    axs = fig.get_axes()
#    cbar = fig.colorbar(mpbl, ax=axs, orientation="horizontal")

    return fig, mpbl

# Input netCDF filename and list of variables to extract. Returns dictionary object keyed appropriately.
def CDFextract(filelist, varnames):
    var_dict = {}
    fnt = (len(filelist),) # file number tuple

    # Instantiate nan array of appropriate size for each variable of interest.
    # I am not entirely sure why this is necessary
    first_file = xr.open_dataset(filelist[0]) # Load first file for shared via xarray. f0.keys() gives all variables
    for var in varnames:
        temp_shape = fnt + np.shape(first_file.variables[var]) # Number of files + original variable shape 
        temp_arr = np.zeros(temp_shape) * np.nan
        var_dict[var] = temp_arr
    first_file.close()

    # Move netCDF data into nan arrays    
    for i, files in enumerate(filelist):                          # For each output file
        filedata = xr.open_dataset(files)                         # Open the file
        for temp_var in varnames:                                 # For each variable of interest
            vardata = filedata.variables[temp_var][:]             # pull out data for each variable
            var_dict[temp_var][i] = vardata                       # And store that data within the correct index of the corresponding dictionary value
        filedata.close()       
    
    return var_dict