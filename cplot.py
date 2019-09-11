#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 18:27:32 2019

@author: jonahshaw
"""

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

# assume [dates, one, levs, lat, lon] format

def cplot(var_arr, lat = False, lon = False):
    if lat == False: arr_lat = np.linspace(0, 178.125, 96) # then use the normal latitude array
        
    if lon == False: arr_lon = np.linspace(0, 257.5, 144) # use the normal longitude array
    fig = plt.figure()
    
    cmin_p = np.nanmin(var_arr)
    cmax_p = np.nanmax(var_arr)

    cmap_p = 'bwr'
    nlevels = 41
    cmap2 = plt.get_cmap(cmap_p)
    
    if cmin_p == cmax_p:
       cmax_p = cmax_p + 0.00001

    levels = np.linspace(cmin_p,cmax_p,nlevels)
    
    datap[datap < cmin_p] = cmin_p  # Selects the data within set bounds, likely removing nans
    datap[datap > cmax_p] = cmax_p

    col = np.shape(var_arr)[0] # assuming the zeroeth index is correct
    
    for i, data in enumerate(var_arr):
        i += 1
        sp = fig.add_subplot(col, col, i, projection=ccrs.PlateCarree()) # add subplot with correct projection

        ax = fig.get_axes()
        mpbl = sp.contourf(arr_lon, arr_lat, i, levels, cmap=cmap2, transform=ccrs.PlateCarree()) # plot the monthly data
# need to figure out lon and lat inputs

    fig.subplots_adjust(hspace = 0.55)
    cbar = fig.colorbar(mpbl, ax=ax, orientation="horizontal",shrink=0.4)
#    cbar.set_label('LWC')
#    cbar.ax.tick_params(labelsize=4) # Could round to 3 digits instead


    return fig