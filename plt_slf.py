#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:46:45 2019

@author: jonahks
Plot annual avg slf as a function of isobar. Based on ts_attempt2.py

There is currently a probelm with averaging! Do not forget!
"""

import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import mpl_toolkits
import matplotlib
#matplotlib.use('Agg')                   # I am unsure of this
import matplotlib.pyplot as plt
import os
from datetime import datetime

bool_online = False

cwd = os.getcwd()

dir_ext = cwd
day = datetime.now()
daystamp = day.strftime("/%m%d%Y")
tstamp = day.strftime("%H%M%S")
output_dir = dir_ext + daystamp
#print(output_dir)

# create output directory if it does not already exist. Update daily for organizational purposes
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

filepaths0 = []
#filepaths00 = os.listdir(rootdirs[0])
filepaths00 = os.listdir(cwd)  # cwd for now, will need to change

filepaths0 = [x for x in filepaths00 if 'cam.h0.2004' in x]   # List of files within the root directory containing string "cam.h0."
#print(filepaths0)

ind_vars = ['lat', 'lon', 'lev'] # Values shared by all output files
dep_vars = ['TGCLDCWP', 'TGCLDLWP', 'TGCLDIWP', 'LWC', 'IWC', 'TS'] # dependent variables (unique for each output)

f0 = xr.open_dataset(filepaths0[0]) # Load via xarray
#print("f0 keys: ", f0.keys())  # like an ncdump
ivar_dict = {}
for var in ind_vars:
    ivar_dict[var] = f0.variables[var][:]
#    print(var, ' ', np.shape(temp_arr))
f0.close()

lat_cam = ivar_dict['lat'][:]
lon_cam = ivar_dict['lon'][:] # selects all longitude 'lon' values from the netcdf instance
lev_cam = ivar_dict['lev'][:]

# JShaw's plan
# for each variable in dep_vars:
dvar_dict = {} # empty dictionary to store data 
fnt = (len(filepaths0),) # file number tuple
#
#for var in dep_vars:
#    # get variable shape
#    temp_arr = np.zeros((len(filepaths0),len(lat_cam),len(lon_cam))) * np.nan
    # Create empty array of nans with dimensions # of output files # internal dims to variable shape = (#files) + np.shape(varfile)
    # then go through files themselves and store the data appropriately, will need to keep track of file #
for var in dep_vars:
    temp_shape = fnt + np.shape(f0.variables[var])
    temp_arr = np.zeros(temp_shape) * np.nan
    dvar_dict[var] = temp_arr
#    print(var, ' ', np.shape(temp_arr))
f0.close()

for i, files in enumerate(filepaths0):                         # For each output file
    filedata = xr.open_dataset(files)
    for temp_var in dep_vars:                                # For each variable of interest
        vardata = filedata.variables[temp_var][:]             # pull out data for each variable
        dvar_dict[temp_var][i] = vardata                     # And store that data within the correct index of the corresponding dictionary value
    filedata.close()                                        # close file before opening the next one
        
# Create list of variables. Go through and get shape for each variable, create empty nan array for each variable by shape. Use xarray correctly, maybe a dictionary is a good object
weights0 = np.cos(np.deg2rad(lat_cam))  # Create list of weights by latitude. But this is only a 1-d array, and the data is...
weights_cam = weights0 / np.nansum(weights0)  # Normalizes weights, so that final value will make sense

vars_cam_tavg = np.zeros((len(filepaths0),len(lat_cam),len(lon_cam))) * np.nan  # Create an empty numpy array of nans for storing output (month, lat, lon)
vars_cam_t = np.zeros(len(filepaths0)) * np.nan # For global temperature averages

"""
This following order of operations is important so that each latitude slice is weighted by it's area independent of the data sampling. 
We use nanmean and nansum to avoid deflating composite values when nans are present
1. Compute average variable value along longitude lines, excluding nans
2. Weight each of these averages by area and sum for a global average independent of sampling efficiency
"""

annualavg = np.nanmean(dvar_dict['TS'], axis = 0)  # average by zeroeth index, which is time

#quit()
for i, files in enumerate(filepaths0):
#    onview = xr.open_dataset(files)
#    ts_cam = onview.variables['TS'][:]
#    onview.close()
#    dictdat = dvar_dict['TS'][i]  # Get temperature map by month
#    print('Onview: '); print(np.shape(ts_cam)); print(ts_cam[0])
#    print('Dict: '); print(np.shape(dictdat)); print(dictdat[0])
#    dictdat_2d = dictdat[0,:,:]

#    lat_avgs2 = np.nanmean(dictdat_2d, axis = 1) # Step 1
#    print('Onview: '); print(np.shape(lat_avgs)); print(lat_avgs)
#    print('Dict: '); print(np.shape(lat_avgs2)); print(lat_avgs2)
#    break
    ts_cam = dvar_dict['TS'][i]  # Get temperature map by month
    ts_cam_2d = ts_cam[0,:,:]  # Eliminate level from surface temperature shape
    lat_avgs = np.nanmean(ts_cam_2d, axis = 1) # Step 1
    vars_cam_tavg[i,:,:] = ts_cam_2d            # Store the month that has just been processed
    vars_cam_t[i] = np.nansum(lat_avgs*weights_cam) # Step 2, then store by month

#print(vars_cam_t - 273.15)
vars_cam_yt = np.nanmean(vars_cam_t) # Annual single average
vars_cam_yavg = np.nanmean(vars_cam_tavg, axis = 0) #lat-lon annual average

cmin_p = np.nanmin(vars_cam_tavg)
cmax_p = np.nanmax(vars_cam_tavg)
#print("Min: ", str(cmin_p), " Max: ", str(cmax_p))

cmap_p = 'bwr'

nlevels = 41

datap = np.copy(vars_cam_tavg)
datap = datap - 273.15 # Convert to C

if cmin_p == cmax_p:
   cmax_p = cmax_p + 0.00001

datap[datap < cmin_p] = cmin_p  # Selects the data within set bounds, likely removing nans
datap[datap > cmax_p] = cmax_p

cmap2 = plt.get_cmap(cmap_p)

nlevels = 41 # Probably relic from 0C to 40C range earlier
levels1 = np.linspace(cmin_p,cmax_p,nlevels)

newfig = plt.figure(3)
ax3 = newfig.add_subplot(111, projection=ccrs.PlateCarree())


ct = ax3.contourf(lon_cam, lat_cam, vars_cam_tavg[i,:,:], levels1, cmap=cmap2, transform=ccrs.PlateCarree())
tempstr = str(np.round(vars_cam_yt - 273.15, 2)) + " C"
ax3.set_title(tempstr)

if bool_online: ax3.coastlines()
#ax3.drawcountries()

cbar = newfig.colorbar(ct, orientation="horizontal",shrink=0.4)
cbar.ax.tick_params(labelsize=7) 

cbar.ax.set_xlabel('$^\circ$C')

filename = "/carto" + tstamp 
newfig.savefig(output_dir + filename  + '.pdf')
newfig.clf()

fig4 = plt.figure(4) #subplots(nrows=3, ncols=4)

#for i in range(len(filepaths0)):
for i, month in enumerate(dvar_dict['TS']):
    i+=1
    sp = fig4.add_subplot(3, 4, i, projection=ccrs.PlateCarree()) # add subplot with correct projection
    i+=-1
#    mpbl = sp.contourf(lon_cam, lat_cam, vars_cam_tavg[i,:,:], levels1, cmap=cmap2, transform=ccrs.PlateCarree()) # plot the monthly data
#    mpbl = sp.contourf(lon_cam, lat_cam, dvar_dict['TS'][i,0,:,:], levels1, cmap=cmap2, transform=ccrs.PlateCarree()) # plot the monthly data
    mpbl = sp.contourf(lon_cam, lat_cam, month[0,:,:], levels1, cmap=cmap2, transform=ccrs.PlateCarree()) # plot the monthly data

#    axi = ax_flat[i]
    temp_string = filepaths0[i][-10:-3] + ' ' + str(round(vars_cam_t[i] - 273, 2)) + ' C'
#    axi.set_title(temp_string,fontsize=10)
    sp.set_title(temp_string,fontsize=10)
    if bool_online: sp.coastlines()
    
ax4 = fig4.get_axes()

cbar = fig4.colorbar(mpbl, ax=ax4, orientation="horizontal",shrink=0.4)
cbar.set_label('$^\circ$C')
cbar.ax.tick_params(labelsize=4) # Could round to 3 digits instead

filename = "/cartomon" + tstamp 
fig4.savefig(output_dir + filename  + '.pdf')
fig4.clf()

fig5 = plt.figure(5)

cmin_p = np.nanmin(dvar_dict['LWC'])
cmax_p = np.nanmax(dvar_dict['LWC'])
#print("Min: ", str(cmin_p), " Max: ", str(cmax_p))

cmap_p = 'bwr'

nlevels = 41

#datap = np.copy(vars_cam_tavg)
#datap = datap - 273.15 # Convert to C

if cmin_p == cmax_p:
   cmax_p = cmax_p + 0.00001

datap[datap < cmin_p] = cmin_p  # Selects the data within set bounds, likely removing nans
datap[datap > cmax_p] = cmax_p

cmap2 = plt.get_cmap(cmap_p)

nlevels = 41 # Probably relic from 0C to 40C range earlier
levels1 = np.linspace(cmin_p,cmax_p,nlevels)


#temp_var = ivar_dict['lev']
right_levs = [x for x in np.array(ivar_dict['lev']) if (np.float(x) > 300 and np.float(x) < 700)]
col = np.int(len(right_levs)**0.5) + 1
#print("Num right levs: ", len(right_levs))
#print("Fig 5 dim: ", col)
#print("WP shape: ", np.shape(dvar_dict['LWC']))
temp_arr = np.nanmean(dvar_dict['LWC'], axis = 0)[0,:,:,:]
# This code is ugly and I don't like it!

#filepaths0 = [x for x in filepaths00 if 'cam.h0.2004' in x] 
j = 1
for i, lev in enumerate(np.array(ivar_dict['lev'])):
    if (np.float(lev) > 300 and np.float(lev) < 700):  # select the pressure range of interest
#        print(lev)
        sp = fig5.add_subplot(col, col, j, projection=ccrs.PlateCarree()) # add subplot with correct projection
        j += 1
        mpbl = sp.contourf(lon_cam, lat_cam, temp_arr[i,:,:], levels1, cmap=cmap2, transform=ccrs.PlateCarree()) # plot the monthly data
    
    #    axi = ax_flat[i]
        temp_string = str(np.round(lev,1)) + ' mbar'
    #    axi.set_title(temp_string,fontsize=10)
        sp.set_title(temp_string,fontsize=10)
        # plt the variable globally? add to a map?
ax5 = fig5.get_axes()

cbar = fig5.colorbar(mpbl, ax=ax5, orientation="horizontal",shrink=0.4)
cbar.set_label('$^\circ$C')
cbar.ax.tick_params(labelsize=4) # Could round to 3 digits instead

filename = "/slf_alt" + tstamp 
fig5.savefig(output_dir + filename  + '.pdf')
fig5.clf()

plt.show()

def var_avg(arr_var): # general function to average across monthly data sets
    # assume three dimensions, could specify dim later. Maybe always assume format (blah, blah,..., lat, lon)
    # iterate over time dimension
    ind_lon = -1 # assume longitude is the last dimension
    arr_shape = np.shape(arr_var)
    for i in range(arr_shape[0]):
        np.nanmean(arr_var, axis = ind_lon) # average across longitudes
    return False
    