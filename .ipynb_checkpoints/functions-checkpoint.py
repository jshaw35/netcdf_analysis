from imports import *

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy as cy

# Functions from Abisko 2019 examples
def masked_average(xa:xr.DataArray,
                   dim=None,
                   weights:xr.DataArray=None,
                   mask:xr.DataArray=None):
    """
    This function will average
    :param xa: dataArray
    :param dim: dimension or list of dimensions. e.g. 'lat' or ['lat','lon','time']
    :param weights: weights (as xarray)
    :param mask: mask (as xarray), True where values to be masked.
    :return: masked average xarray
    """
    #lest make a copy of the xa
    xa_copy:xr.DataArray = xa.copy()

    if mask is not None:
        xa_weighted_average = __weighted_average_with_mask(
            dim, mask, weights, xa, xa_copy
        )
    elif weights is not None:
        xa_weighted_average = __weighted_average(
            dim, weights, xa, xa_copy
        )
    else:
        xa_weighted_average =  xa.mean(dim)

    return xa_weighted_average



    # %% [markdown]


def __weighted_average(dim, weights, xa, xa_copy):
    '''helper function for unmasked_average'''
    _, weights_all_dims = xr.broadcast(xa, weights)  # broadcast to all dims
    x_times_w = xa_copy * weights_all_dims
    xw_sum = x_times_w.sum(dim)
    x_tot = weights_all_dims.where(xa_copy.notnull()).sum(dim=dim)
    xa_weighted_average = xw_sum / x_tot
    return xa_weighted_average


def __weighted_average_with_mask(dim, mask, weights, xa, xa_copy):
    '''helper function for masked_average'''
    _, mask_all_dims = xr.broadcast(xa, mask)  # broadcast to all dims
    xa_copy = xa_copy.where(np.logical_not(mask))
    if weights is not None:
        _, weights_all_dims = xr.broadcast(xa, weights)  # broadcast to all dims
        weights_all_dims = weights_all_dims.where(~mask_all_dims)
        x_times_w = xa_copy * weights_all_dims
        xw_sum = x_times_w.sum(dim=dim)
        x_tot = weights_all_dims.where(xa_copy.notnull()).sum(dim=dim)
        xa_weighted_average = xw_sum / x_tot
    else:
        xa_weighted_average = xa_copy.mean(dim)
    return xa_weighted_average

def polarCentral_set_latlim(lat_lims, ax):
    ax.set_extent([-180, 180, lat_lims[0], lat_lims[1]], ccrs.PlateCarree())
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)

def add_map_features(ax):
    '''
    Single line command for xarray plots
    '''
    ax.coastlines()
    gl = ax.gridlines()
    ax.add_feature(cy.feature.BORDERS);
    gl = ax.gridlines()#draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    
    
def interpretNS(stringin):
    '''
    Interprets the string name of CALIOP longitudinally 
    averaged SLF by isotherm data and returns a weight 
    for the latitude range
    '''
    lat_range = stringin[11:-4]
    first = lat_range[:3]
    second = lat_range[-3:]
    
    if first[-1] is 'N':
        low = int(first[:2])
    else: low = -1 * int(first[:2])
    if second[-1] is 'N':
        high = int(second[:2])
    else: high = -1 * int(second[:2])
    avg_lat = np.abs(np.mean([low,high]))
    weight = np.cos(np.pi/180*avg_lat)
    
    return weight, min([low, high]), max([low,high])

def sp_map(*nrs, projection = ccrs.PlateCarree(), **kwargs):
    return plt.subplots(*nrs, subplot_kw={'projection':projection}, **kwargs)


def select_loc_to_pandas(dataset, coords):
    '''
    This function takes an xarray dataset and a (lat, lon) coordinate iterable.
    It selects the data for the location and returns a pandas datafram object.
    '''
    _xr_ds = xr.Dataset() # empty xarray Dataset
    for vals in dataset:
        _da = dataset[vals]        
        _da = _da.sel(lat=coords[0], lon=coords[1], method='nearest')
        _xr_ds[vals]=_da
    _df = _xr_ds.to_dataframe()
    return _df


def process_caliop(files, obs_dir):
    all_caliop = pd.DataFrame()
    weights = 0; avg = np.zeros(5);
    for file in files:
        _path = obs_dir + file # Get full file path
        _name = 'CALIOP_' + file[11:-4]   # Pick out latitude data from name
        _weight, _, _ = interpretNS(file)
        _slice = pd.read_table(_path, sep="\s+", names=['Isotherm', _name])
        all_caliop = pd.concat([all_caliop, _slice[_name]], axis=1, sort=False)

        # Do the averaging
        avg += _weight * _slice[_name]
        weights += _weight

    # Add the Isotherm colum and set it as an index
    all_caliop = pd.concat([all_caliop, _slice['Isotherm']], axis=1, sort=False)
    all_caliop = all_caliop.set_index('Isotherm')
    all_caliop['CALIOP Average'] = np.array(avg / weights)
    
    return all_caliop

def plot_slf_isotherms(ds):
    slf_isotm = ds['SLF_ISOTM_AVG']    
    
    fig1, axes1 = plt.subplots(nrows=3,ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=[20,10]);

    cmin_p = np.nanmin(slf_isotm)
    cmax_p = np.nanmax(slf_isotm)

    cmap_p = 'bwr'
    nlevels = 41
    cmap2 = plt.get_cmap(cmap_p)

    if cmin_p == cmax_p:
       cmax_p = cmax_p + 0.00001

    levels = np.linspace(cmin_p,cmax_p,nlevels)

    for data, ax in zip(slf_isotm, axes1.flatten()):
        iso = data['isotherms_mpc'].values - 273.15
        map = data.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='bwr', 
                        robust=True, add_colorbar = False, levels=levels)

        ax.set_title('SLF at %s' % str(iso), fontsize=18)
        ax.coastlines()

    cb_ax = fig1.add_axes([0.325, 0.05, 0.4, 0.04])
    cbar = plt.colorbar(map, cax=cb_ax, extend='both', orientation='horizontal', fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.set_xlabel('Supercooled Liquid Fraction', fontsize=16)

    fig1.suptitle('SLF distribution across isotherms', fontsize=28)

    return fig1

def add_weights(ds):
    '''
    And variable to ds for weighting of lat,lon variables
    '''
    gw = ds['gw']    

    _wgs = ds['TS'].copy().mean(dim = 'time', skipna=True)
    _wgs = (_wgs * 0 + 1) * gw # copy gw into the 2d array
    _wgs = _wgs / np.sum(_wgs)  # Normalize
    _wgs.name = 'cell_weight'

    ds['cell_weight'] = _wgs
    
    return ds

def process_for_slf(in_path, out_vars):
    '''
    Add SLF-relevant variables to netcdf file
    return a xr dataset with just variables of interest
    '''
    
    ds = xr.open_dataset(in_path + '.nc')
    ds = add_weights(ds)

    # Create new variable by dividing out the cloud fraction near each isotherm
    ds['SLF_ISOTM'] = ds['SLFXCLD_ISOTM'] / ds['CLD_ISOTM']

    # Select dates after a 3 month wind-up and average slf
    ds['SLF_ISOTM_AVG'] = ds['SLF_ISOTM'].sel(time=slice('0001-04-01', '0002-03-01')).mean(dim = 'time', skipna=True)

    stdev = np.std(ds['SLF_ISOTM_AVG'], axis=2)    
    
    ds_out = ds[out_vars]
    ds.close()
    
    return ds_out

def noresm_slf_to_df(ds, slf_files):
    '''
    Applies appropriate latitude masks to NorESM SLF based on CALIOP file names
    '''
    df = pd.DataFrame()

    df['Isotherm'] = ds['isotherms_mpc'].values - 273.15
    df['NorESM_Average'] = 100*masked_average(ds['SLF_ISOTM_AVG'], dim=['lat','lon'], weights=ds['cell_weight'])
    
    df['NorESM_Average_STD'] = 100*np.std(ds['SLF_ISOTM_AVG'], axis=(1,2))

    # Add each latitude range from NorESM, and the models stdev range
    for i in slf_files:
        _, _lowlat, _highlat = interpretNS(i)
        _mask = np.bitwise_or(ds['lat']<_lowlat, ds['lat']>_highlat)
        
        zone_mean = masked_average(ds['SLF_ISOTM_AVG'], dim=['lat','lon'], weights=ds['cell_weight'], mask=_mask)
        df['NorESM' + i[10:-4]] = 100*zone_mean

        # Add Standard Deviation
        df['NorESM' + i[10:-4] + '_STD'] = 100*np.std(ds['SLF_ISOTM_AVG'].sel(lat=slice(_lowlat,_highlat)), axis=(1,2)) 
        
        
    df = df.set_index('Isotherm')
    
    return df

def regress_1d(xdata, ydata):
    '''
    Returns an sklearn regression object trained on the passed data.
    Might be generalizable to higher dimensions.
    '''
    x = np.array(xdata).reshape(-1,1)
    y = np.array(ydata).reshape(-1,1)
    
    regressor = LinearRegression().fit(x, y)
    
    return regressor