# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Compare observations of cloudtop SLF from Olimpia to those produced by NorESM
# Using a modified micro_mg_cam.F90 that can be found in noresm_slf repo.

# ## Necessary Imports

# +
import sys
# Add common resources folder to path
sys.path.append("/mnt/mcc-ns9600k/jonahks/git_repos/netcdf_analysis/Common/")

from imports import (
    pd, np, xr, mpl, plt, sns, os, 
    datetime, sys, crt, gridspec,
    polyfit, ccrs, LinearRegression, metrics
    )

from functions import (
    masked_average, interpretNS, plot_slf_isotherms, 
    add_weights, process_caliop, process_for_slf,
    noresm_slf_to_df, regress_1d
    )

# %matplotlib inline
# -

# Set up directories based on where the program is being run from

# +
host = os.uname()[1]
if 'jupyter' in host.split('-'): # Check if running on NIRD through the Jupyter Hub
    print('Running through MC2 Jupyter Hub')
    op_dir = '/mnt/mcc-ns9600k/jonahks/'
    os.chdir(op_dir)

else:  # Assume that we're running on a local machine and mounting NIRD
    print('Running on %s, attempting to mount ns9600k/jonahks/ from NIRD' % str(host))
    os.system('fusermount -zu ~/drivemount/')  # unmount first
    os.system('sshfs jonahks@login.nird.sigma2.no:"p/jonahks/" ~/drivemount/')    # Calling mountnird from .bashrc doesn't work
    os.chdir('/home/jonahks/drivemount/')
    save_dir = '~/DATAOUT/'
    save_to = os.path.expanduser(save_dir)

obs_dir = 'caliop_slfs/'
output_dir = 'figures/'
case_dir = 'mnth15runs/' # inconsistent label compared to jupy_test
    
# Check that each important directory can be accessed:    
access_paths = os.path.exists(obs_dir) and os.path.exists(output_dir) and os.path.exists(case_dir)
print('Can access all directory paths:', access_paths)
# -

# List case to choose one of interest

all_cases = os.listdir('mnth15runs/')
all_cases

inp_cases = os.listdir('inp_validation/')
inp_cases

# Pick run to analyze

# +
case = '20200207_145043_singleparam_icenucmod_wbf_1_inp_10'

#run_dir = 'mnth15runs/%s/' % case # inconsistent label compared to jupy_test
run_dir = 'inp_validation/%s' % case
print(run_dir, os.path.exists(run_dir))
# -

# Ongoing list of variables of interest

relevant_vars = [
     'CLDFREE', 'CLDHGH','CLDICE', 'CLDLIQ', 'CLDLOW','CLDMED',
     'CLDTAU','CLDTOT','CLD_ISOTM','CLD_ISOTM_NONSIM','CLD_SLF',
     'CLD_SLF_NONSIM','CLOUD','CLOUDCOVER_CLUBB','CLOUDFRAC_CLUBB',
     'CONCLD', 'BERGO','BERGOXCLD_ISOTM','BERGOXCLD_ISOTM_NONSIM',
     'MG_SADICE','MG_SADLIQ','MNUCCCO','MNUCCDO','MNUCCDOhet',
     'MNUCCRO','MNUCCTO','NUMICE','NUMLIQ','NUMRAI','NUMSNO',
     'N_AER','PRECIPBINOCC_CC','PRECIPBINOCC_CL','PRECIPBINOCC_CT',
     'PRECIPBINRATE_CC','PRECIPBINRATE_CL','PRECIPBINRATE_CT', 
     'SADICEXCLD_ISOTM','SADICEXCLD_ISOTM_NONSIM','SADLIQXCLD_ISOTM',
     'SADLIQXCLD_ISOTM_NONSIM','SLFXCLD_ISOTM','SLFXCLD_ISOTM_NONSIM',
     'cell_weight', 'TS', 'CT_CLD_ISOTM',
     'CT_SLF', 'CT_SLFXCLD_ISOTM', 
     'AREI','FREQI','NUMICE','NUMICE10s','DSTFREZIMM', 'DSTFREZCNT',
     'DSTFREZDEP','NNUCCTO', 'NNUCCCO', 'NNUDEPO', 'NIHOMOO',
     'HOMOO','NIMIX_CNT','NIMIX_IMM','RELHUM','T','RHO_CLUBB'
    ]

# This currently doesn't really work. WATT

# Load cloudtop SLFs and construct global and regional averages and stds
# This should be more flexible for arbitrary latitude ranges.

# +
# Load NorESM data
try:
    _ds = xr.open_dataset('%s/%s.nc' % (run_dir,case))
except:
    _ds = xr.open_dataset('%s/atm/hist/%s.cam.h0.2000-01.nc' % (run_dir,case))
if (len(_ds['time']) > 1):
    try:
        ds = _ds.sel(time=slice('0001-04-01', '0002-03-01'))
    except:
        ds = _ds.sel(time=slice('2000-04-01', '2001-03-01'))
else:
    ds = _ds
ds = add_weights(ds) # still has time here

ds['CT_SLF'] = ds['CT_SLFXCLD_ISOTM']/ds['CT_CLD_ISOTM']
ct_slf_noresm = ds['CT_SLF']

ds['CT_SLF_ISOTM_AVG'] = ds['CT_SLF'].mean(dim = 'time', skipna=True)

# Load CALIOP data
ct_slf_caliop = xr.open_dataset('caliop_cloudtop/cloudtop_slfs.nc')
# -

# Define latitude ranges of interest.

slfvars = ['cell_weight', 'gw', 'TS', 'CT_SLF','CT_SLF_ISOTM_AVG','CT_SLFXCLD_ISOTM',
           'CT_CLD_ISOTM','SLFXCLD_ISOTM','CLD_ISOTM',
           'AREI','FREQI','NUMICE','NUMICE10s','DSTFREZIMM', 'DSTFREZCNT',
           'DSTFREZDEP','NNUCCTO', 'NNUCCCO', 'NNUDEPO', 'NIHOMOO',
           'HOMOO','NIMIX_CNT','NIMIX_IMM','RELHUM','T','RHO_CLUBB']
doop = ds[slfvars]
del doop.attrs['_NCProperties'] # fixes a bug for some reasons: https://github.com/pydata/xarray/issues/2822
doop.to_netcdf(path='%s/%s_slfvars.nc' % (run_dir,case))

bands = {'Global':[-90,90],'Arctic':[66.667,90],'Antarctic':[-90,-66.667],'CALIOP Arctic':[66.667,82]}
df = pd.DataFrame()

#df['isotherm'] = slf1['isotherm']
#df = df.set_index('isotherm')
for i in bands:
    _rng = bands[i]
    
    # Different resolutions, needs different masks. Remember weird mask sign convention: true=not included
    mask1 = np.bitwise_or(ct_slf_caliop['lat']<_rng[0], ct_slf_caliop['lat']>_rng[1])
    mask2 = np.bitwise_or(ct_slf_noresm['lat']<_rng[0], ct_slf_noresm['lat']>_rng[1])
    
    weight1 = ct_slf_caliop['cell_weight']
    weight2 = ds['cell_weight'] #*ds['CT_CLD_ISOTM'] #Not sure about this
    
    slf1 = 100*masked_average(ct_slf_caliop['SLF'], dim=['lat','lon'],weights=weight1, mask=mask1)
    slf2 = 100*masked_average(ct_slf_noresm, dim=['lat','lon','time'],weights=weight2, mask=mask2)
    
    stdev1 = 100*np.std(ct_slf_caliop['SLF'].sel(lat=slice(_rng[0],_rng[1])), axis=(0,1))
    stdev2 = 100*np.std(ct_slf_noresm.sel(lat=slice(_rng[0],_rng[1])), axis=(0,2,3))
    
    df['CALIOP %s SLF' % i] = slf1
    df['CALIOP %s StDev' % i] = stdev1
    df['NorESM %s SLF' % i] = slf2
    df['NorESM %s StDev' % i] = stdev2
df['isotherm'] = slf1['isotherm']
df

data_string = '%s%scloudtop_slf_comparison.csv' % (run_dir, case)
df.to_csv(path_or_buf = data_string)

# +
fig1 = plt.figure(figsize=(10,6))#constrained_layout=True)
spec1 = gridspec.GridSpec(ncols=3, nrows=1, figure=fig1)#, hspace=0.4)
f1_ax1 = fig1.add_subplot(spec1[0, :-1])
f1_ax2 = fig1.add_subplot(spec1[0, -1], sharey=f1_ax1)
axes = [f1_ax1, f1_ax2]
plt.setp(f1_ax2.get_yticklabels(), visible=False)

#isos = np.array(all_slf_clean.index).reshape(-1,1)

fig1.gca().invert_yaxis()
f1_ax1.set_title('Supercooled Liquid Fraction Comparison'); f1_ax1.set_ylabel('Isotherm (C)'); f1_ax1.set_xlabel('SLF (%)')
f1_ax2.set_title('NorESM error'); f1_ax2.set_xlabel('SLF Error (%)')

colors = ['blue', 'orange', 'red', 'purple']

for b,c in zip(bands, colors):
    if (b == 'Arctic' or b == 'CALIOP Arctic'):
        f1_ax1.errorbar(df['CALIOP %s SLF' % b], df['isotherm'], xerr=2*df['CALIOP %s StDev' % b], label='CALIOP %s' % b, color = c, fmt='o', marker='D')
        f1_ax1.plot(df['NorESM %s SLF' % b], df['isotherm'], color = c) #, label=b)
        f1_ax1.fill_betweenx(df['isotherm'], df['NorESM %s SLF' % b] - 2*df['NorESM %s StDev' % b], df['NorESM %s SLF' % b] + 2*df['NorESM %s StDev' % b], alpha=0.2, color=c)

        slf_error = df['NorESM %s SLF' % b] - df['CALIOP %s SLF' % b]
        f1_ax2.scatter(slf_error, df['isotherm'], color=c, label=b)

#_r = regress_1d(isos, all_slf_clean[error])
#_s = _r.score(isos, all_slf_clean[error])
#f1_ax2.plot(_r.predict(isos), isos, color=color, label = ('$R^2 = %f$' % _s))

f1_ax1.set_xlim((0,105))
f1_ax1.legend()
f1_ax2.legend()

fig1.suptitle(case, fontsize=16)

# +
fig1 = plt.figure(figsize=(10,6))#constrained_layout=True)
spec1 = gridspec.GridSpec(ncols=3, nrows=1, figure=fig1)#, hspace=0.4)
f1_ax1 = fig1.add_subplot(spec1[0, :-1])
f1_ax2 = fig1.add_subplot(spec1[0, -1], sharey=f1_ax1)
axes = [f1_ax1, f1_ax2]
plt.setp(f1_ax2.get_yticklabels(), visible=False)

#isos = np.array(all_slf_clean.index).reshape(-1,1)

fig1.gca().invert_yaxis()
f1_ax1.set_title('Supercooled Liquid Fraction Comparison'); f1_ax1.set_ylabel('Isotherm (C)'); f1_ax1.set_xlabel('SLF (%)')
f1_ax2.set_title('NorESM error'); f1_ax2.set_xlabel('SLF Error (%)')

colors = ['blue', 'orange', 'red','purple']
#print(bands)
for b,c in zip(bands, colors):
    if b == 'CALIOP Arctic':
        f1_ax1.errorbar(df['CALIOP %s SLF' % b], df['isotherm'], xerr=2*df['CALIOP %s StDev' % b], label='CALIOP %s' % b, color = c, fmt='o', marker='D')
        f1_ax1.plot(df['NorESM %s SLF' % b], df['isotherm'], color = c) #, label=b)
        f1_ax1.fill_betweenx(df['isotherm'], df['NorESM %s SLF' % b] - 2*df['NorESM %s StDev' % b], df['NorESM %s SLF' % b] + 2*df['NorESM %s StDev' % b], alpha=0.2, color=c)

        slf_error = df['NorESM %s SLF' % b] - df['CALIOP %s SLF' % b]
        f1_ax2.scatter(slf_error, df['isotherm'], color=c, label=b)

#_r = regress_1d(isos, all_slf_clean[error])
#_s = _r.score(isos, all_slf_clean[error])
#f1_ax2.plot(_r.predict(isos), isos, color=color, label = ('$R^2 = %f$' % _s))

f1_ax1.set_xlim((0,105))
f1_ax1.legend()
f1_ax2.legend()

fig1.suptitle(case, fontsize=16)
# -

filename = '%s_slf_comparison.png' % case
filename
if not os.path.exists(filename):
    fig1.savefig(run_dir + filename,format = 'png', dpi = 200)
    fig1.clf()

# Plot the global SLF at each isotherm

iso_fig = plot_slf_isotherms(ds)

filename = '%s_noresm_slf_isotherms.png' % case
filename
if not os.path.exists(filename):
    iso_fig.savefig(run_dir + filename,format = 'png', dpi = 200)
    iso_fig.clf()

# # TRash!
