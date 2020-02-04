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
    model_dir = '/mnt/mcc-ns9600k/jonahks/'
    os.chdir(model_dir)

else:  # Assume that we're running on a local machine and mounting NIRD
    print('Running on %s, attempting to mount ns9600k/jonahks/ from NIRD' % str(host))
    os.system('fusermount -zu ~/drivemount/')  # unmount first
    os.system('sshfs jonahks@login.nird.sigma2.no:"p/jonahks/" ~/drivemount/')    # Calling mountnird from .bashrc doesn't work
    os.chdir('/home/jonahks/drivemount/')
    save_dir = '~/DATAOUT/'
    save_to = os.path.expanduser(save_dir)

obs_dir = 'caliop_slfs/'
output_dir = 'figures/'
model_dir = 'mnth15runs/' # inconsistent label compared to jupy_test
    
# Check that each important directory can be accessed:    
access_paths = os.path.exists(obs_dir) and os.path.exists(output_dir) and os.path.exists(model_dir)
print('Can access all directory paths:', access_paths)
# -

all_cases = os.listdir('mnth15runs/')
all_cases

# +
#specific_model = '20191219_151155_singleparam_cttest_wbf_1_inp_1.cam.h0.0001-01.nc'
specific_model = '20200112_002538_singleparam_nudge_wbf_1_inp_0.nc'
#ct_val = '20191219_151155_singleparam_cttest_wbf_1_inp_1.cam.h0.0001-01.nc'
case = specific_model[:-3]

model_dir = 'mnth15runs/%s/' % case # inconsistent label compared to jupy_test
#model_dir = 'NorESM_validation/%s' % ct_val
os.path.exists(model_dir)
# -

relevant_vars = [
     'CLDFREE', 'CLDHGH','CLDICE', 'CLDLIQ', 'CLDLOW','CLDMED',
     'CLDTAU','CLDTOT','CLD_ISOTM','CLD_ISOTM_NONSIM','CLD_SLF',
     'CLD_SLF_NONSIM','CLOUD','CLOUDCOVER_CLUBB','CLOUDFRAC_CLUBB',
     'CONCLD', 'BERGO','BERGOXCLD_ISOTM','BERGOXCLD_ISOTM_NONSIM',
     'BERGSO','BERGSOXCLD_ISOTM','BERGSOXCLD_ISOTM_NONSIM',
     'MG_SADICE','MG_SADLIQ','MNUCCCO','MNUCCDO','MNUCCDOhet',
     'MNUCCRO','MNUCCTO','NUMICE','NUMLIQ','NUMRAI','NUMSNO',
     'N_AER','PRECIPBINOCC_CC','PRECIPBINOCC_CL','PRECIPBINOCC_CT',
     'PRECIPBINRATE_CC','PRECIPBINRATE_CL','PRECIPBINRATE_CT', 
     'SADICEXCLD_ISOTM','SADICEXCLD_ISOTM_NONSIM','SADLIQXCLD_ISOTM',
     'SADLIQXCLD_ISOTM_NONSIM','SLFXCLD_ISOTM','SLFXCLD_ISOTM_NONSIM',
     'cell_weight','SLF_ISOTM','SLF_ISOTM_AVG', 'TS', 'AREI'
    ]

# ### Process CALIOP data
#
# First, check if the output file already already exists.

# +
data_string = obs_dir + 'MPC_ISO_CALIOP_NorESM.csv'
caliop_processed = os.path.exists(data_string)

#Pick out the right files
file_str = '.dat'
obs_files = os.listdir(obs_dir) # All files in directory
slf_files = [x for x in obs_files if file_str in x]   # files with the CALIOP string
slf_files.sort()

if caliop_processed:
    print('Grabbing data from %s' % data_string)
    all_caliop = pd.read_csv(data_string)

else:
    print('Writing data to %s' % data_string)
    all_caliop = process_caliop(slf_files, obs_dir)
    
    all_caliop.to_csv(path_or_buf = data_string)

all_caliop = all_caliop.set_index('Isotherm')
all_caliop
# -

# This currently doesn't really work

model_dir

# +
data_path = model_dir + case + '_slf_processed.nc'
noresm_processed = os.path.exists(data_path)

# Work around
noresm_processed = False

if noresm_processed:
    print('Grabbing data from %s' % data_path)
    ds = xr.open_dataset(data_path)

else:
#    print('Processing data from %s' % (model_dir + specific_model))
    print('Processing data from %s' % (model_dir + case))
#    ds = process_for_slf(model_dir + specific_model, relevant_vars)
    ds = process_for_slf(model_dir + case, relevant_vars)
#    ds.to_netcdf(data_path)

# +
data_string = '%s%s_slf_caliop_comparison.csv' % (model_dir, case)

all_slf_processed = os.path.exists(data_string)
all_slf_processed=False # just to reset and add stdev
if all_slf_processed:
    print('Loading from %s' % data_string)
    all_slf = pd.read_csv(data_string)
    
else:
    print('Processing from netcdf-based xarray object and saving to %s' % data_string)
    df = noresm_slf_to_df(ds, slf_files)
    
    all_slf = pd.concat([all_caliop, df], axis=1, sort=False)
    all_slf['Arctic Error'] = all_slf['NorESM_90N-70N'] - all_slf['CALIOP_90N-70N']
    all_slf['Global Error'] = all_slf['NorESM_Average'] - all_slf['CALIOP Average']

    all_slf.to_csv(path_or_buf = data_string)

# Catch an annoying error
try:
    all_slf = all_slf.set_index('Isotherm')
except: pass    

all_slf_clean = all_slf.dropna()
# -

all_slf

all_slf_clean

# +
fig1 = plt.figure(figsize=(10,6))#constrained_layout=True)
spec1 = gridspec.GridSpec(ncols=3, nrows=1, figure=fig1)#, hspace=0.4)
f1_ax1 = fig1.add_subplot(spec1[0, :-1])
f1_ax2 = fig1.add_subplot(spec1[0, -1], sharey=f1_ax1)
axes = [f1_ax1, f1_ax2]
plt.setp(f1_ax2.get_yticklabels(), visible=False)

isos = np.array(all_slf_clean.index).reshape(-1,1)

fig1.gca().invert_yaxis()
f1_ax1.set_title('Supercooled Liquid Fraction Comparison'); f1_ax1.set_ylabel('Isotherm (C)'); f1_ax1.set_xlabel('SLF (%)')
f1_ax2.set_title('NorESM error'); f1_ax2.set_xlabel('SLF Error (%)')

colors = ['blue', 'orange']
models = ['NorESM_90N-70N', 'NorESM_Average']
caliops = ['CALIOP_90N-70N', 'CALIOP Average']
errors = ['Arctic Error', 'Global Error']

# Wrap everything up and iterate!
for color, model, caliop, error in zip(colors, models, caliops, errors):
    f1_ax1.plot(all_slf[model], all_slf.index, label=model, color = color) # caliop values
    f1_ax1.fill_betweenx(all_slf.index, all_slf[model] - 2*all_slf[model + "_STD"], all_slf[model] + 2*all_slf[model + "_STD"], alpha=0.2)
    
    f1_ax1.scatter(all_slf[caliop], all_slf.index, label=caliop, color = color)
    f1_ax2.scatter(all_slf[error], all_slf.index, color=color)
    _r = regress_1d(isos, all_slf_clean[error])
    _s = _r.score(isos, all_slf_clean[error])
    f1_ax2.plot(_r.predict(isos), isos, color=color, label = ('$R^2 = %f$' % _s))
    
#f1_ax2.set_xlim((-30,10))
f1_ax1.legend()
f1_ax2.legend()

fig1.suptitle(case, fontsize=16)
# -

filename = '%s_slf_comparison.png' % case
filename
if not os.path.exists(filename):
    fig1.savefig(model_dir + filename,format = 'png', dpi = 200)
    fig1.clf()

# Plot the global SLF at each isotherm

iso_fig = plot_slf_isotherms(ds)

filename = '%s_noresm_slf_isotherms.png' % case
filename
if not os.path.exists(filename):
    iso_fig.savefig(model_dir + filename,format = 'png', dpi = 200)
    iso_fig.clf()


