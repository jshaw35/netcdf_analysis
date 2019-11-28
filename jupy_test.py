#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from imports import (
        pd, np, xr, mpl, plt, sns, os, 
        datetime, sys, crt, gridspec,
        polyfit, ccrs
        )

from functions import (
    masked_average, interpretNS, plot_slf_isotherms, 
    add_weights, process_caliop, process_for_slf,
    noresm_slf_to_df
    )

def main(model = None):
    if model == None:
        specific_model = '20191122_161009_sample_param_set_wbf1_inp1.nc'
        
    else:
        specific_model = model

    case = specific_model[:-3]
    print('Processing %s' % case)

    save_dir = '~/DATAOUT/'
    save_to = os.path.expanduser(save_dir)

    obs_dir = '/home/jonahks/drivemount/caliop_slfs/'
    output_dir = '/home/jonahks/drivemount/figures/'
    model_dir = '/home/jonahks/drivemount/mnth15runs/%s/' % case

    # In[3]:

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
    'cell_weight','SLF_ISOTM','SLF_ISOTM_AVG', 'TS'
    ]

    # In[ ]:

    #os. unmount
    #os. mountnird 'p/jonahks/'

    # ### Process CALIOP data
    # 
    # First, check if the output file already already exists.

    # In[5]:


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


    # This currently doesn't really work

    # In[6]:


    data_path = model_dir + case + 'slf_processed.nc'
    noresm_processed = os.path.exists(data_path)

    # Work around
    noresm_processed = False

    if noresm_processed:
    print('Grabbing data from %s' % data_path)
    ds = xr.open_dataset(data_path)

    else:
    print('Loading data from %s' % (model_dir + specific_model))
    ds = process_for_slf(model_dir + specific_model, relevant_vars)
    #    ds.to_netcdf(data_path)


    # In[7]:


    data_string = model_dir + 'MPC_COMPARE_CALIOP_NorESM.csv'
    all_slf_processed = os.path.exists(data_string)

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


    # In[8]:


    all_slf_clean = all_slf.dropna()
    all_slf


    # In[13]:


    fig1 = plt.figure(figsize=(10,6))#constrained_layout=True)
    spec1 = gridspec.GridSpec(ncols=3, nrows=1, figure=fig1)#, hspace=0.4)
    f1_ax1 = fig1.add_subplot(spec1[0, :-1])
    f1_ax2 = fig1.add_subplot(spec1[0, -1], sharey=f1_ax1)
    axes = [f1_ax1, f1_ax2]
    plt.setp(f1_ax2.get_yticklabels(), visible=False)

    fig1.gca().invert_yaxis()
    f1_ax1.invert_yaxis()
    f1_ax1.plot(all_slf['NorESM_90N-70N'], all_slf.index, label='Arctic >70N')
    f1_ax1.plot(all_slf['NorESM_Average'], all_slf.index, label='Global Average')
    f1_ax1.scatter(all_slf['CALIOP_90N-70N'], all_slf.index, label='CALIOP 70N-90N', c='b')
    f1_ax1.scatter(all_slf['CALIOP Average'], all_slf.index, label='CALIOP Global Avg.', c='orange')

    f1_ax1.legend()
    f1_ax1.set_title('Supercooled Liquid Fraction Comparison'); f1_ax1.set_ylabel('Isotherm (C)'); f1_ax1.set_xlabel('SLF (%)')

    f1_ax2.invert_yaxis()
    f1_ax2.scatter(all_slf['Arctic Error'], all_slf.index, c='b')
    f1_ax2.scatter(all_slf['Global Error'], all_slf.index, c='orange')
    f1_ax2.set_xlim((-30,10))
    f1_ax2.set_title('NorESM error'); f1_ax2.set_xlabel('SLF Error (%)')

    x1 = all_slf_clean['Arctic Error']; y1 = all_slf_clean.index;
    coef = np.polyfit(x1,y1,1)
    poly1d_fn = np.poly1d(coef) 
    f1_ax2.plot(x1, poly1d_fn(x1), color='b')

    x2 = all_slf_clean['Global Error']; y2 = all_slf_clean.index;
    coef = np.polyfit(x2,y2,1)
    poly1d_fn = np.poly1d(coef) 
    f1_ax2.plot(x2, poly1d_fn(x2), color='orange')


    # In[16]:


    filename = '%s_slf_comparison.png' % case
    filename
    if not os.path.exists(filename):
    fig1.savefig(model_dir + filename,format = 'png', dpi = 200)
    fig1.clf()


    # Plot the global SLF at each isotherm

    # In[10]:


    iso_fig = plot_slf_isotherms(ds)

    filename = '%s_noresm_slf_isotherms.png' % case
    filename
    if not os.path.exists(filename):
    iso_fig.savefig(model_dir + filename,format = 'png', dpi = 200)
    iso_fig.clf()

if __name__ == "__main__":
    main(sys.argv[1])
