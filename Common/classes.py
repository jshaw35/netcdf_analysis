from imports import *
from functions import *

class SLF_Metric(object):
    
    def __init__(self, casedir):
        for cls in reversed(self.__class__.mro()):
            if hasattr(cls, 'init'):
                cls.init(self, casedir)

    def init(self, casedir):
        # Catch directory issues from the get-go
        if not os.path.exists(casedir):
            print('The case directory %s does not exist' % casedir)
#            sys.exit() # Is this a problem, does it need to be here?
            
        self.case_dir = casedir
        self.origin = None
        self.cases = {}

    def set_origin(self, case): # case is a string of the case name
        if case not in self.cases: # check if it is already in the dictionary
            self.add_case(case) 
        self.origin = self.cases[case] # should have a check for the object type
        
    def add_case(self, case): # Dictionary can be overwritten, should be fine
        
        self.cases[case] = SLF_case(self.case_dir, case) # Add an slf_case object to the cases dictionary
    
    def get_case(self, case):
        if case not in self.cases: # check if it is already in the dictionary
            print('Could not find %s in the SLF_Metric object' % case)
            return None
        else:
            return self.cases[case] # should have a check for the object type
    
    def get_cases(self):
        return self.cases
#        return list(self.cases.keys()) # This is a dict_keys object and should be an iterable
    
    def get_origin(self):
        return self.origin
    
    def plot_vectors(self):
        vector_plot = plt.figure(figsize=[10,10])

        origin = [self.origin.error_slope, self.origin.error_20c] # origin point (x,y)    
        for i in self.cases:
            _case = self.cases[i]
            # xy is the endpoint, xytext is the origin
            label = 'WBF: %.3f INP: %.3f, %.3f SLF RMS error' % (_case.wbf_mult, _case.inp_mult, _case.rms_error)
            plt.scatter(_case.error_slope, _case.error_20c, label=label)
            plt.annotate('', xy=(_case.error_slope, _case.error_20c), xytext=(origin[0], origin[1]), 
                        arrowprops=dict(facecolor='black', shrink=0., width=0.5, lw=1.0),
                        label = label)

        plt.scatter([[0]], [[0]], color = 'black')
        plt.xlabel('Error Slope (SLF % / C)')
        plt.ylabel('0C Isotherm SLF Error')
        plt.grid()
        plt.legend()

        return vector_plot
    
    def plot_isos(self):
        return None

    def plot_isos_bulk(self):
        isos_plot = plt.figure(figsize=[10,10])
        plt.gca().invert_yaxis()
        
        for i, j in enumerate(self.cases):
            _case = self.cases[j]
            if i==0:
                plt.plot(_case.case_ds['CALIOP_70S-90S'], _case.case_ds.index, color='black', label='CALIOP')
                

#            plt.scatter(slf, slf['isotherms_mpc'] - 273, label=_case.case)
            try:
                plt.errorbar(_case.case_ds['NorESM_90N-70N'], _case.case_ds.index, xerr=_case.case_ds['NorESM_90N-70N_STD'], label=_case.case)
            except:
                plt.plot(_case.case_ds['NorESM_90N-70N'], _case.case_ds.index, label=_case.case)
        plt.xlabel('SLF Error (%)')
        plt.ylabel('Isotherm (C)')
        plt.legend()
        plt.title('Bulk cloud SLF trends with INP modifications', fontsize=16)
        
        return isos_plot
    
    def plot_isos_error(self):
        isos_error_plot = plt.figure(figsize=[10,10])
        plt.gca().invert_yaxis()
        isos_error_plot.axes[0].axvline(color='grey', linestyle='--')
        for j, i in enumerate(self.cases):
            _case = self.cases[i]
            label = 'WBF: %.3f INP: %.3f, %.3f SLF RMS error' % (_case.wbf_mult, _case.inp_mult, _case.rms_error)
            plt.scatter(_case.error, _case.isos, label=label)
            _r = regress_1d(_case.isos, _case.error)
            plt.plot(_r.predict(_case.isos), _case.isos)
            
        plt.xlabel('SLF Error (%)')
        plt.ylabel('Isotherm (C)')
        plt.legend()
                               
        return isos_error_plot
    
    def plot_parameterspace(self):
        parameterspace_plot = plt.figure(figsize=[10,10])       
#        origin = [test.get_origin().wbf_mult, test.get_origin().inp_mult] # origin point (x,y)

        plt.yscale('log')
        plt.xscale('log')

        for i in self.cases:
            _case = self.cases[i]
            label = 'WBF: %.3f INP: %.3f, %.3f SLF RMS error' % (_case.wbf_mult, _case.inp_mult, _case.rms_error)
            
#            if _case.inp_mult == 0: Not sure how to catch this error
            plt.annotate('', xy=(_case.wbf_mult, _case.inp_mult), xytext=(self.origin.wbf_mult, self.origin.inp_mult), 
                        arrowprops=dict(facecolor='black', shrink=0., width=0.5, lw=1.0),
                        )
            plt.scatter(_case.wbf_mult, _case.inp_mult, label=label)

        plt.scatter([self.origin.wbf_mult], [self.origin.inp_mult], label = 'Initial Parametrization')
        plt.scatter([[0]], [[0]], color = 'black')
        #plt.xlim((1e-3,1e2));
        #plt.ylim((1e-3,1e2))
        plt.xlabel('WBF Multiplier')
        plt.ylabel('INP Multiplier')
        plt.grid()
#        plt.legend()
        
        return parameterspace_plot
    
class SLF_case:
    
    def __init__(self, casedir, case):
        self.case_dir = casedir # General directory where all cases are stored
        self.case = case # The origin case name (wbf = 1, slf = 1)
        self.case_ds = pd.read_csv('%s/%s_slf_caliop_comparison.csv' % (casedir + case, case)).dropna()

        # Catch an annoying error
        try:
            self.case_ds = self.case_ds.set_index('Isotherm')
        except: pass
        self.case_ds = self.case_ds.drop(0.0) # remove the OC isotherm because it sucks, optional*
        
        _temp_parsed = self.case.split('_') # Parse string to get model parameters
        self.date, self.time, self.paramfile = _temp_parsed[:3]
        self.wbf_mult = np.float(_temp_parsed[-3]); self.inp_mult = np.float(_temp_parsed[-1])
        self.isos = np.array(self.case_ds.index).reshape(-1,1)
        self.error = self.case_ds['Arctic Error']

        # calculate the two parameters of interest, the error slope and the -20C isotherm error
        self.error_20c = np.float(self.error.loc[[-20]]) # SLF error at the -20C isotherm (metric 2)        
        self.error_slope = np.float(regress_1d(self.isos, self.error).coef_)
        self.rms_error = np.sqrt(np.nanmean(np.power(self.error,2)))
        
class CT_SLF_Metric(object):
    
    def __init__(self, casedir):
        for cls in reversed(self.__class__.mro()):
            if hasattr(cls, 'init'):
                cls.init(self, casedir)

    def init(self, casedir):
        # Catch directory issues from the get-go
        if not os.path.exists(casedir):
            print('The case directory %s does not exist' % casedir)
#            sys.exit() # Is this a problem, does it need to be here?
            
        self.case_dir = casedir
        self.origin = None
        self.cases = {}
        try:
            self.ct_caliop_slf = xr.open_dataset('caliop_cloudtop/cloudtop_slfs.nc')
        except:
            print('Could not load cloudtop cloud CALIOP slfs from caliop_cloudtop/cloudtop_slfs.nc')
        try:
            self.bulk_caliop_slf = pd.read_csv('caliop_slfs/MPC_ISO_CALIOP_NorESM.csv')
        except: 
            print('Could not load bulk cloud CALIOP slfs from caliop_slfs/MPC_ISO_CALIOP_NorESM.csv')
            

    def set_origin(self, case): # case is a string of the case name
        if case not in self.cases: # check if it is already in the dictionary
            self.add_case(case) 
        self.origin = self.cases[case] # should have a check for the object type
        
    def add_case(self, case): # Dictionary can be overwritten, should be fine
        self.cases[case] = CT_SLF_case(self.case_dir, case) # Add a ct_slf_case object to the cases dictionary
    
    def get_case(self, case):
        if case not in self.cases: # check if it is already in the dictionary
            print('Could not find %s in the SLF_Metric object' % case)
            return None
        else:
            return self.cases[case] # should have a check for the object type
    
    def get_cases(self):
        return self.cases
#        return list(self.cases.keys()) # This is a dict_keys object and should be an iterable
    
    def get_origin(self):
        return self.origin
    
    def plot_vectors(self):
        vector_plot = plt.figure(figsize=[10,10])

        origin = [self.origin.error_slope, self.origin.error_20c] # origin point (x,y)    
        for i in self.cases:
            _case = self.cases[i]
            # xy is the endpoint, xytext is the origin
            label = 'WBF: %.3f INP: %.3f, %.3f SLF RMS error' % (_case.wbf_mult, _case.inp_mult, _case.rms_error)
            plt.scatter(_case.error_slope, _case.error_20c, label=label)
            plt.annotate('', xy=(_case.error_slope, _case.error_20c), xytext=(origin[0], origin[1]), 
                        arrowprops=dict(facecolor='black', shrink=0., width=0.5, lw=1.0),
                        label = label)

        plt.scatter([[0]], [[0]], color = 'black')
        plt.xlabel('Error Slope (SLF % / C)')
        plt.ylabel('0C Isotherm SLF Error')
        plt.grid()
        plt.legend()

        return vector_plot
    
    def plot_isos_bulk(self):
        isos_plot = plt.figure(figsize=[10,10])
        plt.gca().invert_yaxis()

        plt.plot(self.bulk_caliop_slf['CALIOP_70S-90S'], self.bulk_caliop_slf['Isotherm'], 
                 color='black', label='CALIOP', linestyle='-', marker='D')
        
        for i in self.cases:
            _run = self.cases[i]
            _case = _run.case_da
            try:
                _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM']).sel(time=slice('0001-04-01',
                                    '0002-03-01'))
            except:
                _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM']).sel(time=slice('2000-04-01',
                                    '2001-03-01'))
                
            weight = _case['cell_weight']
            mask = np.bitwise_or(_case['lat']<70, _case['lat']>90)
            slf = 100*masked_average(_case['SLF_ISOTM'], dim=['lat','lon', 'time'], weights=weight, mask=mask)

            _case['SLF_ISOTM_TAVG'] = _case['SLF_ISOTM'].mean(dim = 'time', skipna=True)
#            weight = _case['cell_weight']#*_case['CT_CLD_ISOTM'] #Not sure about this
#            mask = np.bitwise_or(_case['lat']<70, _case['lat']>90)
            
            stdev = 100*np.std(_case['SLF_ISOTM_TAVG'].sel(lat=slice(70,90)), axis=(1,2))

#            slf = 100*masked_average(_case['SLF_ISOTM_TAVG'], dim=['lat','lon'], weights=weight, mask=mask)
#            plt.plot(slf, slf['isotherms_mpc'] - 273, label=_case.case, linestyle='-', marker='o')
            plt.scatter(slf, slf['isotherms_mpc'] - 273, label=_run.label)
#            plt.errorbar(slf, slf['isotherms_mpc'] - 273, xerr=stdev, label=_case.case)
        
        plt.xlabel('SLF Error (%)')
        plt.ylabel('Isotherm (C)')
        plt.legend()
        plt.title('Bulk cloud SLF trends with INP modifications', fontsize=16)
        
        return isos_plot
    
    def plot_isos_ct(self):
        isos_plot = plt.figure(figsize=[10,10])
        plt.gca().invert_yaxis()

        caliop_weight = self.ct_caliop_slf['cell_weight']
        caliop_mask = np.bitwise_or(self.ct_caliop_slf['lat']<66.667, self.ct_caliop_slf['lat']>90)
        caliop_slf = 100*masked_average(self.ct_caliop_slf['SLF'], dim=['lat','lon'], weights=caliop_weight, mask=caliop_mask)
        caliop_stdev = 100*np.std(self.ct_caliop_slf['SLF'].sel(lat=slice(66.667,90)), axis=(0,1))
        
        plt.errorbar(caliop_slf, caliop_slf['isotherm'], xerr=caliop_stdev, label='CALIOP SLF', fmt='o', marker='D', color = 'black')

        for j, i in enumerate(self.cases):
            _run = self.cases[i]
            _case = _run.case_da
            _case['CT_SLF_TAVG'] = _case['CT_SLF'].mean(dim = 'time', skipna=True)
#            weight = _case['cell_weight']#*_case['CT_CLD_ISOTM'] #Not sure about this weighting
#            mask = np.bitwise_or(_case['lat']<70, _case['lat']>90)
#            slf = 100*masked_average(_case['CT_SLF_TAVG'], dim=['lat','lon'], weights=weight, mask=mask)

            weight = _case['cell_weight']*_case['CT_CLD_ISOTM'] #Not sure about this weighting
            mask = np.bitwise_or(_case['lat']<70, _case['lat']>90)
            slf = 100*masked_average(_case['CT_SLF_TAVG'], dim=['lat','lon', 'time'], weights=weight, mask=mask)

            plt.scatter(slf, slf['isotherms_mpc'] - 273, label=_run.label)
            
        plt.xlabel('SLF Error (%)')
        plt.ylabel('Isotherm (C)')
        plt.legend()
        plt.title('Cloudtop SLF trends with INP modifications', fontsize=16)
                               
        return isos_plot
    
    def plot_parameterspace(self):
        parameterspace_plot = plt.figure(figsize=[10,10])       
#        origin = [test.get_origin().wbf_mult, test.get_origin().inp_mult] # origin point (x,y)

        plt.yscale('log')
        plt.xscale('log')

        for i in self.cases:
            _case = self.cases[i]
            label = 'WBF: %.3f INP: %.3f, %.3f SLF RMS error' % (_case.wbf_mult, _case.inp_mult, _case.rms_error)
            
#            if _case.inp_mult == 0: Not sure how to catch this error
            plt.annotate('', xy=(_case.wbf_mult, _case.inp_mult), xytext=(self.origin.wbf_mult, self.origin.inp_mult), 
                        arrowprops=dict(facecolor='black', shrink=0., width=0.5, lw=1.0),
                        )
            plt.scatter(_case.wbf_mult, _case.inp_mult, label=label)

        plt.scatter([self.origin.wbf_mult], [self.origin.inp_mult], label = 'Initial Parametrization')
        plt.scatter([[0]], [[0]], color = 'black')
        #plt.xlim((1e-3,1e2));
        #plt.ylim((1e-3,1e2))
        plt.xlabel('WBF Multiplier')
        plt.ylabel('INP Multiplier')
        plt.grid()
#        plt.legend()
        
        return parameterspace_plot
    
class CT_SLF_case:
    
    def __init__(self, casedir, case):
        self.case_dir = casedir # General directory where all cases are stored
        self.case = case # The origin case name (wbf = 1, slf = 1)
        
        self.add_ds()
        
        _temp_parsed = self.case.split('_') # Parse string to get model parameters
        self.date, self.time, self.paramfile = _temp_parsed[:3]
        self.wbf_mult = np.float(_temp_parsed[-3]); self.inp_mult = np.float(_temp_parsed[-1])
        self.label = 'WBF: %.3f INP: %.3f' % (self.wbf_mult, self.inp_mult)
        
    def add_ds(self):
        strg = '%s%s/%s_slfvars.nc' % (self.case_dir, self.case, self.case)
        print(strg)
        _ds = xr.open_dataset('%s%s/%s_slfvars.nc' % (self.case_dir, self.case, self.case))
        if (len(_ds['time']) > 1):
            try:
                ds = _ds.sel(time=slice('0001-04-01', '0002-03-01'))
            except:
                ds = _ds.sel(time=slice('2000-04-01', '2001-03-01'))
        #        ds = _ds.isel(time=slice(3,15))
        ds = add_weights(ds) # still has time here

        ds['CT_SLF'] = ds['CT_SLFXCLD_ISOTM']/ds['CT_CLD_ISOTM']
        ct_slf_noresm = ds['CT_SLF']

        ds['CT_SLF_ISOTM_AVG'] = ds['CT_SLF'].mean(dim = 'time', skipna=True)
        self.case_da = ds
        