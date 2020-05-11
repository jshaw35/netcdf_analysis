# netcdf_analysis
General files for processing netcdf format from ESM output.

cplot.py contains functions for grabbing and storing variables (CDFextract) and plotting (cplot) that data.

Other files give examples of plotting.

'cloudtop' files work with CALIOP observations from Olimpia Bruno at KIT. All other slf-related files use data from Trude.

Commonly used shared files are (names are self-explanatory):
functions.py
imports.py
classes.py

I am trying to move old files to 'Depreciated' in order to keep the main directory clean.

Notebooks used for developing a function or technique should be moved to 'Testbeds'.

Also trying to figure out if moving paired .py scripts to a separate folder will prevent them from staying updated with the development notebook. This doesn't seem to work right now and the files move to stay in the same directory as each other.