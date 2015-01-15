#!/usr/bin/env python

# the above shebang uses a virtualenv if one is running, otherwise the
# system python

"""Short example of using python to perform the environment variable
management, etc. needed to use the I/O API fortran library.

Using python has the advantage of providing access to
stem_pytools.ecampbell300_data_paths, where the location of all the
various STEM runs' data files is listed.

Timothy W. Hilton, UC Merced, 13 January 2015
"""

import os, os.path, sys
import subprocess
from stem_pytools import ecampbell300_data_paths as edp

print('python: {}'.format(sys.executable))

# get the GPP I/O API file for the CASA m15 file
runs = edp.get_runs()
fname_GPP = runs['casa_m15_161'].gpp_path
fname_LRU = os.path.join('/home', 'thilton', 'projects', 'COS (ecampbell3)', 
                         'C4_percentage', 'ISLSCP_C4_1DEG_932_regridded',
                         'C4_pct_124x124.nc')
print('GPP file: {}'.format(fname_GPP))
print('LRU file: {}'.format(fname_LRU))

# if the GPP file exists, use I/O API fortran library to open it,
# print out some information about the variables, and create a new
# file with GPP replaced with a dummy variable
if os.path.exists(fname_GPP) and os.path.exists(fname_LRU):
    os.environ['GPP_INPUT'] = fname_GPP
    os.environ['LRU_INPUT'] = fname_LRU
    os.environ['RATIO_FILE'] = './COS_CO2_ratio_const_1.1.nc'
    subprocess.call('make -f calc_fCOS_C4pct.mk clobber', shell=True)
    subprocess.call('make -f calc_fCOS_C4pct.mk', shell=True)
    subprocess.call('./calc_fCOS_C4pct.x', shell=True)
else:
    print('GPP file not found')




