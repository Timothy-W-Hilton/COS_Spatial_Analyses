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
import re
from stem_pytools import ecampbell300_data_paths as edp
import netCDF4

print('python: {}'.format(sys.executable))

# get the GPP I/O API file for the CASA m15 file
runs = edp.get_runs()
this_run = runs['casa_m15_161']
fname_GPP = this_run.gpp_path
fname_C4pct = os.path.join('/home', 'thilton', 'projects', 'COS (ecampbell3)', 
                         'C4_percentage', 'ISLSCP_C4_1DEG_932_regridded',
                         'C4_pct_124x124.nc')
print('GPP file: {}'.format(fname_GPP))
print('C4 file: {}'.format(fname_C4pct))

# if the GPP file exists, use I/O API fortran library to open it,
# print out some information about the variables, and create a new
# file with GPP replaced with a dummy variable

if os.path.exists(fname_GPP) and os.path.exists(fname_C4pct):
    os.environ['GPP_INPUT'] = fname_GPP
    os.environ['C4pct_INPUT'] = fname_C4pct
    os.environ['RATIO_FILE'] = './COS_CO2_ratio_const_1.1.nc'
    os.environ['LRU_FILE'] = './LRU_from_C4pct.nc'
    os.environ['fCOS_FILE'] = './fCOS_from_{}_C4pct.nc'.format(
        re.sub('[\ \-]', '', this_run.model))
    # compile and run the fortran part
    subprocess.call('make -f calc_fCOS_C4pct.mk clobber', shell=True)
    subprocess.call('make -f calc_fCOS_C4pct.mk', shell=True)
    subprocess.call('./calc_fCOS_C4pct.x {}'.format(re.sub(' ', '\ ', this_run.model)), 
                    shell=True)

    nc = netCDF4.Dataset(os.environ['LRU_FILE'])
    lru = nc.variables['LRU'][...].squeeze()
    nc.close()
    nc = netCDF4.Dataset(os.environ['C4pct_INPUT'])
    c4pct = nc.variables['C4pct'][...].squeeze()
    nc.close()

else:
    print('GPP file not found')


