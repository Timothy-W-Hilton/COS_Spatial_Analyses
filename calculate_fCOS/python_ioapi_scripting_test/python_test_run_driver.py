#!/usr/bin/env python

# the above shebang uses a virtualenv if one is running, otherwise the
# system python

"""
Short example of using python to perform the environment variable
management, etc. needed to use the I/O API fortran library.

Timothy W. Hilton, UC Merced, 13 January 2015 
"""

import os, os.path, sys
import subprocess
from stem_pytools import ecampbell300_data_paths as edp

print('using {}'.format(sys.executable))

# get the GPP I/O API file for the CASA m15 file
runs = edp.get_runs()
fname_gpp = runs['casa_m15_161'].gpp_path
print('GPP file: {}'.format(fname_gpp))
# if the GPP file exists, use I/O API fortran library to open it,
# print out some information about the variables, and create a new
# file with GPP replaced with a dummy variable
if os.path.exists(fname_gpp):
    os.environ['GPP_INPUT'] = fname_gpp
    os.environ['OUTPUT'] = './test_out.nc'
    subprocess.call('make -f python_script_ioapi_test.mk clobber', shell=True)
    subprocess.call('make -f python_script_ioapi_test.mk', shell=True)
    subprocess.call('./python_script_ioapi_test.x', shell=True)
else:
    print('GPP file not found')




