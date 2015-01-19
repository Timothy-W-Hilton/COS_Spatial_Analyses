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
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import brewer2mpl
from stem_pytools import ecampbell300_data_paths as edp
from stem_pytools import na_map
from stem_pytools import STEM_parsers as sp
from timutils import midpt_norm

#============================================================
def fCOS_from_C4pct_diagnostics():

    draw_LRU_map()

#============================================================
def draw_LRU_map():
    nc = netCDF4.Dataset(os.environ['LRU_FILE'])
    LRU = nc.variables['LRU'][...].squeeze()
    nc.close()
    stem_lon, stem_lat, topo = sp.parse_STEM_coordinates(
        os.path.join(os.getenv('SARIKA_INPUT'), 'TOPO-124x124.nc'))
    fcos_norm = midpt_norm.MidpointNormalize(midpoint=1.61, vmin=LRU.min(), vmax=LRU.max())
    fcos_cmap = brewer2mpl.get_map('PuOr', 'diverging', 8).mpl_colormap
    lru_map = na_map.NAMapFigure(cb_axis=True,t_str='LRU from C4 veg pct')
    cm = lru_map.map.pcolor(stem_lon, stem_lat, LRU,
                            cmap=fcos_cmap,
                            latlon=True, 
                            norm=fcos_norm)
    plt.colorbar(cm, cax=lru_map.ax_cmap)
    lru_map.fig.savefig('LRU_from_c4pct.png')

run_diagnostics = True

# compile the fortran code
subprocess.call('make -f calc_fCOS_C4pct.mk clobber', shell=True)
subprocess.call('make -f calc_fCOS_C4pct.mk', shell=True)

fname_C4pct = os.path.join('/home', 'thilton', 'projects', 'COS (ecampbell3)', 
                         'C4_percentage', 'ISLSCP_C4_1DEG_932_regridded',
                         'C4_pct_124x124.nc')
# get the GPP I/O API file for the CASA m15 file
runs = edp.get_runs()
for this_run in runs.values()[0:1]:
    fname_GPP = this_run.gpp_path
    print('GPP file: {}'.format(fname_GPP))
    print('C4 file: {}'.format(fname_C4pct))

    if os.path.exists(fname_GPP) and os.path.exists(fname_C4pct):
        os.environ['GRIDDESC'] = os.path.join(os.environ['HOME'], 'Data', 
                                              'STEM', 'input', 'GRIDDESC.txt')
        os.environ['GPP_INPUT'] = fname_GPP
        os.environ['C4pct_INPUT'] = fname_C4pct
        os.environ['RATIO_FILE'] = './COS_CO2_ratio_const_1.1.nc'
        os.environ['LRU_FILE'] = './LRU_from_C4pct.nc'
        os.environ['fCOS_FILE'] = './fCOS_from_{}_C4pct.nc'.format(
            re.sub('[\ \-]', '', this_run.model))
        # run the fortran part
        subprocess.call(
            './calc_fCOS_C4pct.x {}'.format(re.sub(' ', '\ ', this_run.model)),
            shell=True)
    else:
        print('GPP file not found')

if run_diagnostics:
    fCOS_from_C4pct_diagnostics()
