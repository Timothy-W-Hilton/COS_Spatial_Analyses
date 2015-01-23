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
import argparse
import datetime 

from stem_pytools import ecampbell300_data_paths as edp
from stem_pytools import na_map
from stem_pytools import STEM_parsers as sp
from timutils import midpt_norm

#============================================================
def fCOS_from_C4pct_diagnostics():

    #draw_LRU_map()

    runs = edp.get_C3C4runs()
    for k in runs.keys():
        fcos = sp.parse_STEM_var(nc_fname=runs[k].fcos_path,
                                 t0 = datetime.datetime(2008,7,1),
                                 t1=datetime.datetime(2008,8,31,23,59,59),
                                 varname='cos')
        print('{}: min: {:0.2e}    max: {:0.2e}    mean: {:0.2e}'.format(k,
                                                          fcos['data'].min(),
                                                          fcos['data'].max(),
                                                          fcos['data'].mean()))
        fig, ax = plt.subplots()
        cm = ax.pcolor(fcos['data'].squeeze().mean(axis=0), 
                       cmap=plt.get_cmap('Blues'))
        plt.colorbar(cm)
        plt.title(k)
        fig.savefig('/tmp/fcos_{}.png'.format(k))

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

#parse arguments
parser = argparse.ArgumentParser(description=("""script to calculate fCOS for all GPP products defined in stem_pytools.ecampbell300_data_paths."""))
parser.add_argument('--diagnostics',
                    dest='run_diagnostics',
                    action='store_true',
                    help=('if set, some diagnostic plots are created'
                          'after the calculations.'))
args = parser.parse_args()

# compile the fortran code
subprocess.call('make -f calc_fCOS_C4pct.mk clobber', shell=True)
subprocess.call('make -f calc_fCOS_C4pct.mk', shell=True)

fname_C4pct = os.path.join('/home', 'thilton', 'projects', 'COS (ecampbell3)', 
                         'C4_percentage', 'ISLSCP_C4_1DEG_932_regridded',
                         'C4_pct_124x124.nc')
# get the GPP I/O API file for the CASA m15 file
runs = edp.get_C3C4runs()
for this_run in runs.values():
    fname_GPP = this_run.gpp_path
    sys.stdout.write('\n\nGPP file: {}\n'.format(fname_GPP))
    sys.stdout.write('C4 file: {}\n\n'.format(fname_C4pct))
    sys.stdout.flush()

    if os.path.exists(fname_GPP) and os.path.exists(fname_C4pct):
        os.environ['GRIDDESC'] = os.path.join(os.environ['HOME'], 'Data', 
                                              'STEM', 'input', 'GRIDDESC.txt')
        os.environ['GPP_INPUT'] = fname_GPP
        os.environ['C4pct_INPUT'] = fname_C4pct
        os.environ['RATIO_FILE'] = './COS_CO2_ratio_const_1.1.nc'
        os.environ['LRU_FILE'] = './LRU_from_C4pct.nc'
        os.environ['fCOS_FILE'] = './fCOS_{}_2008_124x124_LRUfromC4pct.nc'.format(
            re.sub('[\ \-]', '', this_run.model))
        # run the fortran part
        if this_run.model.lower().find('casa') >= 0:
            t_step = 30000 # 3 hours expressed as HHMMSS
        else:
            t_step = 7320000 # 30.5 days expressed as HHMMSS
        subprocess.call(
            './calc_fCOS_C4pct.x {} {}'.format(
                re.sub(' ', '\ ', this_run.model),
                t_step),
            shell=True)

        #move the new fCOS file to its permanent location
        newname = os.path.join(os.path.dirname(this_run.fcos_path),
                               os.path.basename(os.environ['fCOS_FILE']))
        sys.stdout.write('\n{} --> {}\n\n'.format(os.environ['fCOS_FILE'], 
                                                    newname))
        sys.stdout.flush()
        os.rename(os.environ['fCOS_FILE'], newname)

    else:
        print('GPP file not found')

if args.run_diagnostics:
    print 'creating diagnostic plots'
    fCOS_from_C4pct_diagnostics()
