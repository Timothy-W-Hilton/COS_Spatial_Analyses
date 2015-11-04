import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import numpy as np
import os
import os.path
from datetime import datetime

from timutils import colormap_nlevs
from timutils import midpt_norm
from stem_pytools import calc_drawdown
from stem_pytools import NERSC_data_paths as ndp
from stem_pytools import aqout_postprocess as aqp
from stem_pytools import STEM_mapper

if True:
    print('parsing AQOUT and calculating daily stats')
    GC_bounds_run = ndp.get_Boundaries_runs()['GEOSChem_bounds']
    aqc = aqp.aqout_container(GC_bounds_run.aqout_path)
    aqc.parse(t0=datetime(2008, 7, 1), t1=datetime(2008, 8, 31, 23, 59, 59))
    aqc.sum()
    aqc.calc_stats()
    aqc.stats_to_netcdf(os.path.join(os.getenv('SCRATCH'), 'foo.nc'))

    print('calculating drawdown')
    dd = calc_drawdown.calc_STEM_COS_drawdown(aqc.cos_mean)

print('drawing map')
norm = midpt_norm.MidpointNormalize(midpoint=0.0)

STEM_mapper.Mapper124x124(np.mean(dd[8:, ...], axis=0).squeeze()).draw_map(
    t_str='8 Jul - 31 Aug STEM drawdown',
    cmap=plt.get_cmap('PuOr'),
    norm=norm)
plt.gcf().savefig(os.path.join(os.getenv('HOME'), 'plots', 'dd_map.png'))

cos_t_mean = np.mean(aqc.cos_mean[8:, ...], axis=0).squeeze()

cmap, norm = colormap_nlevs.setup_colormap(
    nlevs=9,
    vmin=cos_t_mean.min(),
    vmax=cos_t_mean.max(),
    cmap=plt.get_cmap('Blues'))

STEM_mapper.Mapper124x124(cos_t_mean[0, ...].squeeze()).draw_map(
    t_str='8 Jul - 31 Aug STEM sfc [COS]',
    cmap=cmap,
    norm=norm)
plt.gcf().savefig(os.path.join(os.getenv('HOME'), 'plots', 'cos_sfc.png'))

STEM_mapper.Mapper124x124(cos_t_mean[-1, ...].squeeze()).draw_map(
    t_str='8 Jul - 31 Aug STEM TOA [COS]',
    cmap=cmap,
    norm=norm)
plt.gcf().savefig(os.path.join(os.getenv('HOME'), 'plots', 'cos_TOA.png'))

plt.close('all')


# from map_grid import assemble_data

# C4_data = os.path.join(os.getenv('HOME'),
#                        'aq_out_data_C4.cpickle')
# constLRU_data = os.path.join(os.getenv('HOME'),
#                              'aq_out_data.cpickle')

# cos_dd, gpp, fCOS = assemble_data(constLRU_data,
#                                   get_GPP=True,
#                                   get_fCOS=True)
# cos_dd_C4, gpp, fCOS = assemble_data(C4_data,
#                                      get_GPP=True,
#                                      get_fCOS=True)
# cos_dd.update(cos_dd_C4)
