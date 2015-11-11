import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import numpy as np
import numpy.ma as ma
import os
import os.path
from datetime import datetime

from timutils import colormap_nlevs
from timutils import midpt_norm
from stem_pytools import calc_drawdown
from stem_pytools import NERSC_data_paths as ndp
from stem_pytools import aqout_postprocess as aqp
from stem_pytools import STEM_mapper
from stem_pytools import noaa_ocs

# I put the Whelan-Kettle hybrid soil fluxes through STEM in pmol m-2
# s-1 when it was expecting mol m-2 s-1.  So the AQOUT concentrations
# are too high be a factor of 1e12.
mol_pmol_correction = 1e-12

if True:
    print('parsing AQOUT and calculating daily stats')
    hybrid_Fsoil_run = ndp.get_runs()['Fsoil_Hybrid5Feb']
    aqc = aqp.aqout_container(hybrid_Fsoil_run.aqout_path)
    aqc.parse(t0=datetime(2008, 7, 1), t1=datetime(2008, 8, 31, 23, 59, 59))


    import pdb; pdb.set_trace()
    aqc.data[0] = aqc.data[0] * mol_pmol_correction
    aqc.sum()
    aqc.calc_stats()
    aqc.stats_to_netcdf(os.path.join(os.getenv('SCRATCH'),
                                     'hybrid_Fsoil_STEM_run_daily.nc'))

    print('calculating drawdown')
    dd = calc_drawdown.calc_STEM_COS_drawdown(aqc.cos_mean)

print('drawing map')
cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(vmin=-40, vmax=20,
                                                     midpoint=0.0,
                                                     bands_above_mdpt=3,
                                                     bands_below_mdpt=7)
f_or_p = 'fast'
mcm3_2_pptv = 1e12  # factor to molecules cm-3 to parts per trillion by volume
STEM_mapper.Mapper124x124(np.mean(dd[:], axis=0).squeeze()).draw_map(
    fast_or_pretty=f_or_p,
    t_str=('1 Jul - 31 Aug 2008 mean STEM drawdown, '
           'Whelan-Kettle "hybrid" Fsoil'),
    cmap=cmap,
    norm=norm)
plt.gcf().savefig(os.path.join(os.getenv('HOME'), 'plots', 'dd_map_Fsoil.png'))


[agl, asl] = noaa_ocs.get_STEMZ_height(
    topo_fname=os.path.join(os.getenv('SARIKA_INPUT'), 'TOPO-124x124.nc'),
    wrfheight_fname=os.path.join(os.getenv('SARIKA_INPUT'),
                                 'wrfheight-124x124-2008-2009-22levs.nc'))
t_steps_per_day = 4
jul1 = (2008183 - 2008060) * t_steps_per_day
aug31 = jul1 + (62 * t_steps_per_day)
agl_mean = np.mean(agl[jul1:aug31, ...], axis=0)
agl_mean = np.broadcast_arrays(agl_mean, aqc.cos_mean[8:, ...])[0]

free_trop = ma.masked_where(agl_mean >= 4000, aqc.cos_mean[8:, ...])
free_trop_mean = mcm3_2_pptv * np.mean(np.mean(free_trop,
                                               axis=0),
                                       axis=0).squeeze()

abl = ma.masked_where(agl_mean <= 2000, aqc.cos_mean[8:, ...])
abl_mean = mcm3_2_pptv * np.mean(np.mean(abl, axis=0), axis=0).squeeze()

cmap, norm = colormap_nlevs.setup_colormap(
    nlevs=9,
    vmin=np.dstack((abl_mean, free_trop_mean)).min(),
    vmax=np.dstack((abl_mean, free_trop_mean)).max(),
    cmap=plt.get_cmap('Blues'))

STEM_mapper.Mapper124x124(abl_mean.squeeze()).draw_map(
    fast_or_pretty=f_or_p,
    t_str='8 Jul - 31 Aug STEM [COS] <= 2000 m; GEOS-Chem boundaries',
    cmap=cmap,
    norm=norm)
plt.gcf().savefig(os.path.join(os.getenv('HOME'), 'plots', 'cos_ABL.png'))

STEM_mapper.Mapper124x124(free_trop_mean.squeeze()).draw_map(
    fast_or_pretty=f_or_p,
    t_str='8 Jul - 31 Aug STEM [COS] >= 4000 m; GEOS-Chem boundaries',
    cmap=cmap,
    norm=norm)
plt.gcf().savefig(os.path.join(os.getenv('HOME'),
                               'plots', 'cos_free_trop.png'))

plt.close('all')
