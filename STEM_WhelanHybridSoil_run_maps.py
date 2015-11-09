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

# I put the Whelan-Kettle hybrid soil fluxes through STEM in pmol m-2
# s-1 when it was expecting mol m-2 s-1.  So the AQOUT concentrations
# are too high be a factor of 1e12.
mol_pmol_correction = 1e-12

if True:
    print('parsing AQOUT and calculating daily stats')
    hybrid_Fsoil_run = ndp.get_runs()['Fsoil_Hybrid5Feb']
    aqc = aqp.aqout_container(hybrid_Fsoil_run.aqout_path)
    aqc.parse(t0=datetime(2008, 7, 1), t1=datetime(2008, 8, 31, 23, 59, 59))

    aqc.data[0] = aqc.data[0] * mol_pmol_correction
    aqc.sum()
    aqc.calc_stats()
    aqc.stats_to_netcdf(os.path.join(os.getenv('SCRATCH'),
                                     'hybrid_Fsoil_STEM_run_daily.nc'))

    print('calculating drawdown')
    dd = calc_drawdown.calc_STEM_COS_drawdown(aqc.cos_mean)

print('drawing map')
norm = midpt_norm.MidpointNormalize(midpoint=0.0)
f_or_p = 'pretty'
mcm3_2_pptv = 1e12  # factor to molecules cm-3 to parts per trillion by volume
STEM_mapper.Mapper124x124(np.mean(dd[:], axis=0).squeeze()).draw_map(
    fast_or_pretty=f_or_p,
    t_str='1 Jul - 31 Aug mean STEM drawdown, Whelan-Kettle "hybrid" Fsoil',
    cmap=plt.get_cmap('PuOr'),
    norm=norm)
plt.gcf().savefig(os.path.join(os.getenv('HOME'), 'plots', 'dd_map_Fsoil.png'))

# cos_t_mean = np.mean(aqc.cos_mean[8:, ...], axis=0).squeeze() * mcm3_2_pptv

# cmap, norm = colormap_nlevs.setup_colormap(
#     nlevs=9,
#     vmin=cos_t_mean.min(),
#     vmax=cos_t_mean.max(),
#     cmap=plt.get_cmap('Blues'))

# STEM_mapper.Mapper124x124(cos_t_mean[0, ...].squeeze()).draw_map(
#     fast_or_pretty=f_or_p,
#     t_str='8 Jul - 31 Aug STEM sfc [COS]',
#     cmap=cmap,
#     norm=norm)
# plt.gcf().savefig(os.path.join(os.getenv('HOME'), 'plots', 'cos_sfc.png'))

# STEM_mapper.Mapper124x124(cos_t_mean[-1, ...].squeeze()).draw_map(
#     fast_or_pretty=f_or_p,
#     t_str='8 Jul - 31 Aug STEM TOA [COS]',
#     cmap=cmap,
#     norm=norm)
# plt.gcf().savefig(os.path.join(os.getenv('HOME'), 'plots', 'cos_TOA.png'))

# plt.close('all')
