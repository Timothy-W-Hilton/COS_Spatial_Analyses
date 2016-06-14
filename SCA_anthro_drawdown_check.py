import matplotlib
matplotlib.use('AGG')

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from datetime import datetime
from stem_pytools import STEM_parsers as sp
from stem_pytools import calc_drawdown
from timutils import colormap_nlevs

# STEM 124x124 domain grid coordinates of NOAA sites
#                   stem_x  stem_y  climatological boundaries  Anthropogenic, Zumkehr  Anthropogenic, Zumkehr, clim
# sample_site_code
# AAO                   48      43                   7.003995               -6.583751                      0.461426
# BNE                   60      36                   5.565293               -5.113157                      0.488139
# CAR                   70      30                  -2.589339               -3.943283                     -6.520686
# CMA                   29      60                   5.623783               -6.152536                     -0.501340
# DND                   70      48                   8.433928               -2.581486                      5.913424
# ESP                  108      39                  11.127015               -1.063064                     10.118439
# ETL                   84      56                   6.070864               -1.976531                      4.149570
# HIL                   46      44                   6.796292               -6.721734                      0.041312
# HIP                  123     123                  18.194332               -0.041563                     18.178090
# LEF                   58      52                   6.819895               -3.918529                      2.920137
# NHA                   33      70                   2.750400               -5.739810                     -2.985367
# OIL                   50      45                   7.749532               -5.901172                      1.836695
# SCA                   23      44                 -12.922313               -7.245397                    -20.285646
# SGP                   56      27                  -1.796552               -5.065788                     -6.812202
# TGC                   42      10                 -16.622076                0.427395                    -15.967099
# THD                  103      21                  14.693239               -1.615453                     13.104427
# ULB                  123     123                  18.194332               -0.041563                     18.178090
# WBI                   54      43                   8.785475               -5.440924                      3.457975


def plot_site_dd(dd, sitename, x, y):
    this_dd = dd[..., x, y].squeeze()
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 8))
    ax.plot(this_dd)
    ax.set_xlabel('hours since 1 July 2008 00:00')
    ax.set_ylabel('surface [COS] drawdown, ppt')
    ax.text(1200, 0.9 * this_dd.max(),
            'mean (10 Jul - 31 Aug): {:0.1f}'.format(
                this_dd[240:].squeeze().mean()))
    ax.text(1200, 0.8 * this_dd.max(),
            'median (10 Jul - 31 Aug): {:0.1f}'.format(
                np.median(this_dd[240:].squeeze())))
    ax.set_title(sitename)
    fig.savefig('/global/homes/t/twhilton/plots/{}_dd.pdf'.format(sitename))

cos_clim = sp.parse_STEM_var(nc_fname='/project/projectdirs/m2319/STEM_Runs/STEM_NAmerica_Climatological_Bounds/output/AQOUT.climatological_bnd.nc',
                             t0=datetime(2008, 7, 1),
                             t1=datetime(2008, 8, 31, 23, 59, 59),
                             varname='CO2_TRACER1')

nha_x = 33
nha_y = 70

cma_x = 29
cma_y = 60

sca_x = 23
sca_y = 44

molecules_m3_to_ppt = 1e12
cos_sca = cos_clim['data'][:, :, sca_x, sca_y] * molecules_m3_to_ppt

cmap, norm = colormap_nlevs.setup_colormap(vmin=cos_sca.min() - 1,
                                           vmax=cos_sca.max() + 1,
                                           nlevs=20,
                                           cmap=plt.get_cmap('Blues'),
                                           extend='neither')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 8))
cm = ax.pcolormesh(np.transpose(cos_sca), cmap=cmap, norm=norm,
                   linewidth=0, rasterized=True)
ax.set_xlim([0, 1488])
ax.set_ylim([0, 22])
ax.set_ylabel('STEM Z level')
ax.set_xlabel('hours since 1 July 2008 00:00')
ax.set_title('SCA [COS], climatological bounds STEM run')
ax_cb = plt.colorbar(cm, cmap=cmap, norm=norm, ax=ax)
ax_cb.set_label('[COS] (ppt)')
ax_cb.solids.set_rasterized(True)
fig.savefig('/global/homes/t/twhilton/plots/SCA_COS.png')
plt.close(fig)

dd = calc_drawdown.calc_STEM_COS_drawdown(cos_clim['data'])
plot_site_dd(dd, 'NHA', nha_x, nha_y)
plot_site_dd(dd, 'CMA', cma_x, cma_y)
plot_site_dd(dd, 'SCA', sca_x, sca_y)
plt.close(fig)
