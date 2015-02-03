"""script to plot the [COS] along the Western boundary of the STEM
domain.  Useful for comparing boundary conditions.
"""


import matplotlib.pyplot as plt
import numpy as np
import os.path
import os
from datetime import datetime
from stem_pytools import STEM_parsers as sp
from stem_pytools.noaa_ocs import get_STEMZ_height
from timutils import midpt_norm
import brewer2mpl

bd_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                        'bdv-124x124-cos-pctm_2008_2009.nc')
topo_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                          'TOPO-124x124.nc')
wrf_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                         'wrfheight-124x124-22levs.nc')

bd = sp.parse_STEM_var(bd_fname,
                       varname='CO2_TRACER1',
                       t0=datetime(2008, 7, 1),
                       t1=datetime(2008, 9, 1))

w_bnd_tavg = np.mean(bd['data'][:, :, -125:-1], axis=0) * 1000

stemlon, stemlat, topo = sp.parse_STEM_coordinates(topo_fname)
agl, asl = get_STEMZ_height(topo_fname, wrf_fname)

cos_norm = midpt_norm.MidpointNormalize(midpoint=450)
cos_cmap = brewer2mpl.get_map('PuOr', 'diverging', 8).mpl_colormap

fig, ax = plt.subplots()
lat = stemlat[:, -1]
ax.set_xlim([lat.min(), lat.max()])
cm = ax.pcolormesh(lat, agl[:, :, -1],
                   w_bnd_tavg,
                   cmap=plt.get_cmap('PuOr'),
                   norm=cos_norm)
ax.set_xlabel('latitude ($^\circ$N)')
ax.set_ylabel('height AGL (m)')
ax.set_title('PCTM Western boundary July/August mean')
ticks = np.round(np.concatenate([np.linspace(w_bnd_tavg.min(), 450, 4),
                                 np.linspace(450, w_bnd_tavg.max(), 6)]))
cbar = fig.colorbar(cm, ax=ax,
                    ticks=ticks)
cbar.ax.set_title('[COS] (ppbv)\n')
fig.savefig('/tmp/STEM_W_bnd_JulAug.pdf')
plt.close('all')
