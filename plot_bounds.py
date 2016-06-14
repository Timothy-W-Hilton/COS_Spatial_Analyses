"""plot climatological boundary ring cell number vs. vertical cell
number along with a map of the physical location of each boundary ring
cell
"""

import matplotlib
matplotlib.use('AGG')

import os
import os.path

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from timutils import colormap_nlevs
from stem_pytools.na_map import NAMapFigure
from stem_pytools import domain

if __name__ == "__main__":

    # plot the boundary ring cell number vs. vertical cell number
    nc = netCDF4.Dataset(
        os.path.join(os.getcwd(),
                     'ClimatologicalBounds',
                     'climatological_COS_bdy_22levs_124x124.nc'))
    ppbv_2_pptv = 1e3
    cos = nc.variables['CO2_TRACER1'][:].squeeze() * ppbv_2_pptv

    cmap, norm = colormap_nlevs.setup_colormap(vmin=cos.min() - 1,
                                               vmax=cos.max() + 1,
                                               nlevs=20,
                                               cmap=plt.get_cmap('Blues'),
                                               extend='neither')

    d = domain.STEM_Domain()
    d.get_STEMZ_height(wrfheight_fname=os.path.join(os.getenv('SARIKA_INPUT'),
                                                    'wrfheight-124x124-22levs.nc'))
    agl_perim = np.array([domain.get_2d_perimeter(d.agl[z, ...]).mean()
                          for z in range(22)])

    print "new settings"
    fontsz = 14
    matplotlib.rcParams.update({'font.size': fontsz})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 6))
    cm = ax.pcolormesh(np.arange(500), agl_perim, cos,
                       cmap=cmap, norm=norm,
                       linewidth=0, rasterized=True)
    ax.set_xlim([0, 500])
    ax.set_ylim([agl_perim.min(), agl_perim.max()])
    ax.set_ylabel('meters above ground', fontdict={'fontsize': fontsz})
    ax.set_xlabel('lateral boundary index', fontdict={'fontsize': fontsz})
    ax_cb = plt.colorbar(cm, cmap=cmap, norm=norm, ax=ax)
    ax_cb.set_label('[COS] (pptv)', fontdict={'fontsize': fontsz})
    ax_cb.solids.set_rasterized(True)

    fig.tight_layout()
    fig.savefig(os.path.join(os.getenv('HOME'),
                             'plots',
                             'climatological_bounds_NEW.pdf'))
    plt.close(fig)
