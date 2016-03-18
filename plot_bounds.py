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
    fontsz = 14
    matplotlib.rcParams.update({'font.size': fontsz})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    cm = ax.pcolormesh(cos, cmap=cmap, norm=norm,
                       linewidth=0, rasterized=True)
    ax.set_xlim([0, 500])
    ax.set_ylim([0, 22])
    ax.set_ylabel('STEM Z level', fontdict={'fontsize': fontsz})
    ax.set_xlabel('lateral boundary index',
                  fontdict={'fontsize': fontsz})
    # ax.set_title('climatological N American boundary',
    #              fontdict={'fontsize': fontsz})
    cb = plt.colorbar(cm, cmap=cmap, norm=norm, ax=ax, format='%d')
    cb.set_label('[COS] (pptv)', fontdict={'fontsize': fontsz})
    cb.ax.tick_params(labelsize=fontsz)
    cb.solids.set_rasterized(True)

    fig.savefig(os.path.join(os.getenv('HOME'),
                             'plots',
                             'climatological_bounds.pdf'))
    plt.close(fig)
