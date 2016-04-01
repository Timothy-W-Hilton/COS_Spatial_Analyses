import matplotlib
matplotlib.use('AGG')

import netCDF4
import os
import os.path
import numpy as np
from mpl_toolkits.basemap import maskoceans

from stem_pytools import na_map
from stem_pytools.domain import STEM_Domain
from timutils import midpt_norm


def draw_map():
    """parse the C3/C4 weighted average LRU and plot it over a map of
    North America.

    returns: stem_pytools.na_map.NAMapFigure object containing the map
    """

    # nc = netCDF4.Dataset(os.path.join(os.environ['WORK'], 'Data',
    #                                   'STEM', 'input', 'LRU_from_C4pct.nc'))
    nc = netCDF4.Dataset(os.path.join(os.environ['PROJ'],
                                      'Data',
                                      'C4_percentage',
                                      'LRU_from_C4pct.nc'))
    LRU = nc.variables['LRU'][...].squeeze()
    nc.close()
    d = STEM_Domain()
    stem_lon = d.get_lon()
    stem_lat = d.get_lat()
    fcos_cmap, fcos_norm = midpt_norm.get_discrete_midpt_cmap_norm(
        vmin=1.12,
        vmax=1.84,
        midpoint=1.61)
    lru_map = na_map.NAMapFigure(label_latlon=True,
                                 cb_axis=True,
                                 t_str='LRU from C4 veg pct')

    cm = lru_map.map.pcolor(stem_lon, stem_lat,
                            maskoceans(stem_lon, stem_lat, LRU),
                            cmap=fcos_cmap,
                            latlon=True,
                            norm=fcos_norm)
    cbar = lru_map.fig.colorbar(cm, cax=lru_map.ax_cmap)
    # ticks = np.concatenate([np.linspace(1.12, 1.84, 5), np.array([1.61])])
    # cbar = lru_map.fig.colorbar(cm, cax=lru_map.ax_cmap, ticks=ticks)
    # get rid of white lines in colorbar?
    # http://stackoverflow.com/questions/15003353/why-does-my-colorbar-have-lines-in-it
    cbar.solids.set_edgecolor("face")
    cbar.ax.set_title('LRU')
    return(lru_map)


if __name__ == "__main__":
    c3c4_map = draw_map()
    c3c4_map.fig.savefig('c3c4_map.png')
