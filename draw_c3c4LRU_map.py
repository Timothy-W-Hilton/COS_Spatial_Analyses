import matplotlib.pyplot as plt
import netCDF4
import os
import os.path
import numpy as np

from stem_pytools import na_map
from stem_pytools import noaa_ocs
from stem_pytools import STEM_parsers as sp
from timutils import midpt_norm
from spatial_analysis_utilities import get_noaa_COS_data_path


def draw_map():
    """parse the C3/C4 weighted average LRU and plot it over a map of
    North America.

    returns: stem_pytools.na_map.NAMapFigure object containing the map
    """

    nc = netCDF4.Dataset(os.path.join(os.environ['WORK'], 'Data',
                                      'STEM', 'input', 'LRU_from_C4pct.nc'))
    LRU = nc.variables['LRU'][...].squeeze()
    nc.close()
    stem_lon, stem_lat, topo = sp.parse_STEM_coordinates(
        os.path.join(os.getenv('SARIKA_INPUT'), 'TOPO-124x124.nc'))
    fcos_norm = midpt_norm.MidpointNormalize(midpoint=1.61)
    lru_map = na_map.NAMapFigure(cb_axis=True, t_str='LRU from C4 veg pct')
    cm = lru_map.map.pcolor(stem_lon, stem_lat, LRU,
                            cmap=plt.get_cmap('PuOr_r'),
                            latlon=True,
                            norm=fcos_norm)
    ticks = np.concatenate([np.linspace(1.12, 1.84, 5), np.array([1.61])])
    cbar = lru_map.fig.colorbar(cm, cax=lru_map.ax_cmap, ticks=ticks)
    # get rid of white lines in colorbar?
    # http://stackoverflow.com/questions/15003353/why-does-my-colorbar-have-lines-in-it
    cbar.solids.set_edgecolor("face")
    cbar.ax.set_title('LRU')
    return(lru_map)
