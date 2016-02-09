"""plot climatological boundary ring cell number vs. vertical cell
number along with a map of the physical location of each boundary ring
cell
"""

import matplotlib.pyplot as plt
import netCDF4
from timutils import colormap_nlevs
from stem_pytools import calc_drawdown
from stem_pytools.na_map import NAMapFigure
from stem_pytools import domain

if __name__ == "__main__":

    # plot the boundary ring cell number vs. vertical cell number
    nc = netCDF4.Dataset('/Users/tim/work/Data/STEM/input/climatological_COS_bdy_22levs_124x124.nc')
    ppbv_2_pptv = 1e3
    cos = nc.variables['CO2_TRACER1'][:].squeeze() * ppbv_2_pptv

    cmap, norm = colormap_nlevs.setup_colormap(vmin=cos.min() - 1,
                                               vmax=cos.max() + 1,
                                               nlevs=20,
                                               cmap=plt.get_cmap('Blues'),
                                               extend='neither')

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
    cm = ax[0].pcolormesh(cos, cmap=cmap, norm=norm,
                          linewidth=0, rasterized=True)
    ax[0].set_xlim([0, 500])
    ax[0].set_ylim([0, 22])
    ax[0].set_ylabel('STEM Z level')
    ax[0].set_xlabel('boundary cell (0 is southernmost corner, proceeds CCW)')
    ax[0].set_title('climatological N American boundary')
    ax_cb = plt.colorbar(cm, cmap=cmap, norm=norm, ax=ax[0])
    ax_cb.set_label('[COS] (pptv)')
    ax_cb.solids.set_rasterized(True)

    # plot boundary ring cell numbers on a map to illustrate their
    # physical locations
    m = NAMapFigure(t_str='boundary cell locations',
                    map_axis=ax[1],
                    mapwidth=9.5e6)
    d = domain.STEM_Domain()
    lat = domain.get_2d_perimeter(d.get_lat())
    lon = domain.get_2d_perimeter(d.get_lon())
    x, y = m.map(lon, lat)
    for i in range(0, lon.size, 50):
        m.ax_map.text(x[i], y[i],  str(i), color='red',
                      size='large', weight='bold')
    fig.savefig('/Users/tim/Desktop/climatological_bounds.pdf')
    plt.close(fig)
