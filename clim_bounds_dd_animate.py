import matplotlib
matplotlib.use('AGG')

import os
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from datetime import datetime

from timutils import midpt_norm
from timutils import colorbar_from_cmap_norm
from stem_pytools.calc_drawdown import calc_STEM_COS_drawdown
from stem_pytools.domain import STEM_Domain
from stem_pytools import STEM_parsers as sp
from stem_pytools.na_map import NAMapFigure


def map_init(cos):

    cmap, norm = get_norm_cmap(plt.get_cmap('Blues'))

    fig = plt.figure()
    ax = plt.subplot2grid((2, 12), (0, 0), rowspan=2, colspan=10)
    cbar_ax = plt.subplot2grid((2, 12), (0, 11), rowspan=2, colspan=1)

    na_map = NAMapFigure(map_axis=ax, fast_or_pretty='pretty',
                         t_str='climatological bounds')

    domain = STEM_Domain()
    lat = domain.get_lat()
    lon = domain.get_lon()
    map_data = na_map.map.pcolormesh(lon, lat,
                                     cos[0, 0, :, :].squeeze(),
                                     cmap=cmap,
                                     norm=norm,
                                     latlon=True)

    t_idx = 0
    colorbar_from_cmap_norm.colorbar_from_cmap_norm(
        cmap, norm, cbar_ax, None, cos[t_idx, ...])
    cbar_ax.set_title('ppt')
    cbar_ax = None
    return(na_map, fig, ax, cbar_ax, map_data)


def get_norm_cmap(cmap_arg=plt.get_cmap('Blues')):
    cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(vmin=-85,
                                                         vmax=85,
                                                         midpoint=0.0,
                                                         bands_above_mdpt=7,
                                                         bands_below_mdpt=7,
                                                         extend='both')
    return(cmap, norm)


def map_update(i, m, fig, ax, cbar_ax, t_idx, cos, map_data):

    z_sfc = 0  # array index of surface
    map_data.set_array(cos[i, z_sfc, :-1, :-1].squeeze().flatten())
    # ax[0].set_xlim([-160.0, -40.0])
    # ax[0].set_ylim([15.0, 85])

    t_str = '{}'.format(t_idx[i])
    ax.set_title(t_str)

    print('plotting for {}'.format(t_str))
    fig.canvas.draw()


if __name__ == "__main__":
    hours_in_day = 24
    data_dir = os.path.join('/project', 'projectdirs', 'm2319',
                            'STEM_Runs',
                            'STEM_NAmerica_Climatological_Bounds', 'output')
    outfile = 'climatological_bounds.mp4'

    cos = sp.parse_STEM_var(nc_fname=os.path.join(
        data_dir, 'AQOUT.climatological_bnd.nc'),
        t0=datetime(2008, 7, 1),
        t1=datetime(2008, 8, 31, 23, 59, 59),
        varname='CO2_TRACER1')
    dd = calc_STEM_COS_drawdown(cos['data'])
    m, fig, ax, cbar_ax, map_data = map_init(dd)

    im_ani = animation.FuncAnimation(fig,
                                     func=map_update,
                                     frames=len(cos['t']),
                                     interval=100,
                                     # blit=True,  # only update parts that changed
                                     fargs=[m, fig, ax, cbar_ax, cos['t'],
                                            dd, map_data])

    im_ani.save(outfile, metadata={'artist': 'STEM'}, bitrate=100000)
    plt.close(fig)
