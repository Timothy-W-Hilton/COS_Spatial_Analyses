"""Module to create an animated map of STEM results for North America
from July and August climatological mean [COS] boundary conditions
with no surface fluxes.

"""

import matplotlib
matplotlib.use('AGG')

import os
import sys
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from datetime import datetime

from timutils import midpt_norm
from timutils import colorbar_from_cmap_norm
from timutils.ffmpeg_tester import check_for_ffmpeg
from stem_pytools.calc_drawdown import calc_STEM_COS_drawdown
from stem_pytools.domain import STEM_Domain
from stem_pytools import STEM_parsers as sp
from stem_pytools.na_map import NAMapFigure


def map_init(cos):
    """
    Initialize a plot showing the STEM simulated [COS] on a map of N
    America for 1 July 2008 00:00 (the first time step of the run).
    The map provides a plotting framework to subsequently update with
    [COS] for later timesteps.

    ARGS:
    cos (numpy.ndarray): 4-D array of [COS] values in parts per
       trillion by volume (ppt).  The dimensions are [t, Z, X, Y].

    RETURNS:
    four-element tuple containing:
    na_map (stem_pytools.na_map.NAMapFigure): NAMapFigure that
       encapsulates the map
    fig (matplotlib.figure.Figure): Figure object containing the plot
    ax (matplotlib.axes.Axes): Axes object containing the map
    cbar_ax (matplotlib.axes.Axes): Axes object containing the colorbar
    map_data ( matplotlib.collections.QuadMesh): output of
       matplotlib.pyplot.pcolormesh containing the data plotted on the
       map
    """
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
    return na_map, fig, ax, cbar_ax, map_data


def get_norm_cmap(cmap_arg=plt.get_cmap('Blues')):
    """returns a colormap and normalizer instance to produce a
    14-element colorbar centered anchored to [-30.0, 0.0, 30.0] based
    on a specified colormap.  It makes the most sense to use a
    diverging colorscale for this because the function is designed to
    produce a colorscale that emphasizes values above and below 0.0.
    The colorscale will have seven intervals above 0.0 and seven
    below.

    ARGS:
    cmap_arg (matplotlib.colors.LinearSegmentedColormap): colormap to
       base the colorscale from.

    RETURNS:
    two-element tuple containing cmap
    (matplotlib.colors.ListedColormap) and norm
    (matplotlib.colors.BoundaryNorm)
    """
    cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(vmin=-30,
                                                         vmax=30,
                                                         midpoint=0.0,
                                                         bands_above_mdpt=7,
                                                         bands_below_mdpt=7,
                                                         extend='both')
    return cmap, norm


def map_update(i, m, fig, ax, cbar_ax, t_idx, cos, map_data):
    """update the map plot with a new two-dimensional field of data

    ARGS:
    i (int): index of the timestep whose data are to be plotted
    m:
    fig (matplotlib.figure.Figure): the figure object containing the plot
    ax: (matplotlib.axes.Axes): axes containing the main map figure
    cbar_ax: (matplotlib.axes.Axes): axes containing the colorbar
    t_idx (numpy.ndarray): array of datetime.datetime objects
       containing the sequence of model timesteps in the animation
    cos (numpy.ndarray): array of [COS] concentrations in parts per
       trillion by volume (ppt), indexed by [t, Z, X, Y].
    map_data (matplotlib.collections.QuadMesh): QuadMesh object
       containing the patches on the map.  This is part of the output of
       map_init()

    RETURNS:
    None
    """
    z_sfc = 0  # array index of surface
    map_data.set_array(cos[i, z_sfc, :-1, :-1].squeeze().flatten())
    # ax[0].set_xlim([-160.0, -40.0])
    # ax[0].set_ylim([15.0, 85])

    t_str = '{}'.format(t_idx[i])
    ax.set_title(t_str)

    print('plotting for {}'.format(t_str))
    fig.canvas.draw()


if __name__ == "__main__":

    ffmpeg_present = check_for_ffmpeg()
    if ffmpeg_present is False:
        sys.exit("ffmpeg not found. Exiting now.")

    hours_in_day = 24
    data_dir = os.path.join('/project', 'projectdirs', 'm2319',
                            'STEM_Runs',
                            'STEM_NAmerica_Climatological_Bounds', 'output')
    outfile = 'climatological_bounds_TEST.mp4'

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
