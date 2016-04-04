import matplotlib
matplotlib.use('AGG')

import matplotlib.pyplot as plt
import netCDF4
import os
import os.path
import numpy as np
from mpl_toolkits.basemap import maskoceans
import matplotlib.gridspec as gridspec

from stem_pytools import na_map
from stem_pytools.domain import STEM_Domain
from timutils import midpt_norm
import timutils.io


def setup_panel_array(nrows=3, ncols=6, figsize=(6, 3)):
    """create a matrix of axes for maps with individual color bars

    create a figure containing a matrix of axes with nrows rows and
    ncols columns, where each "panel" consists of a main axis using
    90% of the horizontal space and a secondary colorbar axis
    occupying 10% of the horizontal space.

    ARGS
    nrows (int): number of rows of axes
    ncols (int): number of columns of axes
    figsize ((float, float)): width and height of figure in inches

    OUTPUTS
    fig: matplotlib.figure.Figure object
    ax: nrows by ncols numpy array of
        matplotlib.axes._subplots.AxesSubplot objects
    cbar_ax: nrows by ncols numpy array of
        matplotlib.axes._subplots.AxesSubplot objects (for colorbars)
    """
    fig = plt.figure(figsize=figsize)
    # two gridspecs - one for maps, one for colorbasr

    gs = gridspec.GridSpec(
        nrows,
        ncols * 2,
        width_ratios=[9, 1] * ncols)
    # gs_maps = gs[range(0, ncols * 2, 2), range(nrows)]
    # gs_maps = gs[range(1, ncols * 2, 2), range(nrows)]
    # gs_maps.update(hspace=0.01, wspace=0.0, left=0.0, right=0.87)
    # gs_cb.update(hspace=0.5, wspace=0.0, left=0.93, right=0.96)

    # arrays to hold axis handles
    ax = np.empty((nrows, ncols), dtype='object')
    cbar_ax = np.empty((nrows, ncols), dtype='object')
    for this_row in range(nrows):
        for this_col in range(ncols):
            ax[this_row, this_col] = plt.subplot(
                gs[this_row, this_col * 2])
            cbar_ax[this_row, this_col] = plt.subplot(
                gs[this_row, (this_col * 2) + 1])

    # make sure we got a unique axis for all maps, colorbars
    assert np.unique(np.concatenate((ax, cbar_ax))).size == (nrows * ncols * 2)
    fig.savefig(timutils.io.get_temp_filename(prefix='image01_'))
    return(fig, ax, cbar_ax)


def draw_map(map_axis=None, cb_axis=None):
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

    lru_map = na_map.NAMapFigure(map_axis=map_axis,
                                 label_latlon=True,
                                 cb_axis=cb_axis,
                                 t_str=None)
    map_axis.figure.savefig(timutils.io.get_temp_filename(prefix='image02_'))

    cm = lru_map.map.pcolor(stem_lon, stem_lat,
                            maskoceans(stem_lon, stem_lat, LRU),
                            cmap=fcos_cmap,
                            latlon=True,
                            norm=fcos_norm)
    map_axis.figure.savefig(timutils.io.get_temp_filename(prefix='image03_'))
    cbar = map_axis.figure.colorbar(cm, cax=cb_axis, format='%0.2f')
    map_axis.figure.savefig(timutils.io.get_temp_filename(prefix='image04_'))
    # ticks = np.concatenate([np.linspace(1.12, 1.84, 5), np.array([1.61])])
    # cbar = lru_map.fig.colorbar(cm, cax=lru_map.ax_cmap, ticks=ticks)
    # get rid of white lines in colorbar?
    # http://stackoverflow.com/questions/15003353/why-does-my-colorbar-have-lines-in-it
    cbar.solids.set_edgecolor("face")
    cbar.ax.set_title('LRU')
    return(lru_map)


if __name__ == "__main__":
    plt.close('all')
    fig, ax, cb_ax = setup_panel_array(nrows=1, ncols=2)
    c3c4_map = draw_map(ax[0, 0], cb_ax[0, 0])

    c3c4_map = draw_map(ax[0, 1], cb_ax[0, 1])
    fig.savefig('c3c4_map.png')
