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
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)

    return(fig, ax)


def draw_c3c4_pct_map(map_axis=None, cb_axis=None,
                      label_lat=False, label_lon=False):
    nc = netCDF4.Dataset(os.path.join(os.environ['PROJ'],
                                      'Data',
                                      'C4_percentage',
                                      'ISLSCP_C4_1DEG_932_regridded',
                                      'C4_pct_124x124.nc'))
    pct = nc.variables['C4pct'][...].squeeze()
    nc.close()

    d = STEM_Domain()
    stem_lon = d.get_lon()
    stem_lat = d.get_lat()

    pct_cmap, pct_norm = midpt_norm.get_discrete_midpt_cmap_norm(
        vmin=0.0,
        vmax=100.0,
        midpoint=50.0,
        this_cmap=plt.get_cmap('Blues'),
        extend='neither')

    pct_map = na_map.NAMapFigure(map_axis=map_axis,
                                 label_latlon=(label_lat, label_lon),
                                 lon_label_interval=30,
                                 cb_axis=cb_axis,
                                 t_str=None)

    cm = pct_map.map.pcolor(stem_lon, stem_lat,
                            maskoceans(stem_lon, stem_lat, pct),
                            cmap=pct_cmap,
                            latlon=True,
                            norm=pct_norm)
    cbar = map_axis.figure.colorbar(cm,
                                    ax=map_axis,
                                    format='%0.1f',
                                    orientation='horizontal')
    ticklabs = cbar.ax.get_xticklabels()
    tickvals = np.array([float(x.get_text()) for x in ticklabs])
    cbar.ax.set_xticklabels(tickvals, rotation=-45)
    # place 'b' label in upper left
    map_axis.text(-0.1, 1.0, 'b', transform=map_axis.transAxes)
    return pct_map


def draw_map(map_axis=None, cb_axis=None, label_lat=False, label_lon=False):
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
        midpoint=1.61,
        extend='neither')

    lru_map = na_map.NAMapFigure(map_axis=map_axis,
                                 label_latlon=(label_lat, label_lon),
                                 lon_label_interval=30,
                                 cb_axis=cb_axis,
                                 t_str=None)

    cm = lru_map.map.pcolor(stem_lon, stem_lat,
                            maskoceans(stem_lon, stem_lat, LRU),
                            cmap=fcos_cmap,
                            latlon=True,
                            norm=fcos_norm)
    cbar = map_axis.figure.colorbar(cm,
                                    ax=map_axis,
                                    format='%0.2f',
                                    orientation='horizontal')
    ticklabs = cbar.ax.get_xticklabels()
    tickvals = np.array([float(x.get_text()) for x in ticklabs])
    cbar.ax.set_xticklabels(tickvals, rotation=-45)
    # place 'b' label in upper left
    map_axis.text(-0.2, 1.0, 'a', transform=map_axis.transAxes)
    return lru_map


if __name__ == "__main__":
    plt.close('all')
    fig, ax = setup_panel_array(nrows=1, ncols=2)
    c3c4_LRU_map = draw_map(ax[0],
                            label_lon=True, label_lat=True)
    c3c4_pct_map = draw_c3c4_pct_map(ax[1],
                                     label_lon=True, label_lat=False)
    fig.tight_layout()
    fig.savefig('c3c4_map.pdf')
