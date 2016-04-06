import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import numpy as np
import os
import os.path
from datetime import datetime
from mpl_toolkits.basemap import maskoceans
from timutils import midpt_norm
from stem_pytools import calc_drawdown
from stem_pytools import STEM_parsers as sp
from stem_pytools import NERSC_data_paths as ndp
from stem_pytools import aqout_postprocess as aqp
from stem_pytools import na_map
from stem_pytools.domain import STEM_Domain

# I put the Whelan-Kettle hybrid soil fluxes through STEM in pmol m-2
# s-1 when it was expecting mol m-2 s-1.  So the AQOUT concentrations
# are too high be a factor of 1e12.
mol_pmol_correction = 1e-12


class MapPanel(object):
    """class to create and populate a grid of maps, each with its own
    colorbar
    """
    def __init__(self, nrows, ncols, figsize=(6, 6)):
        """create a matrix of axes for maps with individual color bars

        create a figure containing a matrix of axes with nrows rows
        and ncols columns, where each "panel" consists of a main axis
        using 90% of the horizontal space and a secondary colorbar
        axis occupying 10% of the horizontal space.

        ARGS
        nrows (int): number of rows of axes
        ncols (int): number of columns of axes
        figsize ((float, float)): width and height of figure in inches
        """
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
        self.fig = fig
        self.ax = ax

    def draw_map(self,
                 data,
                 map_axis_idx=None,
                 cb_axis=None, cb_format='%0.2f',
                 label_lat=False, label_lon=False,
                 vmin=None, vmax=None, midpoint=None,
                 bands_above_mdpt=5,
                 bands_below_mdpt=5,
                 cmap=plt.get_cmap('Blues'),
                 panel_lab='a',
                 extend='neither'):

        """
        returns: stem_pytools.na_map.NAMapFigure object containing the map
        """

        d = STEM_Domain()
        stem_lon = d.get_lon()
        stem_lat = d.get_lat()
        this_cmap, this_norm = midpt_norm.get_discrete_midpt_cmap_norm(
            vmin=vmin,
            vmax=vmax,
            midpoint=midpoint,
            bands_above_mdpt=bands_above_mdpt,
            bands_below_mdpt=bands_below_mdpt,
            extend=extend,
            this_cmap=cmap)

        this_ax = self.ax[map_axis_idx]
        map_obj = na_map.NAMapFigure(map_axis=this_ax,
                                     label_latlon=(label_lat, label_lon),
                                     lon_label_interval=30,
                                     cb_axis=cb_axis,
                                     t_str=None)

        cm = map_obj.map.pcolor(stem_lon, stem_lat,
                                maskoceans(stem_lon, stem_lat, data),
                                cmap=this_cmap,
                                latlon=True,
                                norm=this_norm)
        print 'colorbar format: ' + cb_format
        cbar = this_ax.figure.colorbar(cm,
                                       ax=this_ax,
                                       format=cb_format,
                                       orientation='horizontal')
        ticklabs = cbar.ax.get_xticklabels()
        tickvals = np.array([float(x.get_text()) for x in ticklabs])
        cbar.ax.set_xticklabels(tickvals, rotation=-45)
        # place panel label in upper left
        this_ax.text(-0.2, 1.0, panel_lab, transform=this_ax.transAxes)

    def save(self, fname):
        """save the figure to a file

        ARGS:
        fname (string): full path of file to save figure to.  File
           type is determined by the extension of fname. Supported
           file types are png, pdf.
        """
        self.fig.savefig(fname)


def get_anthro_fCOS(run_key):
    """parse July-August 2008 anthropogenic COS surface flux
    """
    Jul1 = datetime(2008, 7, 1)
    Aug31 = datetime(2008, 8, 31, 23, 59, 59)
    flux = sp.parse_STEM_var(nc_fname=ndp.get_runs()[run_key].fcos_path,
                             varname='Coal_COS',
                             t0=Jul1,
                             t1=Aug31)
    fcos_mean = flux['data'].mean(axis=0).squeeze()
    return fcos_mean


def get_COS_concentration(run_key):
    """parse an AQOUT file and calculate daily mean [COS] for 1 July 2008
    to 31 Aug 2008.  Essentially a wrapper function.

    :param run_key: The key for the NERSC_data_paths.get_runs() dict
        referencing the STEM run desired.

    """
    print('parsing AQOUT and calculating daily stats')
    this_run = ndp.get_runs()[run_key]
    aqc = aqp.aqout_container(this_run.aqout_path)
    aqc.parse(t0=datetime(2008, 7, 1), t1=datetime(2008, 8, 31, 23, 59, 59))
    if "Anthro" in run_key:
        print(('multiplying {} [COS] by 1000 '
               'as per email from Andrew').format(run_key))
        aqc.data[0] = aqc.data[0] * 1000
    raw_data = aqc.data[0]
    # aqc.sum()
    # aqc.calc_stats()
    # aqc.stats_to_netcdf(os.path.join(os.getenv('SCRATCH'),
    #                                  'STEM_run_daily.nc'))

    print('calculating drawdown')
    dd = calc_drawdown.calc_STEM_COS_drawdown(
        aqc.data[0]).mean(axis=0).squeeze()
    return(aqc, dd, raw_data)

if __name__ == "__main__":

    plt.close('all')

#     aqc, dd, raw_data = get_COS_concentration('Fsoil_Hybrid5Feb')

#     print('drawing map')
#     cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(vmin=-10,
#                                                          vmax=25,
#                                                          midpoint=0.0,
#                                                          bands_above_mdpt=3,
#                                                          bands_below_mdpt=7)
#     f_or_p = 'pretty'
#     STEM_mapper.Mapper124x124(np.mean(dd[:], axis=0).squeeze()).draw_map(
#         fast_or_pretty=f_or_p,
#         t_str=('1 Jul - 31 Aug 2008 mean STEM drawdown, '
#                'Whelan-Kettle "hybrid" Fsoil'),
#         cmap=cmap,
#         norm=norm)
#     plt.gcf().savefig(os.path.join(os.getenv('HOME'),
#                                    'plots', 'dd_map_Fsoil.png'))

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#     aqc, dd, raw_data = get_COS_concentration('climatological_bnd')
#     print('drawing map')
#     this_data = np.mean(dd[10:, ...].squeeze(), axis=0)
#     cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(vmin=this_data.min(),
#                                                          vmax=this_data.max(),
#                                                          midpoint=0.0,
#                                                          bands_above_mdpt=7,
#                                                          bands_below_mdpt=3)
#     f_or_p = 'pretty'
#     STEM_mapper.Mapper124x124(this_data).draw_map(
#         fast_or_pretty=f_or_p,
#         t_str=('Jul-Aug mean STEM drawdown, '
#                'climatological boundaries'),
#         cmap=cmap,
#         norm=norm)
#     plt.gcf().savefig(os.path.join(os.getenv('HOME'),
#                                    'plots', 'dd_map_clim_bounds.png'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    maps_anthro = MapPanel(nrows=2, ncols=2)

    aqc, dd_anthro_zumkehr, raw_data = get_COS_concentration('Anthro_Andrew')
    aqc, dd_anthro_kettle, raw_data = get_COS_concentration('Anthro_Kettle')

    vmin = np.array([dd_anthro_zumkehr.min(),
                     dd_anthro_kettle.min()]).min()
    vmin = -20
    vmax = np.array([dd_anthro_zumkehr.max(),
                     dd_anthro_kettle.max()]).max()
    cb_midpt = np.mean((vmin, vmax))
    maps_anthro.draw_map(dd_anthro_zumkehr,
                         map_axis_idx=(0, 0),
                         label_lat=True,
                         vmin=vmin, vmax=vmax,
                         midpoint=cb_midpt,
                         extend='min')
    maps_anthro.draw_map(dd_anthro_kettle,
                         map_axis_idx=(0, 1),
                         vmin=vmin, vmax=vmax,
                         midpoint=cb_midpt,
                         extend='min',
                         panel_lab='b')

    fcos_andrew = get_anthro_fCOS('Anthro_Andrew')
    fcos_kettle = get_anthro_fCOS('Anthro_Kettle')

    vmin = np.dstack((fcos_andrew, fcos_kettle)).flatten().min()
    vmax = np.dstack((fcos_andrew, fcos_kettle)).flatten().max()
    vmin = 0.0
    vmax = 1.0
    maps_anthro.draw_map(fcos_andrew,
                         map_axis_idx=(1, 0),
                         vmin=vmin, vmax=vmax,
                         midpoint=0.5,
                         bands_above_mdpt=10,
                         bands_below_mdpt=10,
                         cmap=plt.get_cmap('PuOr'),
                         extend='both',
                         panel_lab='c')
    maps_anthro.draw_map(fcos_kettle,
                         map_axis_idx=(1, 1),
                         vmin=vmin, vmax=vmax,
                         midpoint=0.5,
                         bands_above_mdpt=10,
                         bands_below_mdpt=10,
                         cmap=plt.get_cmap('PuOr'),
                         extend='both',
                         panel_lab='d')

    maps_anthro.save(fname=os.path.join(os.getenv('HOME'),
                                        'plots', 'map_anthro.pdf'))

    plt.close('all')
