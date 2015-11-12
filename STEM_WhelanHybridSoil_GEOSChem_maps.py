import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import numpy as np
import numpy.ma as ma
import os
import os.path
from datetime import datetime

from timutils import colormap_nlevs
from timutils import midpt_norm
from stem_pytools import calc_drawdown
from stem_pytools import NERSC_data_paths as ndp
from stem_pytools import aqout_postprocess as aqp
from stem_pytools import STEM_mapper
from stem_pytools import noaa_ocs

# I put the Whelan-Kettle hybrid soil fluxes through STEM in pmol m-2
# s-1 when it was expecting mol m-2 s-1.  So the AQOUT concentrations
# are too high be a factor of 1e12.
mol_pmol_correction = 1e-12


def calc_ABL_FreeTrop_Mean_COS(aqc):
    """calculate the mean [COS] for the free troposphere (defined very
    loosely as >= 4000 meters above ground level (AGL) and the
    atmospheric boundary layer (ABL) (defined very loosely as <= 2000
    meters AGL.

    :param aqc: stem_pytools.aqout_postprocess.aqout_container
        instance: contains the [COS] data to be averaged


    """

    # factor to molecules cm-3 to parts per trillion by volume
    mcm3_2_pptv = 1e12

    [agl, asl] = noaa_ocs.get_STEMZ_height(
        topo_fname=os.path.join(os.getenv('SARIKA_INPUT'), 'TOPO-124x124.nc'),
        wrfheight_fname=os.path.join(os.getenv('SARIKA_INPUT'),
                                     'wrfheight-124x124-2008-2009-22levs.nc'))
    t_steps_per_day = 4
    jul8 = (2008191 - 2008060) * t_steps_per_day
    aug31 = jul8 + (62 * t_steps_per_day)
    agl_mean = np.mean(agl[jul8:aug31, ...], axis=0)
    agl_mean = np.broadcast_arrays(agl_mean, aqc.cos_mean[8:, ...])[0]

    free_trop = ma.masked_where(agl_mean >= 4000, aqc.cos_mean[8:, ...])
    free_trop_mean = mcm3_2_pptv * np.mean(np.mean(free_trop,
                                                   axis=0),
                                           axis=0).squeeze()

    abl = ma.masked_where(agl_mean <= 2000, aqc.cos_mean[8:, ...])
    abl_mean = mcm3_2_pptv * np.mean(np.mean(abl, axis=0), axis=0).squeeze()
    return(free_trop_mean, abl_mean)


def get_COS_concentration(run_key):
    """parse an AQOUT file and calculate daily mean [COS] for 1 July 2008
    to 31 Aug 2008.  Essentially a wrapper function.

    :param run_key: The key for the NERSC_data_paths.get_runs() dict
        referencing the STEM run desired.

    """
    print('parsing AQOUT and calculating daily stats')
    this_run = ndp.get_runs()[run_key]
    aqc = aqp.aqout_container(this_run.aqout_path)
    aqc.parse(t0=datetime(2008, 7, 8), t1=datetime(2008, 8, 31, 23, 59, 59))

    raw_data = aqc.data[0]
    aqc.sum()
    aqc.calc_stats()
    # aqc.stats_to_netcdf(os.path.join(os.getenv('SCRATCH'),
    #                                  'STEM_run_daily.nc'))

    print('calculating drawdown')
    dd = calc_drawdown.calc_STEM_COS_drawdown(aqc.cos_mean)
    return(aqc, dd, raw_data)

if __name__ == "__main__":

    aqc, dd, raw_data = get_COS_concentration('Fsoil_Hybrid5Feb')
    free_trop_mean, abl_mean = calc_ABL_FreeTrop_Mean_COS(aqc)

    print('drawing map')
    cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(vmin=-10,
                                                         vmax=25,
                                                         midpoint=0.0,
                                                         bands_above_mdpt=3,
                                                         bands_below_mdpt=7)
    f_or_p = 'pretty'
    STEM_mapper.Mapper124x124(np.mean(dd[:], axis=0).squeeze()).draw_map(
        fast_or_pretty=f_or_p,
        t_str=('1 Jul - 31 Aug 2008 mean STEM drawdown, '
               'Whelan-Kettle "hybrid" Fsoil'),
        cmap=cmap,
        norm=norm)
    plt.gcf().savefig(os.path.join(os.getenv('HOME'),
                                   'plots', 'dd_map_Fsoil.png'))

    cmap, norm = colormap_nlevs.setup_colormap(
        nlevs=9,
        vmin=430,  # np.dstack((abl_mean, free_trop_mean)).min(),
        vmax=455,  # np.dstack((abl_mean, free_trop_mean)).max(),
        cmap=plt.get_cmap('Blues'))

    STEM_mapper.Mapper124x124(abl_mean.squeeze()).draw_map(
        fast_or_pretty=f_or_p,
        t_str='8 Jul - 31 Aug STEM [COS] <= 2000 m; Whelan-Kettle "hybrid" Fsoil',
        cmap=cmap,
        norm=norm)
    plt.gcf().savefig(os.path.join(os.getenv('HOME'),
                                   'plots', 'cos_Fsoil_ABL.png'))

    STEM_mapper.Mapper124x124(free_trop_mean.squeeze()).draw_map(
        fast_or_pretty=f_or_p,
        t_str='8 Jul - 31 Aug STEM [COS] >= 4000 m; Whelan-Kettle "hybrid" Fsoil',
        cmap=cmap,
        norm=norm)
    plt.gcf().savefig(os.path.join(os.getenv('HOME'),
                                   'plots', 'cos_Fsoil_FreeTrop.png'))

    plt.close('all')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    aqc, dd, raw_data = get_COS_concentration('GEOSChem_bounds')
    free_trop_mean, abl_mean = calc_ABL_FreeTrop_Mean_COS(aqc)

    print('drawing map')
    cmap, norm = midpt_norm.get_discrete_midpt_cmap_norm(vmin=-40,
                                                         vmax=25,
                                                         midpoint=0.0,
                                                         bands_above_mdpt=3,
                                                         bands_below_mdpt=7)
    f_or_p = 'pretty'
    STEM_mapper.Mapper124x124(np.mean(dd[:], axis=0).squeeze()).draw_map(
        fast_or_pretty=f_or_p,
        t_str=('8 Jul - 31 Aug 2008 mean STEM drawdown, '
               'GEOS-Chem Boundaries'),
        cmap=cmap,
        norm=norm)
    plt.gcf().savefig(os.path.join(os.getenv('HOME'),
                                   'plots', 'dd_map_GCbounds.png'))

    cmap, norm = colormap_nlevs.setup_colormap(
        nlevs=9,
        vmin=np.dstack((abl_mean, free_trop_mean)).min(),
        vmax=np.dstack((abl_mean, free_trop_mean)).max(),
        cmap=plt.get_cmap('Blues'))

    STEM_mapper.Mapper124x124(abl_mean.squeeze()).draw_map(
        fast_or_pretty=f_or_p,
        t_str='8 Jul - 31 Aug STEM [COS] <= 2000 m; GEOS-Chem boundaries',
        cmap=cmap,
        norm=norm)
    plt.gcf().savefig(os.path.join(os.getenv('HOME'),
                                   'plots', 'cos_GCbounds_ABL.png'))

    STEM_mapper.Mapper124x124(free_trop_mean.squeeze()).draw_map(
        fast_or_pretty=f_or_p,
        t_str='8 Jul - 31 Aug STEM [COS] >= 4000 m; GEOS-Chem boundaries',
        cmap=cmap,
        norm=norm)
    plt.gcf().savefig(os.path.join(os.getenv('HOME'),
                                   'plots', 'cos_GCbounds_FreeTrop.png'))

    plt.close('all')
