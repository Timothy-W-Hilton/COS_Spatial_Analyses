import matplotlib
matplotlib.use('AGG')
# matplotlib.rcParams.update({'font.size': 20})
font = {'family': 'Bitstream Vera Sans',
        'weight': 'normal',
        'size': 20}
matplotlib.rc('font', **font)

import os
import os.path
import socket
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# for parsing STEM grid latitude, longitude, and topography
from stem_pytools import STEM_parsers
# for parsing and plotting NOAA [OCS] observations
from stem_pytools import noaa_ocs
# for finding files
from stem_pytools import ecampbell300_data_paths as edp
# for plotting observations on a map of N America
from stem_pytools import na_map
from map_grid import map_grid_main


def get_aqout_data_path():
    """return the full path to the directory containing (pre-parsed)
    aqout data files on either ecampbell300 or Tim's laptop

    """
    if 'Timothys-MacBook-Air.local' in socket.gethostname():
        aqout_data = (os.path.join(os.getenv('HOME'), 'work', 'Data',
                                   'STEM', 'aq_out_data.cpickle'))
    else:
        aqout_data = os.path.join(os.getenv('HOME'), 'Data', 'STEM',
                                  'aq_out_data_BASC.cpickle')
    return(aqout_data)


def get_noaa_COS_data_path():
    """
    return the full path to the directory containing the NOAA OCS data
    files on either ecampbell300 or Tim's laptop

    """
    if 'Timothys-MacBook-Air.local' in socket.gethostname():
        noaa_dir = os.path.join(os.getenv('HOME'), 'work',
                                'Data', 'NOAA_95244993')
    else:
        noaa_dir = os.path.join('/', 'project', 'projectdirs',
                                'm2319',
                                'Data',
                                'NOAA_95244993')
    return(noaa_dir)


def get_site_mean_cos(data):
    """
    compute mean cos concentration at each site in a NOAA_OCS dataset.

    INPUTS
      data: noaa_ocs.NOAA_OCS object

    OUTPUTS
      ocs_mean: noaa_ocs.NOAA_OCS object containing the mean [OCS]
        value
    """

    data_by_site = data.obs.groupby('sample_site_code')
    ocs_mean = data_by_site.aggregate(np.mean)
    return(ocs_mean)

    # this code calculates only mean OCS
    # ocs_by_site = data.obs.analysis_value.groupby(data.obs.sample_site_code)
    # ocs_mean = ocs_by_site.aggregate(np.mean)
    # return(ocs_mean)


def preprocess_NOAA_airborne_data_for_JA_spatial_analysis(noaa_dir):
    """
    (1) parses all NOAA COS observation files from noaa_dir
    (2) assigns each observation to a STEM x, y, and z cell
    (3) removes all observations with a [COS] value of -999
    (4) removes all observations that are not in July or August
    (5) removes all observations west of 140 deg W longitude (this
        roughly correspondes to the western boundary of the 124x124
        STEM domain)

    INPUTS
    noaa_dir: full path to a directory containing NOAA COS observation files

    OUTPUTS
    a noaa_ocs object

    SEE ALSO
    stem_pytools.noaa_ocs
    """

    data = noaa_ocs.get_all_NOAA_airborne_data(noaa_dir)

    stem_input_dir = os.getenv('SARIKA_INPUT')
    topo_file = os.path.join(stem_input_dir, 'TOPO-124x124.nc')
    wrf_height_file = os.path.join(stem_input_dir,
                                   'wrfheight-124x124-22levs.nc')
    stem_lon, stem_lat, topo = STEM_parsers.parse_STEM_coordinates(topo_file)
    data.get_stem_xy(stem_lon, stem_lat)
    data.get_stem_z(topo_fname=topo_file,
                    wrfheight_fname=wrf_height_file)

    # some observations list longitude of -999: remove those
    keep_idx = data.obs['sample_longitude'].values > -998
    # keep observations from July and August only
    keep_idx = keep_idx & np.in1d(data.obs['sample_month'].values, [7, 8])
    # remove observations in Alaska - this is outside of the STEM
    # domain.  I can do this crudely by cutting it off at 140 deg W
    # longitude.
    keep_idx = keep_idx & (data.obs.sample_longitude > -140)
    data.obs = data.obs.loc[keep_idx]

    return(data)


def get_JA_site_mean_drawdown(noaa_dir):
    """
    calculate July-August [COS] drawdown, both daily and for the full
    62-day period (1 July through 31 Aug)
    """
    # parse the data and calculate drawdown
    print 'filtering for Jul & Aug, etc.'
    ja_data = preprocess_NOAA_airborne_data_for_JA_spatial_analysis(noaa_dir)
    print 'calculating drawdown'
    ja_daily_dd = ja_data.calculate_OCS_daily_vert_drawdown()

    agg_vars = ['sample_latitude', 'sample_longitude',
                'sample_site_code', 'analysis_value']
    ja_mean_ocs = ja_data.obs[agg_vars].groupby(
        ['sample_site_code']).aggregate(np.mean)

    ja_mean_dd = ja_daily_dd.reset_index().groupby(
        ['sample_site_code']).aggregate(np.mean)

    ja_mean = ja_mean_ocs.join(ja_mean_dd)
    return(ja_mean, ja_daily_dd)


def plot_site_altitude_histograms(data, savefig=False):
    """
    plot site-by-site histograms of OCS observation altitudes
    """
    try:
        import seaborn as sns
    except ImportError:
        print('Seaborn not available on this system; exiting')
        return(None)
    else:
        # plot site-by-site histograms of observation altitudes
        g = sns.FacetGrid(data.obs,
                          row="sample_site_code",
                          size=2,
                          aspect=5)
        g.map(plt.hist,
              "sample_altitude",
              bins=np.linspace(0, 15000, 200))
        g.set_xlabels('obseration altitude (m)')
        g.set_ylabels('# of obs')  # check those units
        g.set_titles(row_template='site: {row_name}')
        if savefig:
            plt.savefig(os.path.join('/Users', 'tim', 'work', 'Plots',
                                     'SpatialAnalysisPaper',
                                     'altitude_histograms_by_site.pdf'))
        return(g)


def plot_site_mean_drawdown(dd, cmap=None, norm=None, dd_map=None):

    print('warning - [COS] drawdown < 0.0 reset to 0.0')
    dd.ocs_dd[dd.ocs_dd < 0.0] = 0.0
    print('warning - [COS] drawdown NaNs replaced by -1')
    dd.ocs_dd = dd.ocs_dd.fillna(-1)

    if dd_map is None:
        dd_map = na_map.NAMapFigure(t_str='mean OCS drawdown')

    x, y = dd_map.map(dd.sample_longitude.values,
                      dd.sample_latitude.values)
    dd_map.map.scatter(x, y,
                       c=cmap(norm(dd.ocs_dd.values)),
                       edgecolor='blue',
                       linewidths=1,
                       s=70)
    return(dd_map)


def plot_site_drawdown_timeseries(dd_df):

    try:
        import seaborn as sns
    except ImportError:
        print('Seaborn not available on this system; exiting')
        return(None)
    else:
        dd_df['doy'] = dd_df.index.get_level_values('date').dayofyear
        dd_df = dd_df.reset_index().dropna()
        dd_df.date = pd.to_datetime(dd_df.date)
        dd_df['jdate'] = np.array([t.to_julian_date() for t in dd_df.date])
        # dd_df['dt'] = [datetime.datetime.strptime(np.datetime_as_string(t),
        #                                         '%Y-%m-%dT%H:%M:%S.000000000Z')
        #                for t in dd_df['date'].values]
        g = sns.FacetGrid(dd_df,
                          hue="sample_site_code")
        g.map(plt.scatter,
              'doy',
              'ocs_dd')
        plt.xlabel('day of year')
        plt.ylabel('OCS drawdown (ppt)')
        plt.legend()
        return(g)


if __name__ == "__main__":

    # configuration stuff
    draw_site_locations_map = False
    draw_observation_altitude_histograms = False
    draw_site_drawdown_timeseries = False
    plot_site_mean_drawdown_switch = True

    if plot_site_mean_drawdown_switch:
        ocs_dd, ocs_daily = get_JA_site_mean_drawdown(get_noaa_COS_data_path())
        # c4runs = edp.get_C3C4runs()
        basc_runs = edp.get_BASC_runs()
        fig, map_objs, cos_cmap, cos_norm = map_grid_main(
            aqout_data=os.path.join(os.getenv('HOME'), 'Data', 'STEM',
                                    'aq_out_data_BASC.cpickle'),
            models=[k for k in basc_runs.keys()],
            models_str=[v.model for v in basc_runs.values()])
            # models = ['canibis_161', 'casa_gfed_135'],
            # models_str= ['Can-IBIS', 'CASA-GFED3'])

        for i in range(map_objs.shape[1]):
            dd_map = plot_site_mean_drawdown(ocs_dd,
                                             cmap=cos_cmap,
                                             norm=cos_norm,
                                             dd_map=map_objs[3, i])

        fname = '/tmp/maps_basc.pdf'
        print("saving {}".format(fname))
        fig.savefig(fname)

    if draw_site_drawdown_timeseries:
        # FIX THIS: where does dd come from?  -TWH
        plot_site_drawdown_timeseries(dd)

    if draw_site_locations_map:
        # draw observation sites map
        data = noaa_ocs.get_all_NOAA_airborne_data(get_noaa_COS_data_path())
        location_map = data.plot_obs_site_locations()
        location_map.fig.savefig(os.path.join(os.getenv('PLOTS'),
                                              'SpatialAnalysisPaper',
                                              'noaa_obs_sites.pdf'))
        plt.close(location_map.fig)

    if draw_observation_altitude_histograms:
        # FIX THIS: where does data come from?  -TWH
        draw_observation_altitude_histograms(data)

    if False:

        ## pd.set_option('display.width',160)

        ## to move levels of a multiindex to columns, use reset_index
        ## ocs_mean.head().reset_index(level=0)

        ## to select based on values of a multiindex, use e.g.
        ## ocs_mean.index.get_level_values('sample_site_code') == 'OIL'

        lef_idx = ocs_mean.index.get_level_values('sample_site_code') == 'LEF'
        lef_data = ocs_mean[lef_idx]

        ## now make two dataframes: one with WLEF data < 2000m and the
        ## other with WLEF data > 4000m.  Then align the data based on
        ## jday.
        lef_hi = lef_data[
            lef_data.index.get_level_values('alt_bin') == '(4000, 30000]']
        lef_lo = lef_data[
            lef_data.index.get_level_values('alt_bin') == '(0, 2000]']


        ## to do:
        ## - map of site codes (which sites are where)
        ## - pre process NOAA data for (1) Jul Aug and (2) inside STEM
        ##   domain.  Maybe not (2).
        ##  -
