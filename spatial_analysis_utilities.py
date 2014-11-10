import os, os.path
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

#for parsing STEM grid latitude, longitude, and topography
from stem_pytools import STEM_parsers
# for parsing and plotting NOAA [OCS] observations
from stem_pytools import noaa_ocs
# for plotting observations on a map of N America
from stem_pytools import na_map
from map_grid import map_grid_main

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

    ## this code calculates only mean OCS
    # ocs_by_site = data.obs.analysis_value.groupby(data.obs.sample_site_code)
    # ocs_mean = ocs_by_site.aggregate(np.mean)
    # return(ocs_mean)

def preprocess_NOAA_airborne_data_for_JA_spatial_analysis(noaa_dir):

    data = noaa_ocs.get_all_NOAA_airborne_data(noaa_dir)

    topo_file = '/Users/tim/work/Data/STEM/input/TOPO-124x124.nc'
    wrf_height_file = '/Users/tim/work/Data/STEM/input/wrfheight-124x124-22levs.nc'
    stem_lon, stem_lat, topo = STEM_parsers.parse_STEM_coordinates(topo_file)
    data.get_stem_xy(stem_lon, stem_lat)
    data.get_stem_z(topo_fname=topo_file,
                    wrfheight_fname=wrf_height_file)

    # some observations list longitude of -999: remove those
    keep_idx = data.obs['sample_longitude'].values > -998
    # keep observations from July and August only
    keep_idx = keep_idx & np.in1d(data.obs['sample_month'].values, [7,8])
    # remove observations in Alaska - this is outside of the STEM
    # domain.  I can do this crudely by cutting it off at 140 deg W
    # longitude.
    keep_idx = keep_idx & (data.obs.sample_longitude > -140)
    data.obs = data.obs.loc[keep_idx]

    return(data)

def plot_site_altitude_histograms(data, savefig=False):
    """
    plot site-by-site histograms of OCS observation altitudes
    """

    ## plot site-by-site histograms of observation altitudes
    g = sns.FacetGrid(data.obs,
                      row="sample_site_code",
                      size=2,
                      aspect=5)
    g.map(plt.hist,
          "sample_altitude",
          bins=np.linspace(0,15000,200))
    g.set_xlabels('obseration altitude (m)')
    g.set_ylabels('# of obs')  #check those units
    g.set_titles(row_template='site: {row_name}')
    if savefig:
        plt.savefig(os.path.join('/Users', 'tim', 'work', 'Plots',
                                 'SpatialAnalysisPaper',
                                 'altitude_histograms_by_site.pdf'))
    return(g)

def plot_site_mean_drawdown(dd_df, all_data, cmap=None, norm=None):

    agg_vars = ['sample_latitude', 'sample_longitude', 'sample_site_code']
    data_agg = all_data.obs[agg_vars].groupby(['sample_site_code']).aggregate(np.mean)

    dd_df = dd_df.reset_index().groupby('sample_site_code').mean()

    df = pd.merge(dd_df, data_agg, left_index=True, right_index=True)

    dd_map = na_map.NAMapFigure(t_str='mean OCS drawdown')
    dd_map.map.scatter(df.sample_longitude.values,
                       df.sample_latitude.values,
                       c = cmap(norm(df.ocs_dd.values)),
                       s=200,
                       latlon=True)
    return(dd_map, df)

def plot_site_drawdown_timeseries(dd_df):

    dd_df['doy'] = dd.index.get_level_values('date').dayofyear
    dd_df = dd_df.reset_index().dropna()
    dd_df.date = pd.to_datetime(dd_df.date)
    dd_df['jdate'] = np.array([t.to_julian_date() for t in dd_df.date])
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
    noaa_dir = os.path.join(os.getenv('HOME'), 'work', 'Data', 'NOAA_95244993')

    # parse the data and calculate drawdown
    print 'parsing data'
    data = noaa_ocs.get_all_NOAA_airborne_data(noaa_dir)
    print 'filtering for Jul & Aug, etc.'
    ja_data = preprocess_NOAA_airborne_data_for_JA_spatial_analysis(noaa_dir)
    print 'calculating drawdown'
    dd = ja_data.calculate_OCS_daily_vert_drawdown()

    if plot_site_mean_drawdown_switch:
        map_objs, cos_cmap, cos_norm = map_grid_main()
        # plt.scatter(dd_df.sample_longitude.values, dd_df.sample_latitude.values, c=cos_cmap(cos_norm(dd_df.ocs_dd.values)))
        dd_map, dd_df = plot_site_mean_drawdown(dd, data,
                                                cmap=cos_cmap, norm=cos_norm)

    if draw_site_drawdown_timeseries:
        plot_site_drawdown_timeseries(dd)

    if draw_site_locations_map:
        # draw observation sites map
        location_map = data.plot_obs_site_locations()
        location_map.fig.savefig(os.path.join('/Users', 'tim', 'work',
                                              'Plots',
                                              'SpatialAnalysisPaper',
                                              'noaa_obs_sites.pdf'))
        plt.close(location_map.fig)

    if draw_observation_altitude_histograms:
        draw_observation_altitude_histograms(data)

    if False:

        ##pd.set_option('display.width',160)

        ## to move levels of a multiindex to columns, use reset_index
        ## ocs_mean.head().reset_index(level=0)

        ## to select based on values of a multiindex, use e.g.
        ## ocs_mean.index.get_level_values('sample_site_code') == 'OIL'

        lef_idx = ocs_mean.index.get_level_values('sample_site_code') == 'LEF'
        lef_data = ocs_mean[lef_idx]

        ## now make two dataframes: one with WLEF data < 2000m and the
        ## other with WLEF data > 4000m.  Then align the data based on
        ## jday.
        lef_hi = lef_data[lef_data.index.get_level_values('alt_bin') == '(4000, 30000]']
        lef_lo = lef_data[lef_data.index.get_level_values('alt_bin') == '(0, 2000]']


        ## to do:
        ## - map of site codes (which sites are where)
        ## - pre process NOAA data for (1) Jul Aug and (2) inside STEM
        ##   domain.  Maybe not (2).
        ##  -
        
