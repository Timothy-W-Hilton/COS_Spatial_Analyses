import matplotlib
matplotlib.use('AGG')

import os, os.path
import socket
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import spatial_analysis_utilities as sau
from stem_pytools import noaa_ocs
from stem_pytools import STEM_parsers
from map_grid import assemble_data

def assemble_bar_plot_data():
    noaa_dir = sau.get_noaa_COS_data_path()
    ocs_dd, ocs_daily = sau.get_JA_site_mean_drawdown(noaa_dir)

    stem_input_dir = os.getenv('SARIKA_INPUT')
    topo_file = os.path.join(stem_input_dir, 'TOPO-124x124.nc')
    stem_lon, stem_lat, topo = STEM_parsers.parse_STEM_coordinates(topo_file)

    ocs_dd['stem_x'], ocs_dd['stem_y'] = noaa_ocs.find_nearest_stem_xy(
        ocs_dd.sample_longitude, 
        ocs_dd.sample_latitude,
        stem_lon,
        stem_lat)

    cos_dd, gpp, fCOS = map_grid = assemble_data(sau.get_aqout_data_path(),
                                                 get_GPP=False,
                                                 get_fCOS=False)

    # place model drawdowns into the data frame
    for k, v in cos_dd.items():
        ocs_dd[k] = cos_dd[k][ocs_dd['stem_x'], ocs_dd['stem_y']]

    return(ocs_dd)

def draw_box_plot(df, sites_list):
    g = sns.factorplot(x="sample_site_code", 
                       y="drawdown", 
                       hue='variable', 
                       data=df[df.sample_site_code.isin(sites_list)],
                       kind="bar", 
                       palette="PRGn", 
                       x_order=sites_list,
                       aspect=1.25)
    g.despine(offset=10, trim=True)
    g.set_axis_labels("site", "[OCS] drawdown (ppt)")
    
    return(g)

if __name__ == "__main__":

    ocs_dd = assemble_bar_plot_data()
    
    wet_dry = ['CAR', 'BNE', 'WBI', 'OIL', 'NHA']
    east_coast = ['NHA', 'CMA', 'SCA']
    
    ocs_dd_long = pd.melt(ocs_dd.reset_index(),
                          id_vars=['sample_site_code'],
                          value_vars=['ocs_dd', 'casa_gfed_161',
                                      'canibis_161', 'casa_gfed_187',
                                      'kettle_161', 'casa_m15_161', 'MPI_161',
                                      'casa_gfed_135'], 
                          value_name='drawdown')

    g = draw_box_plot(ocs_dd_long, east_coast)
    g = draw_box_plot(ocs_dd_long, wet_dry)
    plt.gcf().savefig('/tmp/barplots.pdf')
    

