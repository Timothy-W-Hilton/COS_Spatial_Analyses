import matplotlib
matplotlib.use('AGG')

import os
import os.path
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

    C4_data = os.path.join(os.getenv('HOME'), 'Data', 'STEM',
                           'aq_out_data_C4.cpickle')
    constLRU_data = os.path.join(os.getenv('HOME'), 'Data', 'STEM',
                                 'aq_out_data.cpickle')

    cos_dd, gpp, fCOS = assemble_data(constLRU_data,
                                      get_GPP=False,
                                      get_fCOS=False)
    cos_dd_C4, gpp, fCOS = assemble_data(C4_data,
                                         get_GPP=False,
                                         get_fCOS=False)
    cos_dd.update(cos_dd_C4)
    # place model drawdowns into the data frame
    for k, v in cos_dd.items():
        ocs_dd[k] = cos_dd[k][ocs_dd['stem_x'], ocs_dd['stem_y']]

    return(ocs_dd)


def normalize_drawdown(ocs_dd, norm_site='NHA'):
    """
    Within each drawdown "product", normalize to the NHA value.  NHA
    is chosen because it is the maximum observed drawdown in the NOAA
    observations.
    """
    dd_vars = ['ocs_dd',
               'casa_gfed_161',
               'casa_gfed_C4pctLRU',
               'canibis_C4pctLRU',
               'canibis_161']
    for this_var in dd_vars:
        ocs_dd[this_var] = ocs_dd[this_var] / ocs_dd[this_var][norm_site]

    return(ocs_dd)


def draw_box_plot(df, sites_list):
    sns.set_style('ticks')
    sns.set_context('talk')
    pal = sns.color_palette(['#C2A4CE', '#C2A4CE',
                             '#AADEA1', '#AADEA1',
                             '#F4F5F4'])
    g = sns.factorplot(x="sample_site_code",
                       y="drawdown",
                       hue='variable',
                       data=df[df.sample_site_code.isin(sites_list)],
                       kind="bar",
                       palette=pal,
                       x_order=sites_list,
                       aspect=1.25)
    g.despine(offset=10, trim=True)
    g.set_axis_labels("site", "[OCS] drawdown (ppt)")

    return(g)


if __name__ == "__main__":

    ocs_dd = assemble_bar_plot_data()
    ocs_dd = normalize_drawdown(ocs_dd)

    wet_dry = ['CAR', 'BNE', 'WBI', 'OIL', 'NHA']
    east_coast = ['NHA', 'CMA', 'SCA']
    mid_continent = ['ETL', 'DND', 'LEF', 'WBI', 'BNE', 'SGP', 'TGC']

    ocs_dd_long = pd.melt(ocs_dd.reset_index(),
                          id_vars=['sample_site_code'],
                          value_vars=['ocs_dd',
                                      'casa_gfed_161',
                                      'casa_gfed_C4pctLRU',
                                      'canibis_C4pctLRU',
                                      'canibis_161'],
                          value_name='drawdown')

    g = draw_box_plot(ocs_dd_long, east_coast)
    plt.gcf().savefig('/tmp/barplots/barplots_eastcoast_casa_canibis.pdf')

    g = draw_box_plot(ocs_dd_long, wet_dry)
    plt.gcf().savefig('/tmp/barplots/barplots_wetdry.pdf')

    g = draw_box_plot(ocs_dd_long, mid_continent)
    plt.gcf().savefig('/tmp/barplots/barplots_midcontinent.pdf')
