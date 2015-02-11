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
import draw_c3c4LRU_map

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
    g = sns.factorplot(x="sample_site_code",
                       y="drawdown",
                       hue='variable',
                       data=df[df.sample_site_code.isin(sites_list)],
                       kind="bar",
                       palette=sns.color_palette("Paired", 5),
                       x_order=sites_list,
                       aspect=1.25)
    g.despine(offset=10, trim=True)
    g.set_axis_labels("site", "[OCS] drawdown (ppt)")

    return(g)


def rename_columns(df):
    """rename drawdown columns into more human-readable strings.  The
    motivation for this is that these column names eventually become
    the labels in the barplot legends, and changing them here seemed
    easier (though less elegant, perhaps) than digging through the
    seaborn facetgrid object to access and change the legend labels
    and then redrawing the plot.

    """
    columns_dict = {'ocs_dd': 'NOAA obs',
                    'casa_gfed_161': 'CASA-GFED3, LRU=1.61',
                    'MPI_C4pctLRU': 'MPI, LRU=C3/C4',
                    'canibis_161': 'Can-IBIS, LRU=1.61',
                    'casa_gfed_187': 'CASA-GFED3, LRU=1.87',
                    'kettle_C4pctLRU': 'Kettle, LRU=C3/C4',
                    'kettle_161': 'Kettle, LRU=1.61',
                    'casa_m15_161': 'CASA-m15, LRU=1.61',
                    'MPI_161': 'MPI, LRU=1.61',
                    'casa_gfed_C4pctLRU': 'CASA-GFED3, LRU=C3/C4',
                    'casa_gfed_135': 'CASA-GFED3, LRU=1.35',
                    'canibis_C4pctLRU': 'Can-IBIS, LRU=C3/C4',
                    'casa_m15_C4pctLRU': 'CASA-m15, LRU=C3/C4'}
    df = df.rename(columns=columns_dict)
    return(df)


def draw_gradient_map(gradient_sites_dict):
    """draw a NOAA observation site "spatial gradient" over the C3/C4 LRU
    map
    """
    lru_map = draw_c3c4LRU_map.draw_map()
    site_coords = noaa_ocs.get_all_NOAA_airborne_data(
        sau.get_noaa_COS_data_path())
    site_coords = site_coords.get_sites_lats_lons()

    markers = ['x', 'o', 's']
    count = 0
    for grad_name, grad_sites in gradient_sites_dict.items():
        this_gradient = site_coords.loc[grad_sites]
        print(this_gradient)

        black = '#000000'
        # #1b9e77   turquoise color from colorbrewer2 Dark2 palette
        print('drawing: {}'.format(grad_name))
        lru_map.map.plot(this_gradient['sample_longitude'].values,
                         this_gradient['sample_latitude'].values,
                         color=black,
                         linewidth=2.0,
                         latlon=True)
        lru_map.map.scatter(this_gradient['sample_longitude'].values,
                            this_gradient['sample_latitude'].values,
                            latlon=True,
                            color=black,
                            marker=markers[count],
                            linewidth=3,
                            s=160,
                            facecolors='None')
        count = count + 1
    return(lru_map)


if __name__ == "__main__":

    gradients = {'wet_dry': ['CAR', 'BNE', 'WBI', 'OIL', 'NHA'],
                 'east_coast': ['NHA', 'CMA', 'SCA'],
                 'mid_continent': ['ETL', 'DND', 'LEF', 'WBI',
                                   'BNE', 'SGP', 'TGC']}

    ocs_dd = assemble_bar_plot_data()
    ocs_dd = normalize_drawdown(ocs_dd)

    ocs_dd_new = rename_columns(ocs_dd)
    ocs_dd_long = pd.melt(ocs_dd_new.reset_index(),
                          id_vars=['sample_site_code'],
                          value_vars=['NOAA obs',
                                      'CASA-GFED3, LRU=1.61',
                                      'CASA-GFED3, LRU=C3/C4',
                                      'Can-IBIS, LRU=1.61',
                                      'Can-IBIS, LRU=C3/C4'],
                          value_name='drawdown')

    g = draw_box_plot(ocs_dd_long, gradients['east_coast'])
    plt.gcf().savefig('/tmp/barplots/barplots_eastcoast_casa_canibis.pdf')

    g = draw_box_plot(ocs_dd_long, gradients['wet_dry'])
    plt.gcf().savefig('/tmp/barplots/barplots_wetdry.pdf')

    g = draw_box_plot(ocs_dd_long, gradients['mid_continent'])
    plt.gcf().savefig('/tmp/barplots/barplots_midcontinent.pdf')

    gradient_map = draw_gradient_map(gradients)
    gradient_map.fig.savefig('/tmp/gradients_map.pdf')
