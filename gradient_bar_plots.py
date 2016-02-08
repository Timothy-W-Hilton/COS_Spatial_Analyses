import matplotlib
matplotlib.use('AGG')

import os
import os.path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import spatial_analysis_utilities as sau
from stem_pytools import noaa_ocs
from stem_pytools import domain
from stem_pytools import aqout_postprocess as aq
from stem_pytools import calc_drawdown
# from timutils.mpl_fig_joiner import FigJoiner
import map_grid
import draw_c3c4LRU_map


def get_STEM_cos_conc(cpickle_fname=None, const_bounds_cos=4.5e-10):
    """get dict of Jul-Aug mean COS drawdowns for each STEM gridcell.
    Adjust GEOS-Chem boundaries concentration to reflect the *change*
    in [COS] for the dynamic boundaries relative to the constant
    boundaries.

    :param cpickle_fname: path to a cpickle file of [COS] mean,
        standard deviation ,time stamps calculated from AQOUT files.
        That cpickle file will typically be the output of
        stem_pytools.aqout_postprocess.assemble_data.
    :param const_bounds_cos: the COS concentration of the constant
        boundaries; this value will be subtracted out of the GEOS-Chem
        boundaries [COS]
    """
    cos_conc_daily = aq.load_aqout_data(cpickle_fname)
    keys_to_remove = ['casa_gfed_pctm_bnd', 'casa_gfed_KV']
    for k in keys_to_remove:
        if k in cos_conc_daily['cos_mean']:
            del cos_conc_daily['cos_mean'][k]
            del cos_conc_daily['cos_std'][k]
            del cos_conc_daily['t'][k]

    print(('multiplying Anthro [COS] by 1000 '
           'as per email from Andrew'))
    cos_conc_daily['cos_mean']['Anthro_Kettle'] *= 1e3
    cos_conc_daily['cos_mean']['Anthro_Andrew'] *= 1e3

    cos_conc_daily['cos_mean']['CASA-GFED3, Kettle Anthropogenic'] = cos_conc_daily['cos_mean']['casa_gfed_161'] + cos_conc_daily['cos_mean']['Anthro_Kettle']
    cos_conc_daily['cos_mean']['CASA-GFED3, Zumkehr Anthropogenic'] = cos_conc_daily['cos_mean']['casa_gfed_161'] + cos_conc_daily['cos_mean']['Anthro_Andrew']

    # calculate CASA-GFED plus soil fluxes, anthro fluxes
    cos_conc_daily['cos_mean']['CASA-GFED3, Kettle Anthropogenic'] = (cos_conc_daily['cos_mean']['casa_gfed_161'] +
                                                    cos_conc_daily['cos_mean']['Anthro_Kettle'])
    cos_conc_daily['cos_mean']['CASA-GFED3, Zumkehr Anthropogenic'] = (cos_conc_daily['cos_mean']['casa_gfed_161'] +
                                                     cos_conc_daily['cos_mean']['Anthro_Andrew'])
    cos_conc_daily['cos_mean']['CASA-GFED3, Kettle Fsoil'] = (cos_conc_daily['cos_mean']['casa_gfed_161'] +
                                            cos_conc_daily['cos_mean']['Fsoil_Kettle'])
    cos_conc_daily['cos_mean']['CASA-GFED3, Hybrid Fsoil'] = (cos_conc_daily['cos_mean']['casa_gfed_161'] +
                                            cos_conc_daily['cos_mean']['Fsoil_Hybrid5Feb'])

    cos_conc_daily['cos_mean'] = calculate_GCbounds_cos(
        cos_conc_daily['cos_mean'], const_bounds_cos)

    # don't need standard devation for this analysis
    cos_conc = cos_conc_daily['cos_mean']
    # aggregate daily means to a single July-August mean
    cos_conc.update((k, calc_drawdown.calc_STEM_COS_drawdown(v)) for
                    k, v in cos_conc.items())
    cos_conc.update((k, map_grid.daily_to_JulAug(v))
                    for k, v in cos_conc.items())

    return(cos_conc)


def assemble_bar_plot_data(
        cpickle_fname=os.path.join(os.getenv('HOME'),
                                   'STEM_all_runs.cpickle')):
    noaa_dir = sau.get_noaa_COS_data_path()
    noaa_ocs_dd, ocs_daily = sau.get_JA_site_mean_drawdown(noaa_dir)

    d = domain.STEM_Domain()
    stem_lon = d.get_lon()
    stem_lat = d.get_lat()

    (noaa_ocs_dd['stem_x'],
     noaa_ocs_dd['stem_y']) = domain.find_nearest_stem_xy(
        noaa_ocs_dd.sample_longitude,
        noaa_ocs_dd.sample_latitude,
        stem_lon,
        stem_lat)

    stem_ocs_dd = get_STEM_cos_conc(cpickle_fname)

    # place model drawdowns into the data frame
    for k, v in stem_ocs_dd.items():
        noaa_ocs_dd[k] = stem_ocs_dd[k][noaa_ocs_dd['stem_x'],
                                        noaa_ocs_dd['stem_y']]

    return(noaa_ocs_dd)


def calculate_GCbounds_cos(stem_ocs_dd, const_bounds=4.5e-10, verbose=False):
    """Calculate drawdown enhancement or reduction because of dynamic
     boundaries relative to static boundary conditions.  Subtract out
     the 450 pptv static boundary condition from the dynamic
     boundaries STEM [COS] to isolate the impact of the dynamic
     boundaries vs static boundaries.  Static boundary of 450 pptv
     coupled with no surface COS flux must result in [COS] = 450 pptv
     at all places, times.

    ARGS:
    ocs_conc (dict): dict containing `array-lke
        <http://docs.scipy.org/doc/numpy/user/basics.creation.html#converting-python-array-like-objects-to-numpy-arrays>`_
        [COS] values (molecules m-3)
    const_bounds (scalar): [COS] for the constant boundary conditions
        (molecules m-3)
    verbose ({False}|True): if True, display message to stdout for
        each stem_ocs_dd field adjusted
    """
    do_not_adjust = ['sample_latitude', 'sample_longitude', 'analysis_value',
                     'ocs_dd', 'stem_x', 'stem_y', 'climatological_bnd']

    bounds_runs = {}
    for k in stem_ocs_dd.keys():
        if k not in do_not_adjust:
            key_GC = '{}{}'.format(k, ', clim')
            data_GC = (stem_ocs_dd[k][0:-1] - const_bounds +
                       stem_ocs_dd['climatological_bnd'])
            if verbose:
                print 'adding {} to dict'.format(key_GC)
            bounds_runs.update({key_GC: data_GC})
    stem_ocs_dd.update(bounds_runs)
    return(stem_ocs_dd)


def normalize_drawdown(ocs_dd,
                       norm_site='NHA',
                       vars=['ocs_dd',
                             'casa_gfed_161',
                             'casa_gfed_C4pctLRU',
                             'canibis_C4pctLRU',
                             'canibis_161']):
    """
    Within each drawdown "product", normalize to the NHA value.  NHA
    is chosen because it is the maximum observed drawdown in the NOAA
    observations.
    """
    df = ocs_dd.copy()
    # norm_factor = df['NOAA obs'].values.copy()
    # warnings.warn(('normalization factor is set to NOAA_obs; the'
    #                'factor must be *copied* so that when the normalization '
    #                'site is itself normalized subsequent normalizations '
    #                'still work.'))
    for this_var in vars:
        print('Normalizing {}'.format(this_var))
        df[this_var] = df[this_var] / ocs_dd[this_var][norm_site]
        # df[this_var] = df[this_var] / norm_factor

    return(df)


def get_line_styles(df):
    """assign linestyles for gradient plots.  For now, use dashed
    lines for climatological boundaries and solid lines for everything
    else.
    """
    n_vars = df.variable.unique().size
    my_linestyles = list(np.tile(['-'], n_vars))
    is_clim = np.where(["clim" in x for x in df.variable.unique()])[0]
    for idx in is_clim:
        my_linestyles[idx] = '--'
    return my_linestyles


def draw_box_plot(df, sites_list):
    sns.set_style('ticks')
    sns.set_context('paper')

    g = sns.factorplot(x="sample_site_code",
                       y="drawdown",
                       hue='variable',
                       data=df[df.sample_site_code.isin(sites_list)],
                       kind="point",
                       palette=sns.color_palette(
                           "cubehelix",
                           len(df.variable.unique())),
                       x_order=sites_list,
                       aspect=1.25,
                       linestyles=get_line_styles(df))
    # make the left and top axes only extend across the part of the
    # plot where the data are.  That is, make the axis look like "| _"
    # and not "L"
    g.despine(offset=10, trim=True)
    g.set_axis_labels("site", "[OCS] drawdown, normalized to NHA")

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
                    'casa_m15_C4pctLRU': 'CASA-m15, LRU=C3/C4',
                    'Fsoil_Kettle': 'Kettle Fsoil',
                    'Fsoil_Hybrid5Feb': 'Hybrid Fsoil',
                    'GEOSChem_bounds': 'GEOS-Chem boundaries',
                    'climatological_bnd': 'climatological boundaries',
                    'SiB_mech': 'SiB, mechanistic canopy',
                    'SiB_calc': 'SiB, prescribed canopy',
                    'Anthro_Kettle': 'Anthropogenic, Kettle',
                    'Anthro_Andrew': 'Anthropogenic, Zumkehr'}
    df_out = df.copy()
    for this_col in df_out.columns.values:
        if this_col in columns_dict.keys():
            df_out.rename(columns=lambda x: x.replace(this_col,
                                                      columns_dict[this_col]),
                          inplace=True)
            print "replaced {} with {}".format(this_col,
                                               columns_dict[this_col])
    return(df_out)


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


def plot_all_gradients(ocs_dd, plot_vars, fname_suffix):

    ocs_dd_long = pd.melt(ocs_dd.reset_index(),
                          id_vars=['sample_site_code'],
                          value_vars=plot_vars,
                          value_name='drawdown')

    figs = []
    g = draw_box_plot(ocs_dd_long, gradients['east_coast'])
    g.ax.set_title('East Coast N -- S (climatological column mean bounds)')
    plt.gcf().savefig(
        os.path.join(os.getenv('HOME'),
                     'plots',
                     'ECoast_model_components_{}_NHAnorm.svg'.format(fname_suffix)))

    return g
    # figs.append(plt.gcf())

    # g = draw_box_plot(ocs_dd_long, gradients['wet_dry'])
    # g.ax.set_title('West -- East (Dry -- Wet)')
    # figs.append(plt.gcf())

    # g = draw_box_plot(ocs_dd_long, gradients['mid_continent'])
    # g.ax.set_title('Midcontinent N -- S')
    # figs.append(plt.gcf())

    # # save figs to single svg file
    # fj = FigJoiner(
    #     figs,
    #     os.path.join(os.getenv('HOME'),
    #                  'plots',
    #                  'model_components_{}_NHAnorm.svg'.format(fname_suffix)))
    # fj.join()
    # fj.close_figs()

if __name__ == "__main__":

    try:
        gradients = {'wet_dry': ['CAR', 'BNE', 'WBI', 'OIL', 'NHA'],
                     'east_coast': ['NHA', 'CMA', 'SCA'],
                     'mid_continent': ['ETL', 'DND', 'LEF', 'WBI',
                                       'BNE', 'SGP', 'TGC']}

        ocs_dd = assemble_bar_plot_data()

        ocs_dd_renamed = rename_columns(ocs_dd)
        dd_vars = ['NOAA obs', 'GEOS-Chem boundaries',
                   'climatological boundaries', 'CASA-GFED3, LRU=1.61',
                   'MPI, LRU=C3/C4', 'Can-IBIS, LRU=1.61',
                   'CASA-GFED3, LRU=1.87', 'Kettle, LRU=C3/C4',
                   'Kettle, LRU=1.61',
                   'SiB, mechanistic canopy', 'SiB, prescribed canopy',
                   'Hybrid Fsoil',
                   'CASA-m15, LRU=1.61', 'Kettle Fsoil', 'MPI, LRU=1.61',
                   'CASA-GFED3, LRU=C3/C4', 'CASA-GFED3, LRU=1.35',
                   'Can-IBIS, LRU=C3/C4', 'CASA-m15, LRU=C3/C4',
                   'CASA-GFED3, Kettle Anthropogenic',
                   'CASA-GFED3, Zumkehr Anthropogenic',
                   'CASA-GFED3, Kettle Fsoil',
                   'CASA-GFED3, Hybrid Fsoil']
        dd_vars_bounds = [''.join([k, ', clim']) for k in dd_vars[2:]
                          if k is not 'climatological boundaries']
        dd_vars = dd_vars + dd_vars_bounds

        ocs_dd_norm = normalize_drawdown(ocs_dd_renamed, vars=dd_vars)

        vars = ['NOAA obs',
                'SiB, prescribed canopy',
                'SiB, mechanistic canopy',
                'CASA-GFED3, LRU=1.61',
                'CASA-GFED3, LRU=C3/C4',
                'Kettle, LRU=1.61',
                'Can-IBIS, LRU=1.61',
                'CASA-GFED3, Kettle Anthropogenic',
                'CASA-GFED3, Zumkehr Anthropogenic',
                'CASA-GFED3, Kettle Fsoil',
                'CASA-GFED3, Hybrid Fsoil',
                'CASA-GFED3, LRU=1.61, clim',
                'SiB, mechanistic canopy, clim',
                'Can-IBIS, LRU=1.61, clim']
        g = plot_all_gradients(ocs_dd_norm, vars, '02Feb')

        # # show east coast sites
        # ocs_dd_renamed.ix[['NHA', 'SCA', 'CMA']][['analysis_value', 'NOAA obs']]
        # # show east coast mean
        # ocs_dd_renamed.ix[['NHA', 'SCA', 'CMA']][['analysis_value',
        #                                           'NOAA obs']].mean()
        # gradient_map = draw_gradient_map(gradients)
        # gradient_map.fig.savefig(os.path.join(tmpdir, 'gradients_map.pdf'))
    finally:
        plt.close('all')
