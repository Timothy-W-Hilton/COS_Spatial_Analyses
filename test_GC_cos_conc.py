import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import sys
import os
import os.path
from stem_pytools import aqout_postprocess as aq
from stem_pytools.calc_drawdown import calc_STEM_COS_drawdown
import gradient_bar_plots as gbp


def get_cos_conc():
    cpickle_fname = os.path.join(os.getenv('SCRATCH'),
                                 '2015-11-16_all_runs.cpickle')
    cos_conc_daily = aq.load_aqout_data(cpickle_fname)
    keys_to_remove = ['casa_gfed_pctm_bnd', 'casa_gfed_KV']
    for k in keys_to_remove:
        if k in cos_conc_daily['cos_mean']:
            del cos_conc_daily['cos_mean'][k]
            del cos_conc_daily['cos_std'][k]
            del cos_conc_daily['t'][k]

    cos_conc_daily['cos_mean'] = gbp.calculate_GCbounds_cos(
        cos_conc_daily['cos_mean'])
    # aggregate daily means to a single July-August mean
    cos_conc = cos_conc_daily['cos_mean']
    return(cos_conc)


def mirror_gradient_bar_plots_data_loading():

    ocs_dd = gbp.assemble_bar_plot_data()

    dd_vars = ['NOAA obs', 'GEOS-Chem boundaries', 'CASA-GFED3, LRU=1.61',
               'MPI, LRU=C3/C4', 'Can-IBIS, LRU=1.61',
               'CASA-GFED3, LRU=1.87', 'Kettle, LRU=C3/C4',
               'Kettle, LRU=1.61',
               'SiB, mechanistic canopy', 'SiB, prescribed canopy',
               'Hybrid Fsoil',
               'CASA-m15, LRU=1.61', 'Kettle Fsoil', 'MPI, LRU=1.61',
               'CASA-GFED3, LRU=C3/C4', 'CASA-GFED3, LRU=1.35',
               'Can-IBIS, LRU=C3/C4', 'CASA-m15, LRU=C3/C4']
    ocs_dd_new = gbp.rename_columns(ocs_dd)
    ocs_dd_new = gbp.normalize_drawdown(ocs_dd_new, vars=dd_vars)


def draw_histogram(data, xlab_str, ylab_str, title_str, nbins=50):
    plt.hist(data.flatten(), bins=nbins)
    fig = plt.gcf()
    plt.xlabel(xlab_str)
    plt.ylabel(ylab_str)
    plt.title(title_str)
    return(fig)


if __name__ == "__main__":
    try:
        data = get_cos_conc()

        d = data['casa_gfed_161'] - data['casa_gfed_161, GC']
        molecules_cm3_2_ppt = 1e12
        fig = draw_histogram(d * molecules_cm3_2_ppt,
                             xlab_str='[COS], ppt',
                             ylab_str='number of grid cells',
                             title_str=r'$\Delta$[COS]: CASA GFED3 - CASA GFED3, GC')
        fig.savefig(os.path.join(os.getenv('HOME'), 'plots', 'GC_diff.png'))
        plt.close(fig)

        d_dd = (calc_STEM_COS_drawdown(data['casa_gfed_161']) -
                calc_STEM_COS_drawdown(data['casa_gfed_161, GC']))

        fig = draw_histogram(d_dd,
                             xlab_str='COS vertical drawdown (ppt)',
                             ylab_str='number of grid cells',
                             title_str=r'$\Delta$drawdown: CASA GFED3 - CASA GFED3, GC')
        fig.savefig(os.path.join(os.getenv('HOME'),
                                 'plots', 'GC_dd_diff.png'))
        plt.close(fig)
    finally:
        sys.stdout.write('closing all figures')
        sys.stdout.flush
        plt.close('all')
