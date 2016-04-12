import matplotlib
matplotlib.use('AGG')
import pandas as pd
from datetime import datetime
import os.path
import numpy as np
import gc

from stem_pytools import NERSC_data_paths as ndp
from aqout_postprocess_spatial_paper import AqoutContainerSpatialPaper
from stem_pytools import noaa_ocs
from stem_pytools.calc_drawdown import calc_STEM_COS_drawdown
import itertools


def pull_site_xy(site_requested, noaa_dir=None):
    """return STEM grid x and y indices for sitecode"""

    if noaa_dir is None:
        noaa_dir = os.path.join('/', 'project',
                                'projectdirs',
                                'm2319',
                                'Data',
                                'NOAA_95244993')
    sites = noaa_ocs.get_sites_summary(noaa_dir, stemxy=True)
    sites.index = sites.site_code
    return (sites.loc[site_requested, 'x_stem'],
            sites.loc[site_requested, 'y_stem'])


stem_input_dir = os.path.join('/', 'project',
                              'projectdirs',
                              'm2319',
                              'Data',
                              'STEM_124x124_NA_inputs')
topo_file = os.path.join(stem_input_dir,
                         'TOPO-124x124.nc')
wrf_file = os.path.join(stem_input_dir,
                        'wrfheight-124x124-2008-2009-22levs.nc')
x_site, y_site = pull_site_xy('WBI')


t0 = datetime.now()

runs = ndp.get_Spatial_Paper_runs()
Fplant = ['SiB_calc', 'SiB_mech', 'casa_gfed_161',
          'casa_gfed_C4pctLRU', 'canibis_161', 'canibis_C4pctLRU']
Fsoil = ['Fsoil_Kettle', 'Fsoil_Hybrid5Feb']
Fanthro = ['Anthro_Andrew',  'Anthro_Kettle']
Fbounds = [None, 'climatological_bnd']

# gather all permutations and remove None (None is placeholder for
# constant boundaries, which has no explicit AQOUT file)
all_combos = [c for c in itertools.product(Fplant, Fsoil, Fanthro, Fbounds)]
all_combos = [[this_run for this_run in this_set if this_run is not None]
              for this_set in all_combos]

aqcs_all = {"-".join(this_set): AqoutContainerSpatialPaper(
    aqout_paths=[runs[this_run].aqout_path for this_run in this_set],
    aq_keys=this_set,
    key="-".join(this_set)) for this_set in all_combos}


aqcs = dict((k, aqcs_all[k]) for k in aqcs_all.keys() if 'SiB' in k)


for k in [aqcs.keys()[i] for i in [0, 1]]:
    print "processing {}".format(k)
    aqcs[k].parse(t0=datetime(2008, 7, 8),
                  t1=datetime(2008, 8, 31, 0, 0, 0),
                  verbose=True)
    mean_dd = []
    for i, aqdata in enumerate(aqcs[k].data):
        print i
        dd = calc_STEM_COS_drawdown(aqdata,
                                    topo_fname=topo_file,
                                    wrfheight_fname=wrf_file)
        del aqdata
        this_mean_dd = dd[:, :, x_site, y_site].mean(axis=0).squeeze()
        del dd
        gc.collect()
        mean_dd.append(this_mean_dd)
    # site_columns = [arr[:, :, x_site, y_site].mean(axis=0).squeeze()
    #                 for arr in aqcs[k].data]
    # aqcs[k].components = pd.DataFrame(dict(zip(aqcs[k].key.split('-'),
    #                                            site_columns)))

    # this_aqc.sum()
    # # this_aqc.calc_stats()
    # this_aqc.calc_JA_midday_drawdown()
    # this_aqc.calc_JA_midday_drawdown_stderr()
    # this_aqc.extract_noaa_sites(
    #     '/project/projectdirs/m2319/Data/NOAA_95244993/')

# all = pd.concat([this_model.site_vals for this_model in aqcs.values()])
# all.to_csv('./model_components_25Feb.csv')
# print datetime.now() - t0
