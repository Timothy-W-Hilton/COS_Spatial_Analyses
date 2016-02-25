import matplotlib
matplotlib.use('AGG')
import pandas as pd
from datetime import datetime

from stem_pytools import NERSC_data_paths as ndp
from aqout_postprocess_spatial_paper import AqoutContainerSpatialPaper
from stem_pytools import STEM_mapper
import itertools

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

aqcs = {"-".join(this_set): AqoutContainerSpatialPaper(
    aqout_paths=[runs[this_run].aqout_path for this_run in this_set],
    aq_keys=this_set,
    key="-".join(this_set)) for this_set in all_combos}
#aqcs = {k: v for k, v in aqcs_all.items()[0:2]}

# for k, this_aqc in aqcs.items()[(1, -1)]:
# aqcs = dict((k, aqcs_all[k]) for k in ['casa_gfed_161'])
# aqcs = dict((k, aqcs_all[k]) for k in ['casa_gfed_161', 'canibis_161'])

for k, this_aqc in aqcs.items():
    print "processing {}".format(k)
    # this_aqc.parse(t0=None, t1=None, verbose=True)
    this_aqc.parse(t0=datetime(2008, 7, 8),
                   t1=datetime(2008, 8, 31, 0, 0, 0),
                   verbose=True)
    this_aqc.sum()
    # this_aqc.calc_stats()
    this_aqc.calc_JA_midday_drawdown()
    this_aqc.calc_JA_midday_drawdown_stderr()
    this_aqc.extract_noaa_sites(
        '/project/projectdirs/m2319/Data/NOAA_95244993/')

all = pd.concat([this_model.site_vals for this_model in aqcs.values()])
all.to_csv('./model_components_25Feb.csv')
print datetime.now() - t0

# se_map = STEM_mapper.Mapper124x124(aqcs['canibis_161'].dd_se)
# se_map.draw_map(t_str=('Can-IBIS LRU=1.61 [COS] drawdown '
#                        'Jul/Aug mean standard error'))
# se_map.map.fig.savefig('./CanIBIS_se_map.png')

# se_map = STEM_mapper.Mapper124x124(aqcs['canibis_161'].dd_se_neff)
# se_map.draw_map(t_str=('Can-IBIS LRU=1.61 [COS] drawdown '
#                        'Jul/Aug mean standard error, Neff'))
# se_map.map.fig.savefig('./CanIBIS_se_neff_map.png')
