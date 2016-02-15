import matplotlib
matplotlib.use('AGG')

from stem_pytools import NERSC_data_paths as ndp
from aqout_postprocess_spatial_paper import AqoutContainerSpatialPaper
from stem_pytools import STEM_mapper

runs = ndp.get_Spatial_Paper_runs()
aqcs = {k: AqoutContainerSpatialPaper(v.aqout_path) for k, v in runs.items()}
fnames = []

# for k, this_aqc in aqcs.items()[(1, -1)]:
subset = dict((k, aqcs[k]) for k in ['casa_gfed_161', 'casa_m15_C4pctLRU'])
subset = dict((k, aqcs[k]) for k in ['casa_gfed_161'])
for k, this_aqc in subset.items():
    this_aqc.parse(t0=None, t1=None, verbose=True)
    this_aqc.sum()
    # this_aqc.calc_stats()
    this_aqc.calc_JA_midday_drawdown()
    this_aqc.calc_JA_midday_drawdown_stderr()

se_map = STEM_mapper.Mapper124x124(subset['casa_gfed_161'].dd_se)
se_map.draw_map(t_str=('casa-GFED3 LRU=1.61 [COS] drawdown '
                       'Jul/Aug mean standard error'))
se_map.map.fig.savefig('./se_map.png')

se_map = STEM_mapper.Mapper124x124(subset['casa_gfed_161'].dd_se_neff)
se_map.draw_map(t_str=('casa-GFED3 LRU=1.61 [COS] drawdown '
                       'Jul/Aug mean standard error, Neff'))
se_map.map.fig.savefig('./se_neff_map.png')
