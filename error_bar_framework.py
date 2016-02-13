from stem_pytools import NERSC_data_paths as ndp
from aqout_postprocess_spatial_paper import AqoutContainerSpatialPaper

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
