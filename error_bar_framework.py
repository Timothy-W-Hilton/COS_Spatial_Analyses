import os.path

from stem_pytools import NERSC_data_paths as ndp
from aqout_postprocess_spatial_paper import AqoutContainerSpatialPaper
from timutils import std_error

runs = ndp.get_Spatial_Paper_runs()
aqcs = {k: AqoutContainerSpatialPaper(v.aqout_path) for k, v in runs.items()}
fnames = []
dd = []
# for k, this_aqc in aqcs.items()[(1, -1)]:
subset = dict((k, aqcs[k]) for k in ['casa_gfed_161', 'casa_m15_C4pctLRU'])
subset = dict((k, aqcs[k]) for k in ['casa_gfed_161'])
for k, this_aqc in subset.items():
    this_aqc.parse(t0=None, t1=None, verbose=True)
    this_aqc.sum()
    this_aqc.calc_stats()
    dd.append(this_aqc.calc_JA_midday_drawdown())

se = std_error.MeanStdError(dd[0])
se.calc()

#     fnames.append(os.path.join(os.getenv('SCRATCH'), 'AQOUT_stats',
#                                'aqout_{}.nc'.format(k)))
#     print "writing {}".format(fnames[-1])
#     this_aqc.stats_to_netcdf(fnames[-1])

# aqpp.combine_aqout_netcdf_files(fnames,
#                                 fname_out=os.path.join(os.getenv('SCRATCH'),
#                                                        'AQOUT_stats',
#                                                        'AQ_combined.nc'))
