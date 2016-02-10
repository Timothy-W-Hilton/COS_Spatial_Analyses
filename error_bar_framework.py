import os.path
import os

from stem_pytools import NERSC_data_paths as ndp
from stem_pytools import aqout_postprocess as aqpp

runs = ndp.get_Spatial_Paper_runs()
aqcs = {k: aqpp.aqout_container(v.aqout_path) for k, v in runs.items()}
for k, this_aqc in aqcs.items()[0:3]:
    this_aqc.parse(t0=None, t1=None, verbose=True)
    this_aqc.sum()
    this_aqc.calc_stats()
    fname = os.path.join(os.getenv('SCRATCH'), 'aqout_{}.nc'.format(k))
    print "writing {}".format(fname)
    this_aqc.stats_to_netcdf(fname)
