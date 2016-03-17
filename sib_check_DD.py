import matplotlib
matplotlib.use('AGG')
from datetime import datetime

from stem_pytools import NERSC_data_paths as ndp
from aqout_postprocess_spatial_paper import AqoutContainerSpatialPaper

t0 = datetime.now()

runs = ndp.get_Spatial_Paper_runs()
sib_runs = ['SiB_calc', 'SiB_mech']

aqcs = {this_run: AqoutContainerSpatialPaper(
    aqout_paths=runs[this_run].aqout_path,
    aq_keys=this_run,
    key=this_run) for this_run in sib_runs}

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

print datetime.now() - t0
