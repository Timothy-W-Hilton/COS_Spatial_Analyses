"""Calculate daily [COS] mean and standard deviation for:

(1) C4 percentage LRU STEM runs for all five GPP models (CASA-GFED,
CASA-m15, Can-IBIS, Kettle, and MPI).
(2) Fsoil STEM runs for Kettle Fsoil and Whelan-Kettle hybrid Fsoil
(3) all combinations of (1) and (2)

Also serves to demonstrate how to calculate daily [COS] mean and
standard deviation from one or more STEM AQOUT files.
"""

import os.path
import sys
from datetime import datetime
from itertools import product
from stem_pytools import aqout_postprocess as aq


def process(aqc):
    """
    given an aq_container object, parse the AQOUT files, sum them, and
    calculate July-August mid-day [COS] mean and standard deviation.
    """

    jul1 = datetime(2008, 7, 1)
    aug31 = datetime(2008, 8, 31, 23, 59, 59)

    t0 = datetime.now()
    aqc.parse(jul1, aug31, verbose=True)
    sys.stdout.write('done parsing ({}s)\n'.format(
        (datetime.now() - t0).seconds))
    sys.stdout.flush()

    aqc.sum()

    t0 = datetime.now()
    sys.stdout.write('calculating stats')
    sys.stdout.flush()
    aqc.calc_stats()
    sys.stdout.write('done ({}s)\n'.format((datetime.now() - t0).seconds))
    sys.stdout.flush()

    outfile = 'AQOUTagg_{}.nc'.format(aqc.key)
    aqc.stats_to_netcdf(outfile)

if __name__ == "__main__":

    stem_out_root = os.path.join('/',
                                 'Users',
                                 'tim',
                                 'work',
                                 'Data',
                                 'STEM',
                                 'output')

    aq_Fsoil_kettle = os.path.join(
        stem_out_root,
        'STEM_Runs_Fsoil',
        'Kettle_Fsoil',
        'output',
        'AQOUT-124x124-22levs-Kettle_Fsoil.nc')
    aq_Fsoil_hybrid = os.path.join(
        stem_out_root,
        'STEM_Runs_Fsoil',
        'Kettle_Whelan_hybrid_Fsoil',
        'output',
        'AQOUT-124x124-22levs-hybrid_Fsoil_5FebModel.nc')

    fsoil_run_aq = (aq_Fsoil_kettle, aq_Fsoil_hybrid)
    fsoil_run_keys = ('KettleFsoil', 'Whelan5FebKettleFsoil')
    fsoil_run_desc = ('Kettle Fsoil',
                      'Whelan-Kettle hybrid Fsoil model of 5 Feb email')

    aq_casagfed = os.path.join(
        stem_out_root, 'CASA-GFED_C4pctLRU', 'output',
        'AQOUT-124x124-22levs-CASA-GFED_fCOS_C4pctLRU.nc')
    aq_casam15 = os.path.join(
        stem_out_root, 'CASA-m15_C4pctLRU', 'output',
        'AQOUT-124x124-22levs-CASA-m15_fCOS_C4pctLRU.nc')
    aq_kettle = os.path.join(
        stem_out_root, 'Kettle_C4pctLRU', 'output',
        'AQOUT-124x124-22levs-kettle_fCOS_C4pctLRU.nc')
    aq_canibis = os.path.join(
        stem_out_root, 'Can-IBIS_C4pctLRU', 'output',
        'AQOUT-124x124-22levs-CanIBIS_fCOS_C4pctLRU.nc')
    aq_mpi = os.path.join(
        stem_out_root, 'MPI_C4pctLRU', 'output',
        'AQOUT-124x124-22levs-MPI_fCOS_C4pctLRU.nc')

    fplant_run_aq = (aq_casagfed, aq_casam15, aq_canibis,
                     aq_kettle, aq_mpi)
    fplant_run_keys = ('CASAGFEDFplant_C4pctLRU', 'CASAm15Fplant_C4pctLRU',
                       'CanIBISFplant_C4pctLRU', 'KettleFplant_C4pctLRU',
                       'MPIFplant_C4pctLRU')
    fplant_run_desc = ('CASA-GFED C4 pct LRU', 'CASA-m15 C4 pct LRU',
                       'Can-IBIS Fplant C4 pct LRU',
                       'Kettle Fplant C4 pct LRU',
                       'MPI Fplant C4 pct LRU')

    plant_soil_aq = [p for p in product(fplant_run_aq, fsoil_run_aq)]
    plant_soil_keys = [''.join(p) for p in product(fplant_run_keys,
                                                   fsoil_run_keys)]
    plant_soil_desc = [' combined with '.join(p)
                       for p in product(fplant_run_desc,
                                        fsoil_run_desc)]

    plant_soil_aqc = [aq.aqout_container(*args)
                      for args in zip(plant_soil_aq,
                                      plant_soil_desc,
                                      plant_soil_keys)]

    Fsoil_aqc = [aq.aqout_container(*args)
                 for args in zip(fsoil_run_aq,
                                 fsoil_run_desc,
                                 fsoil_run_keys)]

    Fplant_aqc = [aq.aqout_container(*args)
                  for args in zip(fplant_run_aq,
                                  fplant_run_desc,
                                  fplant_run_keys)]

    # for this in (plant_soil_aqc + Fsoil_aqc + Fplant_aqc):
    #     process(this)
    for this in (Fsoil_aqc + Fplant_aqc):
        process(this)
