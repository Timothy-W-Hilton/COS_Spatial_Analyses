import os
import os.path
import numpy as np
import netCDF4

from climatological_bounds import NOAASite
from IOAPIpytools.ioapi_pytools import boundaries_from_csv


if __name__ == "__main__":
    PFA = NOAASite('PFA', 493.54335, 26.31656333)
    ESP = NOAASite('ESP', 495.7402393, 6.118293878)
    THD = NOAASite('THD', 538.1678958, 25.6922375)
    TGC = NOAASite('TGC', 519.1636667, -12.30926667)
    NHA = NOAASite('NHA', 466.444445, 75.521758)
    SCA = NOAASite('SCA', 510.726333, 24.069413)
    CMA = NOAASite('CMA', 490.787018,  56.442238)

    E_mean = NOAASite('E_mean', 489.319265, 52.011136)

    # starting in "lower left" with SW corner of domain and going counter
    # clockwise, pfa could do north and northern pacific, esp a little
    # lower on the pacific, and thd for rest of pacific and southwestern,
    # tgc for rest of south, and maybe an average of nha/sca/cma for the
    # east (which shouldn't matter).

    bounds = np.hstack((np.tile(TGC.column_cos[:, np.newaxis], 31),
                        np.tile(E_mean.column_cos[:, np.newaxis], 93),
                        np.tile(PFA.column_cos[:, np.newaxis], 126 + 42),
                        np.tile(ESP.column_cos[:, np.newaxis], 42),
                        np.tile(THD.column_cos[:, np.newaxis], 42 + 42),
                        np.tile(TGC.column_cos[:, np.newaxis], 82)))

    pptv_2_molecules_m3 = 1e-12
    pptv_2_ppbv = 1e-3
    fname_csv = 'simple_climatological_bounds.csv'
    fname_bdy = 'climatological_COS_bdy_22levs_124x124.nc'
    fname_griddesc = os.path.join(os.environ['HOME'],
                                  'Code', 'Regrid',
                                  'GEOS-Chem_Regrid', 'GRIDDESC_GC')
    nlevs = 22
    fdesc = ("PFA for N and N pacific, ESP a little "
             "lower on the pacific, and THD for rest of pacific and "
             "SW, TGC for rest of S, and NHA/CMA/SCA mean for the E "
             "(which shouldn't matter).")
    # delete_if_exists(fname_csv)
    # np.savetxt(fname_csv,
    #            bounds.reshape([-1, 1]) * pptv_2_molecules_m3,
    #            delimiter=',')

    if True:
        boundaries_from_csv(fname_csv, fname_bdy,
                            fname_griddesc, 'ARCNAGRID', nlevs, fdesc)

        # place the climatological bounds in the dummy boundary file
        nc = netCDF4.Dataset(fname_bdy, 'a')
        nc.variables['CO2_TRACER1'][...] = (bounds[np.newaxis, ...] *
                                            pptv_2_ppbv)
        nc.close()
