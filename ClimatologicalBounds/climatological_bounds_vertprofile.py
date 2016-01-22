import os
import os.path
import numpy as np
import netCDF4
import pandas as pd

from climatological_bounds import NOAASite
from IOAPIpytools.ioapi_pytools import boundaries_from_csv
from stem_pytools import noaa_ocs
from stem_pytools import domain


def rm_nan(arr):
    """remove NaNs from a numpy array.  Return a flattened (i.e. 1-D)
    array containing all non-NaN values from arr.
    """
    return(arr[np.isfinite(arr)].flatten())


class site_clim_mean(object):
    """class to calculate climatological mean vertical OCS profile"""

    def __init__(self, sitecode, noaa_dir=None):
        if noaa_dir is None:
            noaa_dir = os.path.join(os.getenv('PROJ'), 'Data', 'NOAA_95244993')
        fname = os.path.join(
            noaa_dir,
            'ocs_{}_aircraft-pfp_1_hats_event.txt'.format(sitecode))
        self.noaa_site = noaa_ocs.NOAA_OCS.parse_file(fname)
        self.z_obs_mean = None
        self.z_all_agl = None
        self.x_stem = None
        self.y_stem = None

        self.get_jul_aug()
        self.get_z_lev_mean()
        self.get_all_z_agl()

    def get_jul_aug(self):
        jul_aug = self.noaa_site.obs.query('sample_month in [7, 8]')
        self.noaa_site.obs = jul_aug

    def get_z_lev_mean(self):
        d = domain.STEM_Domain()
        topo_dir = os.path.join(os.getenv('PROJ'), 'Data',
                                'STEM_124x124_NA_inputs')
        self.noaa_site.get_stem_xy(d.get_lon(), d.get_lat())
        self.noaa_site.get_stem_z(
            topo_fname=os.path.join(topo_dir, 'TOPO-124x124.nc'),
            wrfheight_fname=os.path.join(topo_dir,
                                         'wrfheight-124x124-22levs.nc'))
        self.z_obs_mean = self.noaa_site.obs.groupby('z_stem').mean()
        self.z_obs_mean = self.z_obs_mean[['x_stem', 'y_stem',
                                           'sample_altitude',
                                           'analysis_value']]
        # change the index (which is the z level after the groupby to
        # a column)
        self.z_obs_mean.reset_index(level=0, inplace=True)
        # a handful of obs are in adjacent STEM cells, resulting in
        # non-integral mean x or y cell locations after the mean is
        # taken.  I think that rounding will pick the the "mode" x and
        # y and convert to an integer in one step
        self.x_stem = np.int(np.unique(np.round(self.z_obs_mean['x_stem']))[0])
        self.y_stem = np.int(np.unique(np.round(self.z_obs_mean['y_stem']))[0])
        self.z_obs_mean['x_stem'][:] = self.x_stem
        self.z_obs_mean['y_stem'][:] = self.y_stem

    def get_all_z_agl(self):
        """get all STEM Z cell heights above ground level (from surface to
        top of domain).
        """
        topo_fname = os.path.join(os.getenv('SARIKA_INPUT'), 'TOPO-124x124.nc')
        wrfheight_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                                       'wrfheight-124x124-22levs.nc')
        dom = domain.STEM_Domain(fname_topo=topo_fname)
        dom.get_STEMZ_height(wrfheight_fname)

        if self.z_obs_mean is None:
            self.get_z_lev_mean()
        n_zlevs = dom.agl.shape[0]

        z_all_agl = pd.DataFrame({'z_stem': np.arange(n_zlevs) + 1,
                                  'z_agl': dom.agl[:, self.x_stem,
                                                   self.y_stem]})
        self.z_obs_mean = pd.merge(z_all_agl, self.z_obs_mean,
                                   how='outer', sort=True)
        self.z_obs_mean['ocs_interp'] = self.z_obs_mean['analysis_value']

        ocs_interp = np.interp(
            self.z_obs_mean.z_agl,
            rm_nan(self.z_obs_mean.sample_altitude.values),
            rm_nan(self.z_obs_mean.analysis_value.values),
            left=self.z_obs_mean.analysis_value.min(),
            right=self.z_obs_mean.analysis_value.max())
        nan_idx = np.where(np.isnan(self.z_obs_mean['analysis_value']))
        self.z_obs_mean['ocs_interp'].iloc[nan_idx] = ocs_interp[nan_idx]

if __name__ == "__main__":

    thd = site_clim_mean('thd')
    pfa = site_clim_mean('pfa')
    esp = site_clim_mean('esp')
    tgc = site_clim_mean('tgc')
    nha = site_clim_mean('nha')
    sca = site_clim_mean('sca')
    cma = site_clim_mean('cma')

    # foo.get_jul_aug()
    # foo.get_z_lev_mean()
    # foo.get_all_z_agl()

    if False:
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
