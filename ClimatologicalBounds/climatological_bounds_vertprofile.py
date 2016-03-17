"""Classes, functions, and main function to create STEM boundary files
in I/O API format containing vertical column mean, July-August
climatological mean, observed [COS] from NOAA airborne observations.
"""

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import os
import os.path
import numpy as np
import netCDF4
import pandas as pd
import brewer2mpl

from IOAPIpytools.ioapi_pytools import boundaries_from_csv, dummy_top_bounds
from stem_pytools import noaa_ocs
from stem_pytools import domain
from stem_pytools import STEM_mapper
from stem_pytools.STEM_parsers import parse_STEM_var
from timutils import colormap_nlevs


def rm_nan(arr):
    """remove NaNs from a numpy array.

    ARGS:
    arr (array-like): an array of values, possible containing np.NaN

    RETURNS:
    a flattened (i.e. 1-D) array containing all non-NaN values from arr.
    """
    return arr[np.isfinite(arr)].flatten()


class Consts(object):
    """Class to contain several constants useful for the climatological
    bounds

    ARGS:

    fname_griddesc (string): full path to the Models-3 I/O API
        GRIDDESC file describing the STEM domain grid.  Default is
        $HOME/Code/Regrid/GEOS-Chem_Regrid/GRIDDESC_GC
    noaa_dir (str): full path to the directory containing NOAA
        observation files.  Default is $PROJ/Data/NOAA_95244993
    topo_fname (string): full path to the Models-3 I/O API topography
        file.  Default is $SARIKA_INPUT/TOPO-124x124.nc
    wrf_height_fname (string): full path to the Models-3 I/O API WRF height
        file.  Default is $SARIKA_INPUT/wrfheight-124x124-22levs.nc

    ATTRIBUTES:

    fname_griddesc (string): full path to the Models-3 I/O API
        GRIDDESC file describing the STEM domain grid.
    noaa_dir (str): full path to the directory containing NOAA
        observation files.  Default is $PROJ/Data/NOAA_95244993
    topo_fname (string): full path to the Models-3 I/O API topography
        file.  Default is $SARIKA_INPUT/TOPO-124x124.nc
    wrf_height_fname (string): full path to the Models-3 I/O API WRF height
        file.  Default is $SARIKA_INPUT/wrfheight-124x124-22levs.nc
    pptv_2_molecules_m3 (real) = conversion factor for converting
        parts per trillion by volume (pptv) to molecules per m^3
    pptv_2_ppbv = conversion factor for converting parts per trillion
        by volume (pptv) to parts per billion by volume (ppbv)
    """

    def __init__(self,
                 fname_griddesc=os.path.join(os.environ['HOME'],
                                             'Code', 'Regrid',
                                             'GEOS-Chem_Regrid',
                                             'GRIDDESC_GC'),
                 noaa_dir=os.path.join(
                     os.getenv('PROJ'), 'Data', 'NOAA_95244993')):
        self.fname_griddesc = fname_griddesc
        self.noaa_dir = noaa_dir
        self.topo_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                                       'TOPO-124x124.nc')
        self.wrfheight_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                                            'wrfheight-124x124-22levs.nc')
        self.pptv_2_molecules_m3 = 1e-12
        self.pptv_2_ppbv = 1e-3


class SiteClimMean(object):
    """class to calculate climatological mean vertical OCS profile for
    a NOAA observation site.
    """

    def __init__(self, sitecode, alt_bin_size=1000, noaa_dir=None):
        """class constructor.

        ARGS:
        alt_bin_size (int): size of bins to group observation
            altitudes into (meters). Default is 1000 m.
        noaa_dir (str): full path to the directory containing NOAA
            observation files.  Default is $PROJ/Data/NOAA_95244993
        """

        if noaa_dir is None:
            noaa_dir = os.path.join(os.getenv('PROJ'), 'Data', 'NOAA_95244993')
        fname = os.path.join(
            noaa_dir,
            'ocs_{}_aircraft-pfp_1_hats_event.txt'.format(sitecode.lower()))
        self.alt_bin_size = alt_bin_size
        self.sitecode = sitecode
        self.noaa_site = noaa_ocs.NOAA_OCS.parse_file(fname)
        self.z_obs_mean = None
        self.z_all_agl = None
        self.x_stem = None
        self.y_stem = None

        self.get_jul_aug()
        self.get_z_lev_mean()
        self.get_all_z_agl()

    def get_jul_aug(self):
        """pare down observations data frame to July and August
        observations only
        """
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

        bin_min = 0  # bottom of bottom altitude bin, meters
        bin_max = 16000   # top of top altitude bin, meters
        bin_width = 1000  # size of altitude bins, meters
        bin_edges = np.arange(bin_min, bin_max, bin_width)
        self.noaa_site.obs['altitude_bin'] = np.digitize(
            self.noaa_site.obs.sample_altitude, bin_edges)
        groups = self.noaa_site.obs.groupby('altitude_bin')
        self.z_obs_mean = groups.mean()
        self.z_obs_mean = self.z_obs_mean[['x_stem', 'y_stem',
                                           'sample_altitude',
                                           'analysis_value']]
        # # change the index (which is the z level after the groupby to
        # # a column)
        self.z_obs_mean.reset_index(drop=True, inplace=True)

        # a handful of obs are in adjacent STEM cells, resulting in
        # non-integral mean x or y cell locations after the mean is
        # taken.  I think that rounding will pick the the "mode" x and
        # y and convert to an integer in one step
        self.x_stem = np.int(np.unique(np.round(self.z_obs_mean['x_stem']))[0])
        self.y_stem = np.int(np.unique(np.round(self.z_obs_mean['y_stem']))[0])
        self.z_obs_mean['x_stem'][:] = self.x_stem
        self.z_obs_mean['y_stem'][:] = self.y_stem

        self.z_obs_mean['z_stem'] = domain.get_stem_z_from_altitude(
            self.z_obs_mean.sample_altitude.values,
            stem_x=self.z_obs_mean.x_stem.values,
            stem_y=self.z_obs_mean.y_stem.values)

        self.z_obs_mean = self.z_obs_mean[['x_stem', 'y_stem',
                                           'z_stem', 'sample_altitude',
                                           'analysis_value']]
        # after binning some z levels end up with multiple
        # observations; average them together.
        self.z_obs_mean = self.z_obs_mean.groupby('z_stem').mean()
        self.z_obs_mean.reset_index(inplace=True)

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

        idx_min = self.z_obs_mean.sample_altitude.idxmin()
        idx_max = self.z_obs_mean.sample_altitude.idxmax()
        ocs_interp = np.interp(
            self.z_obs_mean.z_agl,
            rm_nan(self.z_obs_mean.sample_altitude.values),
            rm_nan(self.z_obs_mean.analysis_value.values),
            left=self.z_obs_mean.analysis_value[idx_min],
            right=self.z_obs_mean.analysis_value[idx_max])
        nan_idx = np.where(np.isnan(self.z_obs_mean['analysis_value']))
        self.z_obs_mean['ocs_interp'].iloc[nan_idx] = ocs_interp[nan_idx]

    def get_col_vals(self):
        """return array containing all July-August mean vertical column
        values, with missing values filled.  Missing values
        between the lowest (in altitude) observed value and the
        highest (in altitude) observed value are interpolated;
        missing values outside that range are set to the
        higest/lowest (in altitude) observed value.
        """
        return(self.z_obs_mean.ocs_interp.values[:, np.newaxis])


class ClimatologicalTopBound(object):
    """class to create a climatological mean July-August top boundary
    file from a list of NOAA observation sites.
    """

    def __init__(self,
                 domain,
                 sites_list,
                 sites_dict=None,
                 noaa_dir=os.path.join(
                     os.getenv('PROJ'), 'Data', 'NOAA_95244993')):
        """create a ClimatologicalTopBound object

        ARGS:
        domain: stem_domain (stem_pytools.domain.STEM_Domain):
            description of the STEM domain to create a top boundary
            for
        sites_list (list): list of three-letter NOAA site codes to be
            included
        sites_dict (dict): dict containing one SiteClimMean object for
            each site in sites_list, with the keys the codes in
            site_list.
        noaa_dir (string): full path to the directory containing NOAA
            observation files.  Default is $PROJ/Data/NOAA_95244993
        """
        self.noaa_dir = noaa_dir
        self.d = domain
        self.sites_dict = sites_dict
        self.find_nearest_noaa_site()
        self.get_top_bound_field()

    def find_nearest_noaa_site(self):
        """Find the nearest NOAA observation site to each horizontal
        stem grid cell.

        populates self.idx, self.noaa_sites, self.result
        self.noaa_sites (pandas.DataFrame): data frame with columns
           site_code, latitude, and longitude describing the sites
           used in the top boundary
        self.idx (numpy.ndarray): a 2-D numpy array of the same size
           as the STEM horizontal domain containing the index in
           self.noaa_sites of the nearest NOAA observation location
        self.nearest_site_array (numpy.ndarray): a 2-D numpy array of
           the same size as the STEM horizontal domain containing the
           site code of the nearest NOAA observation location
        """

        stem_lon = self.d.get_lon()
        stem_lat = self.d.get_lat()
        self.sites_summary = noaa_ocs.get_sites_summary(self.noaa_dir)
        # drop sites with too few Jul/Aug data (see longer comment below)
        self.sites_summary = self.sites_summary[
            self.sites_summary.site_code != 'WGC']
        self.sites_summary = self.sites_summary[
            self.sites_summary.site_code != 'OIL']
        self.sites_summary = self.sites_summary.reset_index(drop=True)

        self.idx = domain.find_nearest_stem_xy(
            stem_lon,
            stem_lat,
            self.sites_summary.longitude.values,
            self.sites_summary.latitude.values)
        self.nearest_site_array = np.ndarray(
            shape=self.idx[0].shape, dtype=object)
        self.nearest_site_array[:] = ''
        for i in range(stem_lon.shape[0]):
            for j in range(stem_lon.shape[1]):
                self.nearest_site_array[i, j] = self.sites_summary.site_code[
                    self.idx[0][i, j]]

    def get_top_bound_field(self):
        """populate the top boundary array (self.top_bnd) with
        climatological mean [COS] from nearest noaa site at Z = top of
        domain (22 for the 60-km N American domain)
        """
        nx = self.nearest_site_array.shape[0]
        ny = self.nearest_site_array.shape[1]
        self.d.get_STEMZ_height()
        nz = self.d.asl.shape[0] - 1
        self.top_bnd = np.zeros([nx, ny])
        for i in range(nx):
            for j in range(ny):
                self.top_bnd[i, j] = self.sites_dict[
                    self.nearest_site_array[i, j]].z_obs_mean['ocs_interp'][nz]

    def write_ioapi(self, fname_bdy='top_bounds.nc'):
        """write a Models-3 I/O API top boundary file from the
        object's top_bnd field.

        ARGS:
        fname_bdy (string): name of the boundary file to create
        """
        fdesc = "climatological mean [COS] from nearest noaa site at Z = 22)"
        dummy_top_bounds(fname_bdy,
                         Consts().fname_griddesc,
                         'ARCNAGRID',
                         fdesc)

        # place the climatological bounds in the dummy boundary file
        nc = netCDF4.Dataset(fname_bdy, 'a')
        nc.variables['CO2_TRACER1'][...] = (
            self.top_bnd[np.newaxis, np.newaxis, ...] * Consts().pptv_2_ppbv)
        nc.close()

    def map_nearest_noaa_site(self):
        """Plot the top boundary using
        stem_pytools.STEM_mapper.Mapper124x124, and save the map to
        ./top_bnd.pdf
        """
        stem_lon = self.d.get_lon()
        stem_lat = self.d.get_lat()

        top_cmap, top_norm = colormap_nlevs.setup_colormap(
            np.floor(self.top_bnd.min()),
            np.ceil(self.top_bnd.max()),
            nlevs=6,
            cmap=plt.get_cmap('Oranges'),
            extend='neither')

        m = STEM_mapper.Mapper124x124(self.top_bnd).draw_map(
            fast_or_pretty='pretty',
            cmap=top_cmap,
            norm=top_norm,
            t_str='climatological top bounds',
            cbar_fmt_str='%d')
        m.map.fig.set_figheight(8)
        m.map.fig.set_figwidth(8)

        sites_col = "#1b9e77"  # http://colorbrewer2.org dark2[1]
        fontsz = 24
        for this_site in np.unique(self.nearest_site_array):
            idx = np.where(self.sites_summary.site_code == this_site)[0][0]
            this_x, this_y = m.map.map(self.sites_summary.longitude[idx],
                                       self.sites_summary.latitude[idx])
            m.map.ax_map.text(this_x,
                              this_y,
                              this_site,
                              fontsize=fontsz,
                              color=sites_col,
                              horizontalalignment='center',
                              verticalalignment='center')
        # draw bounds between NOAA sites' regions of top bound
        vals = np.unique(self.top_bnd)
        levs = vals[0:-1] + (np.diff(vals) / 2.0)
        m.map.map.contour(stem_lon, stem_lat, self.top_bnd,
                          levels=levs, latlon=True, colors=sites_col)
        # label the color bar
        m.map.ax_cmap.set_title('[COS] (pptv)')
        m.map.ax_cmap.tick_params(labelsize=fontsz)
        m.map.fig.savefig('top_bnd.pdf')


class ClimatologicalLateralBoundNAmerica(object):
    """class to create a climatological mean July-August lateral boundary
    file for the 60-km N Pole Stereographic North America STEM domain
    from sites THD, PFA, ESP, TGC, NHA, SCA, and CMA.

    """

    def __init__(self,
                 sites_dict):
        """set up a ClimatologicalLateralBoundNAmerica instance.

        sites_dict (dict): dict containing one SiteClimMean object for
           THD, PFA, ESP, TGC, NHA, SCA, and CMA.  Other sites may be
           present in the dict; they will be ignored.
        """

        # starting in "lower left" with SW corner of domain and going counter
        # clockwise, pfa could do north and northern pacific, esp a little
        # lower on the pacific, and thd for rest of pacific and southwestern,
        # tgc for rest of south, and maybe an average of nha/sca/cma for the
        # east (which shouldn't matter).

        self.bounds = np.hstack(
            (np.tile(sites_dict['TGC'].get_col_vals(), 31),
             np.tile(sites_dict['NHA'].get_col_vals(), 31),
             np.tile(sites_dict['CMA'].get_col_vals(), 31),
             np.tile(sites_dict['SCA'].get_col_vals(), 31),
             np.tile(sites_dict['PFA'].get_col_vals(), 126 + 42),
             np.tile(sites_dict['ESP'].get_col_vals(), 42),
             np.tile(sites_dict['THD'].get_col_vals(), 42 + 42),
             np.tile(sites_dict['TGC'].get_col_vals(), 82)))

    def write_bounds_ioapi_file(
            self,
            fname_bdy='climatological_COS_bdy_22levs_124x124.nc'):

        """write a lateral boundary file in I/O API format

        ARGS:
        fname_bdy (string): name for the boundary file to create.
        """

        fname_csv = 'simple_climatological_bounds.csv'

        nlevs = 22
        fdesc = ("PFA for N and N pacific, ESP a little "
                 "lower on the pacific, and THD for rest of pacific and "
                 "SW, TGC for rest of S, and NHA/CMA/SCA mean for the E "
                 "(which shouldn't matter).")
        # delete_if_exists(fname_csv)
        # np.savetxt(fname_csv,
        #            bounds.reshape([-1, 1]) * pptv_2_molecules_m3,
        #            delimiter=',')

        # write a "dummy" boundary file filled with 0.0
        boundaries_from_csv(fname_csv, fname_bdy,
                            Consts().fname_griddesc,
                            'ARCNAGRID', nlevs, fdesc)

        # place the climatological bounds in the dummy boundary file
        nc = netCDF4.Dataset(fname_bdy, 'a')
        nc.variables['CO2_TRACER1'][...] = (self.bounds[np.newaxis, ...] *
                                            Consts().pptv_2_ppbv)
        nc.close()


def create_sites_dict(sites_list):
    """create a dict of SiteClimMean object from a list of site codes

    ARGS:
    sites_list (list): list of three-letter NOAA site codes to include
       in the dict.  These will be the keys of the dict.

    RETURNS:
    dict containing one SiteClimMean object for each site in
       sites_list, with the keys the codes in site_list.
    """
    sites_dict = {}
    for s in sites_list:
        try:
            sites_dict.update({s: SiteClimMean(s)})
        except IndexError:
            print "unable to process {}".format(s)
    return sites_dict


def plot_vertical_profiles(sites_list, title_suffix=None):
    """plot vertical profiles of each site in the argument sites_list

    ARGS:
    sites_list (list): list of SiteClimMean objects containing the
        vertical profiles to be plotted.
    title_suffix (string): optional string to appear as the second
        line of the plot title

    RETURNS:
    None
    """
    # ax.set_color_cycle(palettable.colorbrewer.qualitative.Dark2_8.mpl_colors)
    ncolors = 8
    bmap = brewer2mpl.get_map('Set2', 'qualitative', ncolors)
    colors = bmap.mpl_colors
    matplotlib.rcParams['axes.color_cycle'] = colors
    fig, ax = plt.subplots(figsize=(10, 10))
    for i, this_site in enumerate(sites_list):
        ax.plot(this_site.z_obs_mean.ocs_interp,
                this_site.z_obs_mean.z_agl,
                label=this_site.sitecode,
                linewidth=2)
        ax.scatter(this_site.z_obs_mean.analysis_value,
                   this_site.z_obs_mean.z_agl,
                   marker='o', s=60, c=colors[i % ncolors])
    ax.set_title(('NOAA sites Jul-Aug climatological mean [COS],'
                  ' {} m altitude bins\n{}'.format(
                      sites_list[0].alt_bin_size,
                      title_suffix)))
    ax.set_ylabel('height above ground (m)')
    ax.set_xlabel('[COS] (pptv)')
    ax.set_ylim([0, 16000])
    ax.legend(loc='best')
    fig.savefig('jul_aug_column_profiles0.pdf')
    plt.close(fig)


def top_bounds_QC(top_fname):
    """plot a map of top boundary file to make sure orientation is correct"""
    cos = parse_STEM_var('upbound_124x124-climatological_124x124.nc',
                         varname='CO2_TRACER1')
    m = STEM_mapper.Mapper124x124(cos['data'].squeeze())
    m.draw_map(fast_or_pretty='pretty',
               t_str='Climatological top bounds I/O API',
               cmap=plt.get_cmap('Blues'))
    m.map.fig.savefig('./top_bounds_ioapi_map.png')
    plt.close(m.map.fig)

if __name__ == "__main__":

    # --
    # get STEM domain parameteres
    d = domain.STEM_Domain(Consts().topo_fname)

    # --
    # read NOAA [COS] observations data
    sites = noaa_ocs.get_all_NOAA_airborne_data(Consts().noaa_dir)
    sites_list = list(sites.obs.sample_site_code.unique())
    # drop WGC because there are no Jul/Aug observations
    if 'WGC' in sites_list:
        sites_list.remove('WGC')
    # drop OIL because of weird-looking column profile from only two
    # days of data, with three nearby sites (AAO, WBI, HIL) with much
    # more data and very different column means.
    if 'OIL' in sites_list:
        sites_list.remove('OIL')
    lateral_bounds_sites_list = ['THD', 'PFA', 'ESP',
                                 'TGC', 'NHA', 'SCA', 'CMA']
    sites_dict = create_sites_dict(sites_list)

    # --
    # create the lateral boundary file for N America
    lat_bnd = ClimatologicalLateralBoundNAmerica(sites_dict)
    lat_bnd.write_bounds_ioapi_file()

    # --
    # create a top boundary file
    top_bnd = ClimatologicalTopBound(d, sites_list, sites_dict)
    top_bnd.write_ioapi(fname_bdy='upbound_124x124-climatological_124x124.nc')
    top_bnd.map_nearest_noaa_site()
    # top_bounds_QC('upbound_124x124-climatological_124x124.nc')

    # plot_vertical_profiles([sites_dict[k] for k in lateral_bounds_sites_list],
    #                        'lateral bounds sites')
    # plot_vertical_profiles([sites_dict[k] for k in sites_list], 'all sites')
    upper_midwest = ['OIL', 'HIL', 'WBI', 'AAO']
    # plot_vertical_profiles([sites_dict[k] for k in upper_midwest],
    #                        'Upper Midwest')
