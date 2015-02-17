import os
import os.path
import numpy as np
import numpy.ma as ma
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import maskoceans
from matplotlib.ticker import FuncFormatter

from timutils import midpt_norm, scinot_format
from stem_pytools import STEM_parsers as sp
from stem_pytools import na_map


def get_WRF_Tsoil_VWC(fname_wrf):
    """obtain WRF soil T and soil moisture.  soil T is masked below 10 C
    because Mary's soil flux model fitting data went no lower..
    """

    vwc = sp.parse_STEM_var(nc_fname=fname_wrf, varname='SMOIS')
    Tsoil = sp.parse_STEM_var(nc_fname=fname_wrf, varname='TSOIL')

    ten_C = 273.15 + 10  # 10 C expressed in Kelvins
    Tsoil['data'] = ma.masked_less(Tsoil['data'], ten_C)

    return(vwc, Tsoil)


def calc_fsoil(vwc, Tsoil):
    """implement Mary Whelan's soil COS flux model described in her email
    of 5 Feb 2015.

    INPUTS
       vwc: soil volumetric water content [fraction]
       Tsoil: soil temperature [k]

    OUTPUTS:
       fsoil: soil COS flux [pmol/m2/sec]
    """

    fsoil = (vwc * -28.77873448 + (Tsoil * 0.88867741) - 252.76497309)

    return(fsoil)


def integrate_mary_fsoil(fsoil, s_per_tstamp):
    """
    convert COS flux from pmol m-2 s-1 to mol m-2
    """
    mol_per_pmol = 1e-12
    n_months = 2  # July and Aug

    fsoil_itgd = np.sum(fsoil * mol_per_pmol * s_per_tstamp,
                        axis=0)
    fsoil_itgd = (fsoil_itgd / n_months).squeeze()
    return(fsoil_itgd)


def get_kettle_soil(fname_kettle_fcos):
    s_per_tstamp = 60 * 60 * 24  # one day expressed as seconds
    n_months = 2  # July and Aug
    fsoil_k = sp.parse_STEM_var(fname_kettle_fcos,
                                t0=datetime(2008, 7, 1),
                                t1=datetime(2008, 8, 31, 23, 59, 59),
                                varname='cos')

    fsoil_k['data'] = ma.masked_invalid(fsoil_k['data'])
    # convert mol m-2 s-1 to mol m-2 mon-1
    fsoil_k_itgd = np.sum(fsoil_k['data'] * s_per_tstamp, axis=0)
    fsoil_k_itgd = (fsoil_k_itgd / n_months).squeeze()
    return(fsoil_k_itgd)


def calc_ratio(fsoil_mary, fsoil_kettle):
    lon, lat, topo = sp.parse_STEM_coordinates(
        os.path.join(os.environ['SARIKA_INPUT'], 'TOPO-124x124.nc'))
    fsoil_mary = maskoceans(lon, lat, fsoil_mary)
    fsoil_kettle = maskoceans(lon, lat, fsoil_kettle)
    ratio = ma.masked_invalid(fsoil_kettle) / ma.masked_invalid(fsoil_mary)
    return(ratio)


def draw_ratio(map, ratio):
    lon, lat, topo = sp.parse_STEM_coordinates(
        os.path.join(os.environ['SARIKA_INPUT'], 'TOPO-124x124.nc'))
    ratio = maskoceans(lon, lat, ratio)
    ratio_norm = midpt_norm.MidpointNormalize(midpoint=1.0)
    cm = map.map.pcolor(lon, lat, ratio,
                        vmin=np.percentile(ratio, 1),
                        vmax=np.percentile(ratio, 99),
                        cmap=plt.get_cmap('PuOr'),
                        norm=ratio_norm,
                        latlon=True)
    cb = plt.colorbar(cm, ax=map.ax_map, extend='both',
                      ticks=np.arange(-7, 9, 2))
    cb.solids.set_edgecolor("face")


def draw_fsoil(map, fsoil, vmin, vmax):
    pmol_per_mol = 1e12
    fsoil = fsoil * pmol_per_mol
    vmin = vmin * pmol_per_mol
    vmax = vmax * pmol_per_mol
    lon, lat, topo = sp.parse_STEM_coordinates(
        os.path.join(os.environ['SARIKA_INPUT'], 'TOPO-124x124.nc'))
    fsoil = maskoceans(lon, lat, fsoil)
    norm = midpt_norm.MidpointNormalize(midpoint=0.0)
    cm = map.map.pcolor(lon, lat, fsoil,
                        vmin=vmin,  # np.nanmin(fsoil),
                        vmax=vmax,  # np.nanmax(fsoil),
                        norm=norm,
                        cmap=plt.get_cmap('RdGy_r'),
                        latlon=True)
    cb = plt.colorbar(cm, ax=map.ax_map, extend='both',
                      format=FuncFormatter(scinot_format.scinot_format))
    cb.solids.set_edgecolor("face")
    cb.ax.set_title('pmol COS m$^{-2}$ mon$^{-1}$',
                    fontdict={'fontsize': 8})


def draw_fsoil_maps(fsoil_mary, fsoil_kettle):

    # NAMapFigure arguments to zoom map on the Eastern USA
    E_USA = {'lon_0': -88.6275,
             'lat_0': 37.0722,
             'mapwidth': 3e6,
             'mapheight': 2.5e6}
    kwargs = E_USA

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18, 4))
    map_m = na_map.NAMapFigure(t_str="Mary's COS F$_{soil}$",
                               map_axis=ax[0],
                               cb_axis=None,
                               **kwargs)
    map_k = na_map.NAMapFigure(t_str="Kettle's COS F$_{soil}$",
                               map_axis=ax[1],
                               cb_axis=None,
                               **kwargs)
    map_r = na_map.NAMapFigure(t_str="Kettle's F$_{soil}$ / Mary's F$_{soil}$",
                               map_axis=ax[2],
                               cb_axis=None,
                               **kwargs)

    vmin = np.nanmin(np.vstack((fsoil_mary, fsoil_kettle)).flatten())
    vmax = np.nanmax(np.vstack((fsoil_mary, fsoil_kettle)).flatten())
    draw_fsoil(map_m, fsoil_mary, vmin, vmax)
    draw_fsoil(map_k, fsoil_kettle, vmin, vmax)

    # draw_fsoil(map_m, fsoil_mary, np.nanmin(fsoil_kettle), 6e-5)
    # draw_fsoil(map_k, fsoil_kettle, np.nanmin(fsoil_kettle), 6e-5)
    ratio = calc_ratio(fsoil_mary, fsoil_kettle)
    draw_ratio(map_r, ratio)
    return(fig, map_m, map_k, map_r)

if __name__ == "__main__":
    s_per_6hrs = 6 * 60 * 60  # six hours expressed in seconds

    fname_wrf = os.path.join(os.environ['SARIKA_INPUT'],
                             'soil_T_moisture_JulAug.nc')
    vwc, Tsoil = get_WRF_Tsoil_VWC(fname_wrf)
    fsoil = calc_fsoil(vwc['data'], Tsoil['data'])
    fsoil_itgd = integrate_mary_fsoil(fsoil, s_per_6hrs)

    fsoil_k_itgd = get_kettle_soil(os.path.join(
        os.environ['SARIKA_INPUT'],
        'surfem-124x124-kettle-soil-cos_2008_2009.nc'))

    plt.close('all')
    fig, map_m, map_k, map_r = draw_fsoil_maps(fsoil_itgd, fsoil_k_itgd)

    fig.savefig('/tmp/soil_model_maps_EUSA.pdf')
