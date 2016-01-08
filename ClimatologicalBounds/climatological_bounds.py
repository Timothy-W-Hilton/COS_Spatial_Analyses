import numpy as np


class NOAASite(object):
    """
    Class to create climatological boundaries from NOAA airborne observations

    Class attributes:
       obs: pandas DataFrame containing NOAA observations
       obs_color: matplotlib color to use for plotting observations.
       grid_color: matplotlib color to use for plotting STEM grid cell points
    """
    def __init__(self, site_name, free_trop_cos, drawdown, nlevs=22):
        self.site_name = site_name
        self.free_trop_cos = free_trop_cos
        self.drawdown = drawdown
        self.nlevs = nlevs
        self.column_cos = np.linspace(start=self.free_trop_cos - self.drawdown,
                                      stop=self.free_trop_cos,
                                      num=self.nlevs)
