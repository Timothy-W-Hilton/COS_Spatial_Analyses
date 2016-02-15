"""Extension to stem_pytools.aqout_postprocess.aqout_container with
july/aug mean COS drawdown calculation with uncertainty (functionality
specific to the spatial paper).
"""

import numpy as np
from stem_pytools import aqout_postprocess as aqpp
from timutils import std_error


class AqoutContainerSpatialPaper(aqpp.aqout_container):
    """Derived class implements July-August-specific functionality

    This functionality is specifically pertinent to the spatial paper;
    it populates fields dd_JA_midday, dd_se, dd_se_neff, dd_neff

    dd_JA_midday (numpy array): dimensions [time, x, y]; contains
       July-August mean COS vertical drawdown in each STEM horizontal
       grid cell
    dd_se: (numpy array): dimensions [time, x, y]; July-August mean
       COS drawdown standard error
    dd_se_neff: (numpy array): dimensions [time, x, y]; July-August mean
       COS drawdown standard error calculated using the "effective
       sample size" adjusted for autocorrelation within the STEM [COS]
    dd_neff: (numpy array): dimensions [time, x, y]; "effective sample
       size" adjusted for autocorrelation within the STEM [COS].  The
       number of effectively independent [COS] data points from STEM.
    """
    def calc_JA_midday_drawdown(self):
        """calculates and populates fields dd_JA_midday (see
        AqoutContainerSpatialPaper docstring)
        """

        is_midday_bool = np.array(map(aqpp.is_midday, self.t))
        is_not_midday = np.ones_like(self.cos_total)
        is_not_midday[is_midday_bool, ...] = 0
        self.cos_total = np.ma.masked_where(is_not_midday == 1, self.cos_total)
        self.dd_JA_midday = self.calc_drawdown().squeeze()

    def calc_JA_midday_drawdown_stderr(self):
        """calculates and populates fields, dd_se, dd_se_neff, dd_neff (see
        AqoutContainerSpatialPaper docstring)
        """

        if self.dd_JA_midday is None:
            self.calc_JA_midday_drawdown()

        se = std_error.MeanStdError(self.dd_JA_midday)
        se.calc()

        self.dd_se = se.std_err
        self.dd_se_neff = se.std_err_neff
        self.dd_neff = se.neff
