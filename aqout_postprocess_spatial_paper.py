"""Extension to stem_pytools.aqout_postprocess.aqout_container with
july/aug mean COS drawdown calculation with uncertainty (functionality
specific to the spatial paper).
"""

import numpy as np
from stem_pytools import aqout_postprocess as aqpp
from stem_pytools import noaa_ocs
from stem_pytools import domain as domain_tools
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
    def parse(self, const_bounds=450, *args, **kwargs):
        """Parse aqout data by calling parent class parse method.  Then apply:

        (1) corrections to anthropogenic fluxes (multiply by 1e3, as
        per email from Andrew Zumkehr).  Anthro STEM runs are
        identified by the substring "Anthro" existing in self.key
        (2) adjustments to climatic boundaries totals: substract the
        constant bounds [COS] of 450 ppt from all non-boundary run
        AQOUT files, then add the bounds STEM run [COS] to each, then
        remove the bounds run AQOUT from self.

        ARGS:
        const_bounds (scalar): the constant [COS] (in pptv) used in
           constant boundaries STEM runs
        """
        # call the parent class parse method
        super(AqoutContainerSpatialPaper, self).parse(*args, **kwargs)

        for i, this_key in enumerate(self.aq_keys):
            print(('multiplying Anthro [COS] by 1000 '
                   'as per email from Andrew'))
            if 'Anthro' in this_key:
                self.data[i] = self.data[i] * 1e3

        for i, this_key in enumerate(self.aq_keys):
            if 'climatological_bnd' in this_key:
                print 'adjusting {} for climatolotical bounds'.format(self.key)
                for j in range(len(self.data)):
                    if (j != i):
                        self.data[j] = (self.data[j] -
                                        const_bounds + self.data[i])

    def calc_JA_midday_drawdown(self):
        """calculates and populates fields dd_JA_midday (see
        AqoutContainerSpatialPaper docstring)
        """

        is_midday_bool = np.array(map(aqpp.is_midday, self.t))
        is_not_midday = np.ones_like(self.cos_total)
        is_not_midday[is_midday_bool, ...] = 0
        self.cos_total = np.ma.masked_where(is_not_midday == 1, self.cos_total)
        self.dd_JA_midday = self.calc_drawdown().squeeze()
        self.dd_JA_midday_mean = self.dd_JA_midday.mean(axis=0)

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

    def extract_noaa_sites(self, noaa_dir):
        """extract drawdown, standard error for each NOAA observation site

        populates field noaa_site_vals with a Pandas data frame
        indexed by three-letter NOAA site code, with columns dd
        (drawdown, pptv) and se (standard error, pptv)

        ARGS:
        noaa_dir (string): fully path the to the directory containing
           the noaa data
        """

        site_vals = noaa_ocs.get_sites_summary(noaa_dir)
        domain = domain_tools.STEM_Domain()

        stem_x, stem_y = domain_tools.find_nearest_stem_xy(
            site_vals['longitude'],
            site_vals['latitude'],
            domain.get_lon(),
            domain.get_lat())
        site_vals['stem_x'] = stem_x
        site_vals['stem_y'] = stem_y

        site_vals['dd'] = self.dd_JA_midday_mean[site_vals.stem_x,
                                                 site_vals.stem_y]
        site_vals['dd_neff'] = self.dd_neff[site_vals.stem_x,
                                            site_vals.stem_y]
        site_vals['dd_se'] = self.dd_se[site_vals.stem_x,
                                        site_vals.stem_y]
        site_vals['dd_se_neff'] = self.dd_se_neff[site_vals.stem_x,
                                                  site_vals.stem_y]
        self.site_vals = site_vals
        site_vals.insert(0, 'model', self.key)
