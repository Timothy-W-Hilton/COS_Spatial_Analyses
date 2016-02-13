"""Extension to stem_pytools.aqout_postprocess.aqout_container with
july/aug mean COS drawdown calculation with uncertainty (functionality
specific to the spatial paper).
"""

import numpy as np
from stem_pytools import aqout_postprocess as aqpp


class AqoutContainerSpatialPaper(aqpp.aqout_container):
    """Derived class in which to implement July-August-specific functionality
    """
    def calc_JA_midday_drawdown(self):

        midday_mask = np.array(map(aqpp.is_midday, self.t))
        self.cos_total = self.cos_total[midday_mask, ...]
        dd = self.calc_drawdown().squeeze()
        return dd
