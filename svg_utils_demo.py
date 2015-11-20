
import matplotlib
matplotlib.use('AGG')


from svgutils.transform import from_mpl
from svgutils.templates import VerticalLayout
import matplotlib.pyplot as plt
import numpy as np
import os
import os.path

try:
    figs = []
    for i in range(2):
        figs.append(plt.figure())
        plt.plot(np.random.random(100))

    layout = VerticalLayout()
    sz = map(int, from_mpl(figs[0]).get_size())
    sz[1] *= 3
    sz = map(str, sz)
    layout.set_size(sz)
    layout.add_figure(from_mpl(figs[0]))
    layout.add_figure(from_mpl(figs[1]))

    print from_mpl(figs[0]).get_size()
    layout.save(os.path.join(os.getenv('HOME'), 'plots', 'stack_plots.svg'))

finally:
    plt.close('all')
