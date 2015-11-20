import subprocess
import matplotlib
matplotlib.use('AGG')


import svgutils.transform as sg
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
    sz = map(int, sg.from_mpl(figs[0]).get_size())
    sz[1] *= 3
    sz = map(str, sz)
    layout.set_size(sz)
    layout.add_figure(sg.from_mpl(figs[0]))
    layout.add_figure(sg.from_mpl(figs[1]))

    txt1 = sg.TextElement(50, 50, "HELLO", size=12)
    layout.append([txt1])
    layout.save(os.path.join('/', 'tmp', 'stack_plots.svg'))
    try:
        print('converting to pdf')
        subprocess.call('/Applications/Inkscape.app/Contents/Resources/bin/inkscape --export-pdf=/tmp/stack_plots.pdf /tmp/stack_plots.svg', shell=True)
    except:
        print('unable to run inkscape')
finally:
    plt.close('all')
