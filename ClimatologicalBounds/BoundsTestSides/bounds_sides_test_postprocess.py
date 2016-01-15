import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import netCDF4
from stem_pytools import STEM_mapper


def read_aqout(fname, tstep=-1):
    """read and return surface "slab" at specified timestep from STEM
    AQOUT file (default timestep is the last step in the file)
    """

    nc = netCDF4.Dataset(fname)
    z_sfc = 0  # z = 0 is the surface
    cos = nc.variables['CO2_TRACER1'][tstep, z_sfc, ...].squeeze()
    return(cos)


if __name__ == "__main__":

    suffixes = ['0_124', '125_249', '250_374', '375_499']
    cos = [read_aqout("./STEM_TestRuns/output/AQOUT_{}.nc".format(sfx))
           for sfx in suffixes]

    molecules_m3_to_pptv = 1e12
    for i, this_cos in enumerate(cos):
        this_cos = this_cos * molecules_m3_to_pptv
        m = STEM_mapper.Mapper124x124(this_cos).draw_map(
            t_str="I/O API bounds test, COS in bounds[{}]".format(suffixes[i]),
            fast_or_pretty='pretty',
            cmap=plt.get_cmap('Blues'))
        m.map.fig.savefig('bounds{}.png'.format(suffixes[i]))
