"""Create four test I/O API external boundary files, each with 450
pptv along one side of the rectangular domain and zeros elsewhere.
Running each file through STEM will then show which cells of the
boundary "ring" end up on which side of the domain during a run.  This
information is crucial for developing a climatalogical boundary file.
"""

import netCDF4
import shutil


def create_sides_bounds_files(outfiles):
    """from an all-zeros boundary file (bounds_zeros.nc) create four new
    bounds files, each with 450 pptv along one side of the rectangular
    domain and zeros elsewhere.  A 124x124 rectangular domain has a
    500-element boundary "ring" (124 x 4 plus four corner cells), so
    the four files will contain 450 in cells 0 to 125, 126 to 250, 251
    to 375, and 376 to 499, respectively.
    """

    for f in outfiles:
        shutil.copyfile('bounds_zeros.nc', f)


def fill_one_side(fname, idx0, idx1):
    """fill idx0:idx1 of variable CO2_TRACER1 in boundary file with
    name fname with 450 pptv (0.45 ppbv)"""

    print("replacing {}:{} in {}".format(idx0, idx1, fname))
    nc = netCDF4.Dataset(fname, 'a')
    nc.variables['CO2_TRACER1'][:, :, idx0:idx1] = 0.45  # 450 pptv = 0.45 ppbv
    nc.close()


if __name__ == "__main__":

    outfiles = ['bounds_0_124.nc', 'bounds_125_249.nc',
                'bounds_250_374.nc', 'bounds_375_499.nc']
    create_sides_bounds_files(outfiles)
    idx0 = 0
    idx1 = 124
    for f in outfiles:
        fill_one_side(f, idx0, idx1)
        idx0 += 125
        idx1 += 125
