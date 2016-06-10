"""Calculate total size of COS spatial paper data files and copy the
data files to a single directory tree for archiving
"""

import os
import os.path
import shutil
from stem_pytools import NERSC_data_paths as ndp


def Get_Human_Readable(size, precision=2):
    """http://stackoverflow.com/questions/5194057/better-way-to-convert-file-sizes-in-python
    """

    suffixes = ['B', 'KB', 'MB', 'GB', 'TB']
    suffixIndex = 0
    while size > 1024 and suffixIndex < 4:
        suffixIndex += 1  # increment the index of the suffix
        size = size/1024.0  # apply the division
    return "%.*f%s" % (precision, size, suffixes[suffixIndex])


def get_spatial_paper_data_total(runs):
    all_data_sum = 0
    for k, this_run in runs.items():
        for this_file in (this_run.aqout_path, this_run.gpp_path,
                          this_run.gppraw_path, this_run.fcos_path):
            if this_file is not None:
                all_data_sum += os.path.getsize(this_file)

    print "Spatial paper data total: " + Get_Human_Readable(all_data_sum)
    return all_data_sum


def make_data_archive(root_dir, runs):
    """Copy all non-regridded GPP, regridded GPP, STEM AQOUT, and fCOS
    netcdf files to a single directory tree for archiving.
    """
    if os.path.exists(root_dir):
        try:
            shutil.rmtree(root_dir)
        except:
            print "unable to delete".format(root_dir)
    try:
        os.makedirs(root_dir)
    except:
        print "unable to create {}".format(root_dir)

    for k, this_run in runs.items():
        print "copying {} files".format(k)
        this_run_dir = os.path.join(root_dir, k)
        os.makedirs(this_run_dir)
        for this_file in (this_run.aqout_path, this_run.gpp_path,
                          this_run.gppraw_path, this_run.fcos_path):
            if this_file is not None:
                print "    copying {}".format(os.path.basename(this_file))
                shutil.copy(this_file, this_run_dir)
        if k is 'climatological_bnd':
            for this_bnd in (runs[k].top_bounds_path,
                             runs[k].lateral_bounds_path):
                print "    copying {}".format(os.path.basename(this_bnd))
                shutil.copy(this_bnd, this_run_dir)


if __name__ == "__main__":
    runs = ndp.get_Spatial_Paper_runs()
    total = get_spatial_paper_data_total(runs)
    archive_dir = os.path.join(os.getenv('SCRATCH'), 'SpatialPaperData')
    make_data_archive(archive_dir, runs)
