import sys
import os
import os.path
import subprocess

if __name__ == "__main__":

    os.environ['GRIDDESC'] = os.path.join(os.environ['HOME'],
                                          'Data',
                                          'STEM',
                                          'input',
                                          'GRIDDESC.txt')
    os.environ['WRF_sfc_met'] = os.path.join(
        os.environ['SARIKA_INPUT'],
        'meteo2d-124x124-18levs-2008-2009.nc')
    os.environ['kettle_fsoil'] = os.path.join(
        os.environ['SARIKA_INPUT'],
        'surfem-124x124-kettle-soil-cos_2008_2009.nc')
    os.environ['crop_pct'] = os.path.join(
        os.environ['HOME'],
        'projects', 'COS (ecampbell3)', 'Fractional_US_Cropland',
        'Ramankutty_etal_Cropland2000_pct_124x124_IOAPI.nc')
    outfile = './whelan_kettle_hybrid_fsoil_124x124.nc'
    os.environ['fsoil_out'] = outfile
    if os.path.exists(outfile):
        sys.stdout.write('removing {}'.format(outfile))
        sys.stdout.flush()
        os.remove(outfile)
    subprocess.call('make -f meld_whelan_kettle_fsoils.mk', shell=True)
    subprocess.call('./meld_whelan_kettle_fsoils.x > log.log', shell=True)
