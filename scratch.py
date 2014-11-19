import socket
import map_grid
import os, os.path

if 'Timothys-MacBook-Air.local' in socket.gethostname():
    aqout_data = (os.path.join(os.getenv('HOME'), 'work', 'Data',
                               'STEM', 'aq_out_data.cpickle'))
else:
    aqout_data = os.path.join(os.getenv('HOME'), 'Data', 'STEM',
                              'aq_out_data.cpickle')
cos_conc, gpp, fCOS = map_grid.assemble_data(aqout_data)
