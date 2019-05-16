import sys
import os 
sys.path.insert(0, '/home/cw494069/pyGSM')

from .de_gsm import *

from .orca import *

lot = Orca.from_options(
        states = [(1,0)], 
        charge = 0, basis = '6-31g(d)', 
        functional = 'B3LYP', 
        #nproc = nproc
        )

filepath = 'examples/tests/benzene.xyz'
geom = manage_xyz.read_xyz(filepath, scale = 1)
e = lot.get_energy(geom, multiplicity = 1, state = 0)
