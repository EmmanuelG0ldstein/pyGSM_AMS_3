
# coding: utf-8

import sys
sys.path.insert(0 ,'/home/cw494069/pyGSM_AMS/')
import numpy as np
from de_gsm import GSM
#from pytc import PyTC
#from orca import Orca
from ams import AMS
from pes import PES
from eigenvector_follow import eigenvector_follow
#from psiw import *
#from rhf_lot import RHF_LOT
from nifty import pvec1d,pmat2d,click,printcool
from molecule import Molecule
import manage_xyz
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
#get_ipython().magic(u'matplotlib inline')

print(sys.path)

##### => Job Data <= #####
charge=0
basis='6-31gs'
filepath1 = 'data/butadiene_ethene.xyz'
filepath2 = 'data/cyclohexene.xyz'



#### => PSIW Obj <= ######
#printcool("Build resources")
#resources = ls.ResourceList.build()
#printcool('{}'.format(resources))

molecule1 = manage_xyz.read_xyz(filepath1, scale = 1)
molecule2 = manage_xyz.read_xyz(filepath2, scale = 1)


####### =>  Build the pyGSM objects <= #########
# level of theory
printcool("Build the pyGSM Level of Theory object (LOT)")
lot1=AMS.from_options(
        states = [(1,0)],
        lot_inp_file = 'ts.in',
        fnm = filepath1
        )
lot2=AMS.from_options(
        states=[(1,0)],
        lot_inp_file = 'ts.in',
        fnm=filepath2
        )

# => Create PES objects <= #
printcool("Building the PES objects")
pes1 = PES.from_options(lot=lot1,ad_idx=0,multiplicity=1)
pes2 = PES.from_options(lot=lot2,ad_idx=0,multiplicity=1)

# => Molecule <= #
printcool("Build the pyGSM Molecule object \n with Translation and Rotation Internal Coordinates (TRIC)")
reactant = Molecule.from_options(fnm=filepath1,PES=pes1,coordinate_type="TRIC")
print(reactant.energy)
product = Molecule.from_options(fnm=filepath2,PES=pes2,coordinate_type="TRIC")

optimizer=eigenvector_follow.from_options(print_level=1)  #default parameters fine here/opt_type will get set by GSM


gsm = GSM.from_options(reactant=reactant,product=product,nnodes=9,optimizer=optimizer,print_level=1)

gsm.go_gsm(rtype = 2, opt_steps = 3)
