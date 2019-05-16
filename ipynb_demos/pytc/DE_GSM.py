
# coding: utf-8

# In[3]:

import sys
sys.path.insert(0,'/home/caldaz/module/pyGSM')
import numpy as np
from de_gsm import GSM
from pytc import PyTC
from pes import PES
from eigenvector_follow import eigenvector_follow
from psiw import *
from rhf_lot import RHF_LOT
from nifty import pvec1d,pmat2d,click,printcool
from molecule import Molecule
import manage_xyz
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
get_ipython().magic(u'matplotlib inline')


# In[4]:

##### => Job Data <= #####
charge=0
basis='6-31gs'
filepath1 = 'data/butadiene_ethene.xyz'
filepath2 = 'data/cyclohexene.xyz'


# In[5]:

#### => PSIW Obj <= ######
printcool("Build resources")
resources = ls.ResourceList.build()
printcool('{}'.format(resources))

molecule1 = ls.Molecule.from_xyz_file(filepath1)
molecule2 = ls.Molecule.from_xyz_file(filepath2)

geom1 = geometry.Geometry.build(
    resources=resources,
    molecule=molecule1,
    basisname=basis,
    )
geom2 = geometry.Geometry.build(
    resources=resources,
    molecule=molecule2,
    basisname=basis,
    )

printcool('{}'.format(geom1))
printcool('{}'.format(geom2))

ref1 = RHF(RHF.default_options().set_values({
    'geometry' : geom1,
    'dft_functional' : 'B3LYP',
    'dft_grid_name' : 'SG0',
    }))
ref2 = RHF(RHF.default_options().set_values({
    'geometry' : geom2,
    'dft_functional' : 'B3LYP',
    'dft_grid_name' : 'SG0',
    }))

psiw1 = RHF_LOT.from_options(rhf=ref1)
psiw2 = RHF_LOT.from_options(rhf=ref2)


# In[6]:

####### =>  Build the pyGSM objects <= #########
# level of theory
printcool("Build the pyGSM Level of Theory object (LOT)")
lot1=PyTC.from_options(states=[(1,0)],psiw=psiw1,fnm=filepath1)
lot2=PyTC.from_options(states=[(1,0)],psiw=psiw2,fnm=filepath2)

# => Create PES objects <= #
printcool("Building the PES objects")
pes1 = PES.from_options(lot=lot1,ad_idx=0,multiplicity=1)
pes2 = PES.from_options(lot=lot2,ad_idx=0,multiplicity=1)

# => Molecule <= #
printcool("Build the pyGSM Molecule object \n with Translation and Rotation Internal Coordinates (TRIC)")
reactant = Molecule.from_options(fnm=filepath1,PES=pes1,coordinate_type="TRIC")
product = Molecule.from_options(fnm=filepath2,PES=pes2,coordinate_type="TRIC")

optimizer=eigenvector_follow.from_options(print_level=1)  #default parameters fine here/opt_type will get set by GSM


# In[7]:

printcool("Primitives before union")
print(reactant.coord_obj.Prims.Internals)


# In[8]:

printcool("Primitives before union")
print(product.coord_obj.Prims.Internals)


# In[9]:

printcool("Topology before union")
plt.plot()
nx.draw(reactant.coord_obj.Prims.topology,with_labels=True,font_weight='bold') 
plt.show()


# In[10]:

printcool("Topology before union")
plt.plot()
nx.draw(product.coord_obj.Prims.topology,with_labels=True,font_weight='bold') 
plt.show()


# In[11]:

printcool("Wilson B-Matrix (dq_i/dx_j) before union")
Bmatp = reactant.coord_obj.Prims.wilsonB(reactant.xyz)
plt.imshow(Bmatp, cmap=plt.cm.get_cmap('RdBu'))
plt.show()


# In[12]:

printcool("Coordinate Basis before union")
plt.imshow(reactant.coord_basis, cmap=plt.cm.get_cmap('RdBu'))
plt.show()


# In[13]:

gsm = GSM.from_options(reactant=reactant,product=product,nnodes=9,optimizer=optimizer,print_level=1)


# In[14]:

printcool("Topology after union")
plt.plot()
nx.draw(gsm.nodes[0].coord_obj.Prims.topology,with_labels=True,font_weight='bold') 
plt.show()


# In[15]:

printcool("Wilson B-Matrix (dq_i/dx_j)")
Bmatp = gsm.nodes[0].coord_obj.Prims.wilsonB(reactant.xyz)
plt.imshow(Bmatp, cmap=plt.cm.get_cmap('RdBu'))
plt.show()


# In[16]:

printcool("Coordinate Basis")
plt.imshow(gsm.nodes[0].coord_basis, cmap=plt.cm.get_cmap('RdBu'))
plt.show()


# In[ ]:



