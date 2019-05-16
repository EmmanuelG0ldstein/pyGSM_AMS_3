import scm.plams as plams

import os 
import sys
print(sys.path)

rkffile = plams.KFFile('ams.results/ams.rkf')
#print(rkffile.get_skeleton())
#print(rkffile.__iter__())
resultpath = rkffile[('EngineResults', 'Files(1)')]
#print(rkffile[('History','maxStep(95)')])

print(resultpath)

resultfile = plams.KFFile(resultpath)
#print(resultfile.get_skeleton())
print(resultfile[('AMSResults', 'Energy')])
print(resultfile[('AMSResults', 'Gradients')])

