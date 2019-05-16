from .base_lot import *
import numpy as np
import os
from .units import *
#import scm.plams as plams

class AMS(Lot):
    def run(self, coords, multiplicity):

        tempfilename = 'tempAMSinp.run'
        tempfile = open(tempfilename, 'w')
        if self.lot_inp_file == False:
            return None
        else:
            with open(self.lot_inp_file) as lot_inp:
                inpstring = ''
                inpstring += '#! /bin/sh \n \n'
                inpstring += '$ADFBIN/ams << eor \n \n'
                tempfile.write(inpstring)
                lot_inp_lines = lot_inp.readlines()
            for line in lot_inp_lines:
                tempfile.write(line)
            inpstring = 'System \n    Atoms \n        '
            for coord in coords:
                for i in coord:
                    inpstring += str(i) + ' '
                inpstring += '\n        '
            inpstring += 'End \nEnd \n\neor'
        tempfile.write(inpstring)
        tempfile.close()
        cmd = 'chmod u+x {}; ./{}'.format(tempfilename, tempfilename)
        os.system(cmd)

        #Energie auslesen
        #rkfpath = 'ams.results/ams.rkf'
        #rkffile = plams.KFFile(rkfpath)
        #engine = rkffile[('General', 'engine')]
        #enginefile = plams.KFFile('{}.rkf'.format(engine))
        #self.E.append(multiplicity, enginefile[('AMSResults', 'Energy')])

        ##Grad auslesen
        #self.grada.append(multiplicity, enginefile[('AMSResults', 'Gradients')])
        return

    def get_energy(self,coords,multiplicity,state):
        #if self.has_nelectrons==False:
        #    for i in self.states:
        #        self.get_nelec(geom,i[0])
        #    self.has_nelectrons==True
        if self.hasRanForCurrentCoords==False or (coords != self.currentCoords).any():
            self.currentCoords = coords.copy()
            geom = manage_xyz.np_to_xyz(self.geom,self.currentCoords)
            self.runall(geom)
            self.hasRanForCurrentCoords=True
        print((self.E))
        print(multiplicity)
        tmp = self.search_tuple(self.E,multiplicity)
        return np.asarray(tmp[state][1])*KCAL_MOL_PER_AU

    def get_gradient(self,coords,multiplicity,state):
        #if self.has_nelectrons==False:
        #    for i in self.states:
        #        self.get_nelec(geom,i[0])
        #    self.has_nelectrons==True
        if self.hasRanForCurrentCoords==False or (coords != self.currentCoords).any():
            self.currentCoords = coords.copy()
            geom = manage_xyz.np_to_xyz(self.geom,self.currentCoords)
            self.runall(geom)
        tmp = self.search_tuple(self.grada,multiplicity)
        return np.asarray(tmp[state][1])*ANGSTROM_TO_AU
