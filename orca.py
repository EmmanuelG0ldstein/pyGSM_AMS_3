from .base_lot import *
import numpy as np
import os
from .units import *

class Orca(Lot):
    
    def run(self,coords,multiplicity):
        
        if self.lot_inp_file == False:
            inpstring = '!'
            inpstring += ' '+self.functional
            inpstring += ' '+self.basis
            inpstring += ' EnGrad\n\n'# SOSCF SlowConv \n\n'
            inpstring += '%scf\nMaxIter 300\nconvergence strong\n sthresh 1e-7\n'
            inpstring += 'thresh 1e-11\n tcut 1e-13 \n directresetfreq 1 \n SOSCFStart 0.00033\nend\n'
            #inpstring += '%scf\nMaxIter 300\nend\n'
            inpstring += '\n%maxcore 1000\n\n'
            inpstring += '%pal\nnproc {}\nend\n\n'.format(self.nproc)
        else:
            inpstring = ''
            with open(self.lot_inp_file) as lot_inp:
                lot_inp_lines = lot_inp.readlines()
            for line in lot_inp_lines:
                inpstring += line

        inpstring += '\n* xyz {} {}\n'.format(self.charge,multiplicity)
        for coord in coords:
            for i in coord:
                inpstring += str(i)+' '
            inpstring += '\n'
        inpstring += '*'
        print(inpstring)
        tempfilename = 'tempORCAinp_{}'.format(multiplicity)
        tempfile = open(tempfilename,'w')
        tempfile.write(inpstring)
        tempfile.close()

        testfilename = 'test'
        testfile = open(testfilename, 'w')
        testfile.write(inpstring)
        testfile.close()
        
#        path2orca = os.popen('which orca').read().rstrip()
        path2orca = '/home/cw494069/Orca/orca_4_0_1_2_linux_x86-64_openmpi202/orca'
        user = os.environ['USER']
        cwd = os.environ['PWD']
        try:
            slurmID = os.environ['SLURM_ARRAY_JOB_ID']
            try:
                slurmTASK = os.environ['SLURM_ARRAY_TASK_ID']
                runscr = '/tmp/'+user+'/'+slurmID+'/'+slurmTASK
            except:
                runscr = '/tmp/'+user+'/'+slurmID
        except:
            orcascr = 'temporcarun'
            #runscr = '/tmp/'+user+'/'+orcascr
            #runscr = '/home/cw494069/tmp/' + user + '/' + orcascr
            runscr = cwd #+ '/' + orcascr 
        os.system('mkdir -p {}'.format(runscr))
        os.system('mv {} {}/'.format(tempfilename,runscr))
        #cmd = 'cd {}; {} {} > {}/{}; cd {}'.format(runscr,path2orca,tempfilename,runscr,tempfilename,cwd)
        cmd = '~/Orca/orca_4_0_1_2_linux_x86-64_openmpi202/orca test'
        print(cmd)
        print((os.system('cd {};ls -l'.format(runscr))))
        print((os.system('cd {}; cat {}'.format(runscr, tempfilename))))
        os.system(cmd)

        #engradpath = runscr+'/{}.engrad'.format(tempfilename) 
        engradpath = runscr+'/{}.engrad'.format(testfilename) 
        with open(engradpath) as engradfile:
            engradlines = engradfile.readlines()
            print(engradlines)
        temp = 100000
        for i,lines in enumerate(engradlines):
            if 'current total energy' in lines:
                temp = i
            if i > temp+1:
                self.E.append((multiplicity,float(lines.split()[0]))) 
                break

        temp = 100000
        tmp = []
        tmp2 = []
        for i,lines in enumerate(engradlines):
            if 'current gradient' in lines:
                temp = i
            if  i> temp+1:
                if "#" in lines:
                    break
                tmp2.append(float(lines.split()[0]))
            if len(tmp2) == 3:
                tmp.append(tmp2)
                tmp2 = []
        self.grada.append((multiplicity,tmp))
        return

    def get_energy(self,coords,multiplicity,state):
        #if self.has_nelectrons==False:
        #    for i in self.states:
        #        self.get_nelec(coords,i[0])
        #    self.has_nelectrons==True
        if self.hasRanForCurrentCoords==False or (coords != self.currentCoords).any():
            self.currentCoords = coords.copy()
            geom = manage_xyz.np_to_xyz(self.geom,self.currentCoords)
            self.runall(geom)
        tmp = self.search_tuple(self.E,multiplicity)
        print((type(self.E)))
        print(tmp)
        return np.asarray(tmp[state][1])*KCAL_MOL_PER_AU

    def get_gradient(self,coords,multiplicity,state):
        #if self.has_nelectrons==False:
        #    for i in self.states:
        #        self.get_nelec(coords,i[0])
        #    self.has_nelectrons==True
        if self.hasRanForCurrentCoords==False or (coords != self.currentCoords).any():
            self.currentCoords = coords.copy()
            geom = manage_xyz.np_to_xyz(self.geom,self.currentCoords)
            self.runall(geom)
        tmp = self.search_tuple(self.grada,multiplicity)
        return np.asarray(tmp[state][1])#*ANGSTROM_TO_AU #ORCA grad is given in AU
