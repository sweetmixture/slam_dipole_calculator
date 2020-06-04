#!/bin/python

import sys,math
import numpy as np
	
class slam_dipole_mod:


    def __init__(self,*args):

        self.pos_integral = float(args[0])                      # position integral val
        self.config = args[1]                                   # main configuration
        self.mo     = args[2]                                   # molecular orbital info
        self.type   = args[3]                                   # type info

        self.mm_cnt = int(args[4])                              # MM_CNT
        self.qm_cnt = int(args[5])                              # QM_CNT

        self.mm_config = [[ 0 for i in range(6)] for j in range(self.mm_cnt)]   # 1,2,3,4,5,6 - name type charge x y z
        self.qm_config = [[ 0 for i in range(6)] for j in range(self.qm_cnt)]   # 1,2,3,4     - name core_q sp_q x y z
        self.mo_config = [[ 0 for i in range(6)] for j in range(self.qm_cnt)]   # 1,2,3,4     - s, px, py, pz 

        try:
            with open(self.config,"r") as c:			# READ CONFIG INFO

                rl = c.readline()

                for i in range(self.mm_cnt):			# READ MM CONFIG
								# LIST ROW STRIDE 6 ... NAME / TYPE(c/s) / CHARGE / x / y / z
                    rl = c.readline()
                    spl= rl.split()

                    self.mm_config[i][0] = spl[0]
                    self.mm_config[i][1] = spl[1]
                    #self.mm_config[i][2] = 
                    self.mm_config[i][3] = float(spl[2])
                    self.mm_config[i][4] = float(spl[3])
                    self.mm_config[i][5] = float(spl[4])

                for i in range(self.qm_cnt):

                    rl = c.readline()
                    spl= rl.split()
                    self.qm_config[i][0] = spl[0]
                    self.qm_config[i][3] = float(spl[1])
                    self.qm_config[i][4] = float(spl[2])
                    self.qm_config[i][5] = float(spl[3])

        except FileNotFoundError:
            print(FileNotFoundError)

        try:							# READ TYPE INFO ... SAVE CHARGE INFO
            with open(self.type,"r") as t:

                rl = t.readline()
        
                for i in range(self.mm_cnt):

                    rl = t.readline()
                    spl= rl.split()

                    if spl[0] == self.mm_config[i][0] and spl[1] == self.mm_config[i][1]:
                        self.mm_config[i][2] = float(spl[2])

                for i in range(self.qm_cnt):

                    rl = t.readline()
                    spl= rl.split()

                    if spl[0] == self.qm_config[i][0]:

                        self.qm_config[i][1] = float(spl[1])
                        self.qm_config[i][2] = float(spl[2])

        except FileNotFoundError:
            print(FileNotFoundError)

        try:							# READ MO INFO ... 
            with open(self.mo,"r") as m:
								# LIST ROW STRIDE 6 ... NAME / ENERGY(eV) / s / px / py / pz
                rl = m.readline()

                for i in range(self.qm_cnt):

                    rl = m.readline()
                    spl= rl.split()
    
                    self.mo_config[i][0] = spl[0]           # mo atom name
                    self.mo_config[i][1] = float(spl[1])    # mo atom energy
                    self.mo_config[i][2] = float(spl[2])    # mo atom s_
                    self.mo_config[i][3] = float(spl[3])    # mo atom px
                    self.mo_config[i][4] = float(spl[4])    # mo atom py
                    self.mo_config[i][5] = float(spl[5])    # mo atom pz

        except FileNotFoundError:
            print(FileNotFoundError)

    def get_mm_dipole(self):

        mm_dip = [[ 0. for i in range(3) ] for j in range(self.mm_cnt)]

        for i in range(self.mm_cnt):	# GET MM DIPOLE ... VALUES MAY CHANGE DEPENDING ON THE GIVEN REFERENCE FRAME ( IF CHARGE IS NOT NEUTRAL )
            
            mm_dip[i][0] = self.mm_config[i][3]*self.mm_config[i][2]          # get dip x elem
            mm_dip[i][1] = self.mm_config[i][4]*self.mm_config[i][2]          # get dip y elem
            mm_dip[i][2] = self.mm_config[i][5]*self.mm_config[i][2]          # get dip z elem

        return mm_dip

    def get_qm_dipole(self):		# GET QM DIPOLE ... VALUES MAY CHANGE DEPENDING ON THE GIVEN REFERENCE FRAME ( IF CHARGE IS NOT NEUTRAL )

        qm_dip = [[ 0. for i in range(3) ] for j in range(self.qm_cnt)]
    
        for i in range(self.qm_cnt):
	#   dip_elem =  r*(qc + qs) + <psi| r | psi>*qs
            qm_dip[i][0]=self.qm_config[i][3]*(self.qm_config[i][1] + self.qm_config[i][2]) + 2.*self.mo_config[i][2]*self.mo_config[i][3]*self.pos_integral*self.qm_config[i][2]
            qm_dip[i][1]=self.qm_config[i][4]*(self.qm_config[i][1] + self.qm_config[i][2]) + 2.*self.mo_config[i][2]*self.mo_config[i][4]*self.pos_integral*self.qm_config[i][2]
            qm_dip[i][2]=self.qm_config[i][5]*(self.qm_config[i][1] + self.qm_config[i][2]) + 2.*self.mo_config[i][2]*self.mo_config[i][5]*self.pos_integral*self.qm_config[i][2]

        return qm_dip

    def get_cluster_dipole(self):	# SUM OF THE CLUSTER DIPOLE MOMENT

        self.mm_dip = self.get_mm_dipole()
        self.qm_dip = self.get_qm_dipole()
    
        self.cluster_dip = [ 0. for i in range(3) ]

        for i in range(self.mm_cnt):
            self.cluster_dip[0] += self.mm_dip[i][0]
            self.cluster_dip[1] += self.mm_dip[i][1]
            self.cluster_dip[2] += self.mm_dip[i][2]
        for i in range(self.qm_cnt):
            self.cluster_dip[0] += self.qm_dip[i][0]
            self.cluster_dip[1] += self.qm_dip[i][1]
            self.cluster_dip[2] += self.qm_dip[i][2]


    def write_intro(self):

	print("#"*90)
	print("")
	print(" Sp-Lone pair involved Atomistic Model ( S L A M ) ")
	print("")
	print(" Dipole Moment Calculations Toolkit ")
	print("")
	print(" Compatibility	: SLAM Version 2.2_snapshot ")
	print("")
	print(" Last edited	: 04 - 06 - 2020 ")
	print("")
	print("#"*90)
	print("")


    def write(self):			# WRITE GENERAL OUTPUT

	self.write_intro()
	print(" Configuration / Dipole Moment Elements ")
	print("")
	print("-"*90)
	print(" Species.  x           y           z               u_x           u_y           u_z")
	print("-"*90)

        for i in range(self.mm_cnt):

            if self.mm_config[i][1] == "c" and self.mm_config[i+1][1] == "s":   # if the type is core & it has shell
                atom_name   = self.mm_config[i][0]
                pos_x       = self.mm_config[i][3]
                pos_y       = self.mm_config[i][4]
                pos_z       = self.mm_config[i][5]
            
                dip_x       = self.mm_dip[i+1][0] - self.mm_dip[i][0]           # Following standard convention ... r(-) - r(+)
                dip_y       = self.mm_dip[i+1][1] - self.mm_dip[i][1]
                dip_z       = self.mm_dip[i+1][2] - self.mm_dip[i][2]

                #print("%3s%12.6f%12.6f%12.6f%18.6e%14.6e%14.6e" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))
                print("%3s%12.6f%12.6f%12.6f%18.6f%14.6f%14.6f" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))

            elif self.mm_config[i][1] == "c" and self.mm_config[i+1][1] != "s": # if the type is core and it does not have shell
                atom_name   = self.mm_config[i][0]
                pos_x       = self.mm_config[i][3]
                pos_y       = self.mm_config[i][4]
                pos_z       = self.mm_config[i][5]
            
                dip_x       = 0.
                dip_y       = 0.
                dip_z       = 0.
            
                #print("%3s%12.6f%12.6f%12.6f%18.6e%14.6e%14.6e" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))
                print("%3s%12.6f%12.6f%12.6f%18.6f%14.6f%14.6f" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))
        for i in range(self.qm_cnt):
            atom_name       = self.qm_config[i][0]
            pos_x           = self.qm_config[i][3]
            pos_y           = self.qm_config[i][4]
            pos_z           = self.qm_config[i][5]

            dip_x           = self.qm_dip[i][0]
            dip_y           = self.qm_dip[i][1]
            dip_z           = self.qm_dip[i][2]

            #print("%3s%12.6f%12.6f%12.6f%18.6e%14.6e%14.6e" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))
            print("%3s%12.6f%12.6f%12.6f%18.6f%14.6f%14.6f" % (atom_name,pos_x,pos_y,pos_z,dip_x,dip_y,dip_z))
	print("-"*90)
        print("")
        print(" Total Cluster Dipole (e*Angstrom) :%12.6f%12.6f%12.6f" % (self.cluster_dip[0], self.cluster_dip[1], self.cluster_dip[2]))
        print("")

	print("#"*90)
	print("%45s" % ("Finalising"))
	print("#"*90)

if __name__=='__main__':

    dip_inst = slam_dipole_mod(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
    dip_inst.get_cluster_dipole()
    dip_inst.write()
