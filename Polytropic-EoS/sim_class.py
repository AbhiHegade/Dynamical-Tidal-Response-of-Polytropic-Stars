#!/usr/bin/env python
# coding: utf-8
from datetime import datetime
import os, sys, time, shutil, subprocess
import numpy as np

class Sim:
    def __init__(self):
        self.home_dir = str(os.getcwd())
    def get_time(self):
        current_time = datetime.now()
        # timestr = current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"  \
        # + str(current_time.day) +"_"+ str(current_time.hour) \
        # + "_"+str(current_time.minute)+"_"+str(current_time.second) \
        # + "_n_" + str(self.n) + "_b_" + str(self.b)

        timestr = "n_" + str(self.n) + "_N_" + str(self.N1) \
        + "_" + "_b_" + str(self.b) \
        + current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"  \
        + str(current_time.day) +"_"+ str(current_time.hour) \
        + "_"+str(current_time.minute)+"_"+str(current_time.second) 

        # timestr = "n_" + str(self.n) + "_b_" + str(self.b) 

        return timestr
#===============================================================================
    def make_output_dir(self,level=-1):
        timestr = self.get_time()
        if level == -1:
            self.output_dir = self.out_dir + "/" + timestr 
            #     + \
            # "_{}_".format(self.ic)+ "_{}_".format(self.type)+ "_{}_".format(self.nr) + "_{}_".format(self.rhoc) + \
            #  "_{}_".format(self.tolnum) + "_{}_".format(self.thrnum) + "_{}_".format(self.cfl) + "_{}_".format(self.diss)

            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
        else:
            self.output_dir = self.convergence_dir +"/{}".format(level)
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
#===============================================================================
    def write_sim_params(self):
        with open(self.output_dir+'/sim_params.txt','w') as f:
            attrs= vars(self)
            for param in attrs:
                if param == "b_vals":
                    pass
                else:
                    f.write('{} {}\n'.format(param,attrs[param]))
#===============================================================================
    def make_output_file(self):
        self.output_file= self.output_dir+'/'+'output.out'
        with open(self.output_file, 'w') as f:
            pass
    def copy_anim_script(self,level=0):
        if(level ==0):
            subprocess.call('cp {} {}/'.format(self.animscript,self.output_dir), shell=True)
        else:
            subprocess.call('cp {} {}/'.format(self.animscript_adv,self.output_dir), shell=True)
#===============================================================================
    def output_b_vals(self):
        if self.output_bvals:
            np.savetxt(self.output_dir+ "/bvals.txt",self.b_vals)
        else:
            pass
#===============================================================================
    def launch(self,level=-1):
        self.make_output_dir()
        self.make_output_file()
        self.copy_anim_script(0)
        # self.copy_anim_script(1)
        self.write_sim_params()
        self.output_b_vals()

        subprocess.call('cp {} {}/'.format(self.setup_script,self.out_dir), shell=True)
        subprocess.call('\n./bin/default.run {} > {}/output.out'.format(self.output_dir,self.output_dir), shell=True)
   
#===============================================================================
