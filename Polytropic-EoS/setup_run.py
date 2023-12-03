#!/usr/bin/env python
import sys,subprocess
import os
#===============================================================================
home_path = os.getcwd()
out_path = home_path+ "/output"
#===============================================================================
from sim_class import Sim
import numpy as np
from multiprocessing import Pool
import time
from datetime import datetime
def get_time():
    current_time = datetime.now()
    timestr = current_time.strftime("%a")+"_"+current_time.strftime("%b")+"_"  \
    + str(current_time.day) +"_"+ str(current_time.hour) \
    + "_"+str(current_time.minute)+"_"+str(current_time.second)

    return timestr

def get_xi_out(n, b, Omega, val):
    ans = val/Omega
    ans /= pow(n+1,0.5)*pow(b,0.5)
    return ans

#===============================================================================
# input_data  = np.array([[0.5,0.0162962963],[0.5,0.0651851852], [0.5,0.4074074074], [0.5,0.7985185185], [0.5,1.0429629630], [0.5,1.3200000000],
# [0.75,0.0092592593],[0.75,0.0370370370], [0.75,0.1481481481], [0.75,0.4537037037], [0.75,0.5925925926], [0.75,0.7500000000]])

# input_data  = np.array([ [0.5,0.5,1.0429629630], [0.75,0.75,0.5925925926]])


input_data = np.array([[1.0,0.75,0.44]])
#===============================================================================
sim = Sim()
sim.output_bvals = 0
# sim.solve_visc=1
sim.xil = 0.05
sim.N_interp_LE = 5000
sim.N_write = 2999
sim.steps_num = 3
sim.factor_match = 0.7
#----------------------
sim.Omegalow = 1e-3
sim.Omegahigh = 0.5
sim.div_Omega = 50
sim.out_dir = out_path + "/Runs"
sim.animscript = home_path + "/plot-vals.ipynb"
#----------------------
def launch_sim(vals):
    sim.n = vals[0]
    sim.N1 = vals[1]
    sim.b = vals[2]
    sim.launch()


#===============================================================================
if __name__ == '__main__':
    if len(input_data) >=6:
        pool_nums = 6
    else :
        pool_nums = len(input_data)


    print("pool_nums = ", pool_nums)

    t_start = time.time()
    print("Starting multiprocessing pool..")
    print("Data saved at:{}".format(sim.out_dir))
    pool = Pool(pool_nums)
    result = pool.map_async(launch_sim, input_data)
    pool.close()

    while True:
        if not result.ready():
            print('We\'re not done yet, %s tasks to go!' % result._number_left)
            time.sleep(2)
        else:
            break

    pool.join()

    t_end = time.time()
    print("Finished process. \nTime = ",t_end-t_start," s")
    print("Data saved at:{}".format(sim.out_dir))
# #===============================================================================


