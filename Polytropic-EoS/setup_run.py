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
def get_input_array(n,N1,bhigh, blow, Nn):
    input_data  = []
    delta = (bhigh-blow)/(Nn)
    for j in range(Nn+1):
        input_data.append([n,N1, blow+ j*delta])

    input_data = np.array(input_data)
    return input_data
def get_input_array_from_b(n,N1,barr):
    input_data  = []
    for j in range(len(barr)):
        input_data.append([n,N1, barr[j]])
    input_data = np.array(input_data)
    return input_data
#===============================================================================
input_data = get_input_array_from_b(1.0,0.75, [0.0054320988,
0.1358024691,0.2661728395,0.347654321])
#===============================================================================
sim = Sim()
sim.output_bvals = 0
sim.solve_tides = 1
sim.solve_visc=0
sim.write_vec=0
sim.xil = 0.05
sim.N_interp_LE = 5000
sim.N_write = 99
sim.steps_num = 5
sim.factor_match = 0.8
#----------------------
sim.Omegalow = 0.06
sim.Omegahigh = 0.16
sim.div_Omega = 1000
sim.out_dir = out_path + "/Runs"
sim.animscript = home_path + "/plot-vals.ipynb"
sim.setup_script = home_path + "/setup_run.py"
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


