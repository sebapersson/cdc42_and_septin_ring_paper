#!/usr/bin/env python


'''
This file runs the reaction diffusion experiments. The experiment to run 
is provided by the user as a command line argument. Also, the user can 
provide the number of threads to use, and the number of times to repeat 
each condition (optional argument). 
Args (note order is important):
    experiment (str), the experiment to run 
    n_threads (int), the number of threads to use
    n_times_run = 5 (int), the number of times to run each condition 

The available experiments are:
    find_mesh_zero_bud_scar, finding the mesh for zero bud-scars 
    find_d_zero_bud_scar, finding d for no bud-scars 
    no_bud_scars, the experimental design not using bud-scars 
'''

import sys
import FEM
#import Solve_rd_system as srd
from Solve_rd_system import *
from Model_parameters import ModelParametersWithDim
import Model_parameters

import multiprocessing as mp
mp.set_start_method("spawn", force=True)

if __name__ == '__main__':
    mp.freeze_support()
    
    # Sanity check the input arguments
    if len(sys.argv) - 1 < 2:
        print("Error: To few input arguments")
        print("Provide two or more inputs")
        sys.exit(1)
    elif len(sys.argv) - 1 > 3:
        print("Error: To many input arguments")
        print("Provide a maximum of three arguments")
        sys.exit(1)
    
    
    # Read the arguments 
    experiment = sys.argv[1]
    n_threads = int(sys.argv[2])
    if len(sys.argv) == 4:
        n_times_run = int(sys.argv[3])
    else:
        n_times_run = 10
    
    if experiment == "validate_pole_size":
        bud_scar_list = ["zero"]
        r_list = np.linspace(3.0, 4.35, num=20)
        param_list = []
        for r in r_list:
            param_list.append(Model_parameters.ModelParametersWithDim(d=110, k3=1*1.0, R = r))
        print("Length of param_list = {}".format(len(param_list)))
        run_rd_simulations(param_list, bud_scar_list, times_run=n_times_run, delete_old_pvd=True,
                           n_threads=n_threads, nucleus=False, tag_exp="Pole_size_pos/",
                           fem_opt_list = [FEM.FemOpt(check_termination=False, max_it=30000)])
        sys.exit(0)

