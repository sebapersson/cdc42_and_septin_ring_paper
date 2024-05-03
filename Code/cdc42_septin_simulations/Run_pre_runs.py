#!/usr/bin/env python

# For running the pre Cdc42 simulations for the septin ring simulations (to 
# make ring simulations more efficient)
# Input for running is 
# ARGS[1] - experiment to run 
# ARGS[2] - number of processes (running in parallel)
# ARGS[3] - number of reactions (simulations are typically stochastic)
# Simulations should be run on a cluster, and they take typically at least 
# 100h to complete. First run Base_runs_main_sim_fine followed by 
# Smaller_pole_fine. Results from these runs are then read for the septin 
# models.

import Common
import Model_2d
import Solve_steday_state

import multiprocessing as mp
import sys

def number_to_str(number):
    return str(round(number, 1)).replace(".", "P") 

mp.set_start_method("spawn", force=True)

if __name__ == '__main__':
    mp.freeze_support()

    # Sanity check the input arguments
    if len(sys.argv) < 4:
        print("Error: To few input arguments")
        print("Provide three input arguments")
        print("on form experiment, n_threads")
        print("and number of repeats")
        sys.exit(1)

    experiment = sys.argv[1]
    n_threads = int(sys.argv[2])
    n_times_run = int(sys.argv[3])

    if experiment == "Smaller_pole_fine":
        r_list, l_list = [2.0, 2.5, 3.0], [2.0, 2.5, 3.0]
        mod_param_list = []

        # Actin cables parameters 
        for i in range(10):
            for j in range(len(r_list)):
                r, L = r_list[j], l_list[j]
                if r == 2.0:
                    name_read = "Base_fine_R2P0" + str(i)
                    name_save = "Small_pole/Actin" + "R2P0" + str(i)
                elif r == 2.5:
                    name_read = "Base_fine_R2P5" + str(i)
                    name_save = "Small_pole/Actin" + "R2P5" + str(i)
                elif r == 3.0:
                    name_read = "Base_fine_R3P0" + str(i)
                    name_save = "Small_pole/Actin" + "R3P0" + str(i)

                mod_param_list.append(Common.ModelParametersDim(r=r, L=L, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144, k5b=20.8, 
                    read_run=name_read, save_run=name_save))

        sol_opt_list = [Common.SolverOption(term_max_it=350000, print_pwd_it=1000, term_check=False, save_run=True, read_old_run=True)]
        file_loc = Common.FileLocations("Unit_mesh", "2d_pos_save_fine/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d.run_2d_model(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)

    if experiment == "Base_runs_main_sim_fine":
        n_runs = 10
        mod_param_list = []
        for i in range(n_runs):
            i_save = str(i)
            mod_param_list.append(Common.ModelParametersDim(r=2.0, L=2.0, Dm=0.0025, Dm_gdp=0.0025, k2b=0.65, save_run="Base_fine_R2P0" + i_save))
            mod_param_list.append(Common.ModelParametersDim(r=2.5, L=2.5, Dm=0.0025, Dm_gdp=0.0025, k2b=0.65, save_run="Base_fine_R2P5" + i_save))
            mod_param_list.append(Common.ModelParametersDim(r=3.0, L=3.0, Dm=0.0025, Dm_gdp=0.0025, k2b=0.65, save_run="Base_fine_R3P0" + i_save))
            mod_param_list.append(Common.ModelParametersDim(r=4.0, L=4.0, Dm=0.0025, Dm_gdp=0.0025, k2b=0.65, save_run="Base_fine_R4P0" + i_save))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state(mod_param)
        sol_opt_list = [Common.SolverOption(term_max_it=250000, print_pwd_it=10000, term_check=False, save_run=True)]

        file_loc = Common.FileLocations("Unit_mesh", "Base_runs_fine/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d.run_2d_model(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)  

