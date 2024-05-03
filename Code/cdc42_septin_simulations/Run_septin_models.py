#!/usr/bin/env python

# For running the pre Cdc42 simulations for the septin ring simulations
# were explored there are many options, but the name should be clear. For 
# example experiment == "Axl2_v4_exo_nf_diluted_septin" is the Axl2_v4 
# model (SBER in MS) with NF conditions, and lower septin concentration.
# Input for running is 
# ARGS[1] - experiment to run 
# ARGS[2] - number of processes (running in parallel)
# ARGS[3] - number of reactions (simulations are typically stochastic)
# Simulations should be run on a cluster, and they take typically at least 
# 100h to complete.

import Common 
import Solve_steday_state
import Exo_endo_cytosis
import Axl2_model_v1
import Axl2_model_v4

import multiprocessing as mp
import numpy as np 
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

    if len(sys.argv) == 5:
        cell_size = sys.argv[4]
    else:
        cell_size = "standard"


    if experiment == "Axl2_v4_exo_nf_base":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5, 2.0, 2.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for k13 in k13_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole_medium_GAP/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole_weak_GAP/Actin" + "R2P50"
                                    elif r == 3.0:
                                        read_run = "Small_pole_weak_GAP/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole_medium_GAP/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.35, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=k13, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=0.00045))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=20000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_nf_cond/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)

    if experiment == "Axl2_v4_exo_nf_diluted_septin":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        Stot_list = [0.15, 0.10]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for S_tot in Stot_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole_medium_GAP/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole_weak_GAP/Actin" + "R2P50"
                                    elif r == 3.0:
                                        read_run = "Small_pole_weak_GAP/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole_medium_GAP/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.35, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=0.00045, S_tot=S_tot))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=20000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_nf_cond_Stot/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)


    if experiment == "Axl2_v4_diluted_septin":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        Stot_list = [0.15, 0.10]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for S_tot in Stot_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole_medium_GAP/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P50"
                                    elif r == 3.0:
                                        read_run = "Small_pole_weak_GAP/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole_medium_GAP/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=0.00045, S_tot=S_tot))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=20000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_exo_Stot/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)

    if experiment == "Axl2_v4_k20_strong":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.4, 0.5, 0.6]
        k22_list = [10.5]
        k23_list = [26.0]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for i in range(len(k23_list)):
                        for k24 in k24_list:
                            for k25 in k25_list:
                                if r == 2.0:
                                    read_run = "Small_pole_medium_GAP/Actin" + "R2P0"
                                elif r == 2.5:
                                    read_run = "Small_pole/Actin" + "R2P50"
                                elif r == 3.0:
                                    read_run = "Small_pole_weak_GAP/Actin" + "R3P0"
                                elif r == 4.0:
                                    read_run = "Small_pole_medium_GAP/Actin" + "R3P0"
                                        
                                mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=0.00045))
                                exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=250000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=20000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_k20/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)                        

    if experiment == "Axl2_v4_exo_stronger_GAPs":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [2.0, 2.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for k13 in k13_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole_medium_GAP/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P50"
                                    elif r == 3.0:
                                        read_run = "Small_pole_weak_GAP/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole_medium_GAP/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=k13, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=0.00045))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=20000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_stronger_GapS/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)        

    if experiment == "Axl2_v4_exo_nf_normal_base":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P5"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=300000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=20000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_normal_cond/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)

    if experiment == "Axl2_v4_cdc24kk":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k15_list = [1.0, 1.0*0.75]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5, 4.5*0.8]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k19 in k19_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P5"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=k19, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=0.1, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=250000, print_pwd_it=20000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=15000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_cdc24kk/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)

    if experiment == "Axl2_v1_noexo":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k13 in k13_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P5"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V1(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=k13, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 9,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=1/40000, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                         gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=2000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=15000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v1/No_exo/")

        Axl2_model_v1.run_axl2_v1(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)        

    if experiment == "Axl2_v4_cdc24kk_sat_rec":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k15_list = [1.0, 1.0*0.75]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5, 4.5*0.8]
        DmFOO_list = [0.00045]
        GapSlim_list = [1.5, 1.4]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for GapSlim in GapSlim_list:
                        for i in range(len(k23_list)):
                            for k19 in k19_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P5"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=k19, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=0.1, k25=k25, DmFOO=0.00045, GapSlim=GapSlim))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=20000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=15000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_cdc24kk_GapS_limit/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)        
        
    # Fun experiment showing that for subset of parameters we can collapse the ring
    if experiment == "Axl2_v4_exo_collapse_ring":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k2b_list = [0.63*1.2, 0.63*1.5]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k2b in k2b_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P5"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=k2b, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=0.1, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=100000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=20000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_collapse_ring/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)                

    if experiment == "Axl2_v4_faster_with_FOO":
        
        r_list = [2.5]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [0.0, 5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P5"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=2000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=5000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/Test_faster_with_FOO/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)        

    if experiment == "Axl2_v4_get_nodes_accepted":
        
        r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        eps_list = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for eps_hit in eps_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P5"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=False,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=eps_hit, only_print_eps_hit=True))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=2000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_nodes_accepted/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)
        
    if experiment == "Axl2_v4_exo_wide1":

        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P51"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                         gamma_P=1.0, cutoff_p_move=100, eps_hit=0.3))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide1" + cell_size + "/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)

    if experiment == "Axl2_v4_exo_wide2":
        
        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P51"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                         gamma_P=1.0, cutoff_p_move=100, eps_hit=0.4))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide2" + cell_size + "/") 

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)        

    if experiment == "Axl2_v4_exo_wide3":
        
        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P51"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.5))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide3" + cell_size + "/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

    if experiment == "Axl2_v4_exo_wide4":
        
        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P51"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.6))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide4" + cell_size + "/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

    if experiment == "Axl2_v4_exo_wide5":
        
        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P51"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.7))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide5" + cell_size + "/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)        

    if experiment == "Axl2_v4_exo_wide6":
        
        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P51"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.8))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide6" + cell_size + "/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)        

    if experiment == "Axl2_v4_exo_wide7":
        
        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P51"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.9))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide7" + cell_size + "/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)        

    if experiment == "Axl2_v4_exo_wide8":
        
        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P51"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.99))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=5000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide8" + cell_size + "/")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

    if experiment == "Axl2_v4_exo_wide9":
        
        if cell_size == "standard":
            r_list = [2.5]
        elif cell_size == "medium":
            r_list = [3.0]
        elif cell_size == "large":
            r_list = [4.0]
        elif cell_size == "small":
            r_list = [2.0]
        crowding_list = [False]
        k16_list = [0.05]
        k25_list = [5.5]
        k24_list = [0.1]
        k20_list = [0.2, 0.3]
        k22_list = [10.5]
        k23_list = [26.0]
        k13_list = [1.5]
        k19_list = [4.5]
        DmFOO_list = [0.00045]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k20 in k20_list:
                    for DmFOO in DmFOO_list:
                        for i in range(len(k23_list)):
                            for k24 in k24_list:
                                for k25 in k25_list:
                                    if r == 2.0:
                                        read_run = "Small_pole/Actin" + "R2P0"
                                    elif r == 2.5:
                                        read_run = "Small_pole/Actin" + "R2P5"
                                    elif r == 3.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                    elif r == 4.0:
                                        read_run = "Small_pole/Actin" + "R3P0"
                                        
                                    mod_param_list.append(Common.ModelParametersAxl2V4(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                                       k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                                       k21=0.65, k19=4.5, gamma_k=0.005, k16=0.05, gamma_d=0.25*4.0, k14=0.0,
                                                                                       use_gamma=False, k13=1.5, k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                                       Kt=1.0, k20=k20, k22=k22_list[i], k23=k23_list[i], read_run=read_run,
                                                                                       k24=k24, k25=k25, DmFOO=DmFOO))
                                    exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 10,
                                                                                         np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 20.0]),
                                                                                         r, L = r, n_cables=10, use_cables=False,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                             gamma_P=1.0, cutoff_p_move=100, eps_hit=0.2))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2_v4(mod_param, True)

        sol_opt_list = [Common.SolverOption(term_max_it=150000, print_pwd_it=2000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v4/No_cable_wide9 + cell_size + "/"")

        Axl2_model_v4.run_axl2_v4(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)

    # Old model main simulations 
    if experiment == "Axl2_v1_exo_size_no_crowd_default_main":
        r_list = [2.5, 3.0]
        crowding_list = [False]
        k16_list = [0.05]
        k22_list = [3.0]
        k23_list = [25.0]
        k13_list = [5.0, 6.0]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k22 in k22_list:
                    for k23 in k23_list:
                        for k13 in k13_list:
                            if r == 2.0:
                                read_run = "Small_pole/Actin" + "R2P01"
                            elif r == 2.5:
                                read_run = "Small_pole/Actin" + "R2P51"
                            elif r == 3.0:
                                read_run = "Small_pole/Actin" + "R3P01"
                            elif r == 4.0:
                                read_run = "Small_pole/Actin" + "R3P01"
                            mod_param_list.append(Common.ModelParametersAxl2V1(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                               k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                               k21=0.65, k19=2.5, gamma_k=0.005, k16=0.05,
                                                                               gamma_d=0.25*4.0, k14=0.0, use_gamma=True, k13=k13,
                                                                               k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                               Kt=1.0, k20=0.5, k22=k22, k23=k23, read_run=read_run))
                            exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 9, np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                                                                                         r, L = r, n_cables=10, use_cables=True, centroid_cables=True, plot_cables=True,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                         i_S_states=[5,7, 9], gamma_P=1.0))
                                                    
        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2(mod_param, False)

        sol_opt_list = [Common.SolverOption(term_max_it=300000, print_pwd_it=2000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v1/Different_size_exo_default_main/")

        Axl2_model_v1.run_axl2_v1(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)

    if experiment == "Axl2_v1_exo_size_no_crowd_default_main_wide":

        r_list = [2.5, 3.0]
        crowding_list = [False]
        k16_list = [0.05]
        k22_list = [3.0]
        k23_list = [25.0]
        k13_list = [5.0, 6.0]
        mod_param_list = []
        exocyte_data_list = []
        for r in r_list:
            for crowding_p in crowding_list:
                for k22 in k22_list:
                    for k23 in k23_list:
                        for k13 in k13_list:
                            if r == 2.0:
                                read_run = "Small_pole/Actin" + "R2P01"
                            elif r == 2.5:
                                read_run = "Small_pole/Actin" + "R2P51"
                            elif r == 3.0:
                                read_run = "Small_pole/Actin" + "R3P01"
                            elif r == 4.0:
                                read_run = "Small_pole/Actin" + "R3P01"
                            mod_param_list.append(Common.ModelParametersAxl2V1(L=r, r=r, Dm=0.0045, Dm_gdp=0.0045, k2b=0.63, k5a=144,
                                                                               k5b=20.8, recruited=True, k17_alt=True, Dms=0.0045, k17=2.5*0.05,
                                                                               k21=0.65, k19=2.5, gamma_k=0.005, k16=0.05,
                                                                               gamma_d=0.25*4.0, k14=0.0, use_gamma=True, k13=k13,
                                                                               k12a=10.0, k12b=10.0, crowding_p=crowding_p,
                                                                               Kt=1.0, k20=0.5, k22=k22, k23=k23, read_run=read_run))
                            exocyte_data_list.append(Exo_endo_cytosis.ExocytosisData(0.05 * (r / 2.5), 9, np.array([0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                                                                                         r, L = r, n_cables=10, use_cables=True, centroid_cables=True, plot_cables=True,
                                                                                         lambda_hit=0.4, alpha=0.5, endocytosis=False, cables_window=True, i_P=6,
                                                                                         i_S_states=[5,7, 9], gamma_P=1.0, eps_hit=0.90))
                                                    
        for mod_param in mod_param_list:
            Solve_steday_state.solve_steady_state_axl2(mod_param, False)

        sol_opt_list = [Common.SolverOption(term_max_it=300000, print_pwd_it=2000, adaptive_rel_tol=5e-2, read_old_run=True, save_dist=10000)]
        file_loc = Common.FileLocations("Unit_mesh", "Axl2_v1/Different_size_exo_default_main_wide/")

        Axl2_model_v1.run_axl2_v1(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads, exocyte_data=exocyte_data_list, couple_mod_param_and_exo=True)

        sys.exit(0)      
