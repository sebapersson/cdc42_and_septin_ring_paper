#!/usr/bin/env python

# For running Cdc42 simulations. Because many different conditions 
# were explored there are many options, but the name should be clear. For 
# example experiment == "2d_model_nf_test_GAP" is the negative feedback model 
# where the effect of GAP activity is explored.
# Input for running is 
# ARGS[1] - experiment to run 
# ARGS[2] - number of processes (running in parallel)
# ARGS[3] - number of reactions (simulations are typically stochastic)
# Simulations should be run on a cluster, and they take typically at least 
# 100h to complete.

import Common 
import Model_2d
import Model_2d_nf
import Solve_steday_state

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
    
    if experiment == "2d_model_nf_large":
        mod_param_list = [Common.ModelParametersNf(L=7.0, r=7.0)]
        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state_nf(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=100000, print_pwd_it=200)]
        
        file_loc = Common.FileLocations("Unit_mesh", "2d_model_nf_large")
        
        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d_nf.run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)

    if experiment == "2d_model_nf_test_GAP":
        r_list = np.concatenate([np.linspace(2.5, 4.0, 20)])
        k8max_list = [0.0063*0.5]
        k9max_list = [0.0044*0.5]
        k2b_list = [0.35, 0.16]
        mod_param_list = []
        for r in r_list:
            for i in range(len(k9max_list)):
                for k2b in k2b_list:
                    mod_param_list.append(Common.ModelParametersNf(r=r, L=r, k8max=k8max_list[i], k9max=k9max_list[i], k2b=k2b))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state_nf(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=150000, print_pwd_it=2000, save_data_t=True, 
            save_data_t_it=10000, save_dist=True, term_check=False)]
        
        file_loc = Common.FileLocations("Unit_mesh", "2d_model_nf_test_GAP/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d_nf.run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)

    if experiment == "2d_model_pos_test_GAP":
        r_list = np.concatenate([np.linspace(2.5, 4.0, 20)])
        k2b_list = [0.63, 0.35]
        r_ref = 2.5
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        beta = 1.0
        mod_param_list = []
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref*beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            for k2b in k2b_list:
                mod_param_list.append(Common.ModelParametersDim(r=r, L=r, k2b=k2b, Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=600000, print_pwd_it=5000, save_data_t=True,
                                            save_data_t_it=5000, save_dist=True, term_check=False, term_time=6000.0,
                                            adaptive_rel_tol=5e-3)]
        
        file_loc = Common.FileLocations("Unit_mesh", "2d_model_pos_test_GAP/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d.run_2d_model(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)        

    if experiment == "2d_model_nf_second_feedback":
        r_list = np.linspace(3.0, 5.0, 12)
        k8max_list = [0.0063*17, 0.0063*20, 0.0063*23, 0.0063*30]
        k9max_list = [0.0044*17, 0.0044*20, 0.0044*23, 0.0044*30]
        mod_param_list = []
        for r in r_list:
            for i in range(len(k8max_list)):
                mod_param_list.append(Common.ModelParametersNf(r=r, L=r, k8max=0.0, k8max2=k8max_list[i], k9max=k9max_list[i], k8h2=0.2, second_feedback=True))                

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state_nf(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=150000, print_pwd_it=2000, save_data_t=True, 
            save_data_t_it=10000, save_dist=True, term_check=False)]

        file_loc = Common.FileLocations("Unit_mesh", "2d_model_nf_second_feedback/")
        
        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d_nf.run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)

    if experiment == "2d_pos_amount_weak_beta100":
        r_list = np.linspace(2.5, 4.0, 20)
        r_ref = 2.5
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        mod_param_list = []
        beta = 1.0
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref*beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            mod_param_list.append(Common.ModelParametersDim(r=r, L=r, k2b=0.35, Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=600000, print_pwd_it=10000, save_data_t=True,
                                            save_data_t_it=5000, save_dist=True, term_check=False, term_time=6000.0,
                                            adaptive_rel_tol=5e-3)]

        file_loc = Common.FileLocations("Unit_mesh", "2d_pos_weak_inc_beta100/")
        Model_2d.run_2d_model(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)
        
    if experiment == "2d_pos_amount_weak_beta75":
        r_list = np.linspace(2.5, 4.0, 20)
        r_ref = 2.5
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        mod_param_list = []
        beta = 0.75
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref*beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            mod_param_list.append(Common.ModelParametersDim(r=r, L=r, k2b=0.35, Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=1000000, print_pwd_it=10000, save_data_t=True,
                                            save_data_t_it=10000, save_dist=True, term_check=False, term_time=6000.0)]

        file_loc = Common.FileLocations("Unit_mesh", "2d_pos_weak_inc_beta75/")
        Model_2d.run_2d_model(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)

    if experiment == "2d_pos_amount_weak_beta50":
        r_list = np.linspace(2.5, 4.0, 20)
        r_ref = 2.5
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        mod_param_list = []
        beta = 0.50
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref *beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            mod_param_list.append(Common.ModelParametersDim(r=r, L=r, k2b=0.40, Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=1000000, print_pwd_it=10000, save_data_t=True,
                                            save_data_t_it=10000, save_dist=True, term_check=False, term_time=6000.0)]

        file_loc = Common.FileLocations("Unit_mesh", "2d_pos_weak_inc_beta50/")
        Model_2d.run_2d_model(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)

    if experiment == "2d_model_nf_beta100_strong_feedback":
        r_list = np.concatenate([np.linspace(2.5, 4.0, 20)])
        k8max_list = [0.0063*2.0]
        k9max_list = [0.0044*2.0]
        k2b_list = [0.35]
        r_ref = 2.5
        beta=1.0
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        mod_param_list = []
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref *beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            for i in range(len(k9max_list)):
                for k2b in k2b_list:
                    mod_param_list.append(Common.ModelParametersNf(r=r, L=r, k8max=k8max_list[i], k9max=k9max_list[i], k2b=k2b,
                                                                   Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state_nf(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=200000, print_pwd_it=2000, save_data_t=True, 
            save_data_t_it=10000, save_dist=True, term_check=False)]
        
        file_loc = Common.FileLocations("Unit_mesh", "2d_model_nf_beta100_strong_feedback/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d_nf.run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)        
        
    if experiment == "2d_model_nf_amount_weak_beta100":
        r_list = np.concatenate([np.linspace(2.5, 4.0, 20)])
        k8max_list = [0.0063*0.5]
        k9max_list = [0.0044*0.5]
        k2b_list = [0.35]
        r_ref = 2.5
        beta=1.0
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        mod_param_list = []
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref *beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            for i in range(len(k9max_list)):
                for k2b in k2b_list:
                    mod_param_list.append(Common.ModelParametersNf(r=r, L=r, k8max=k8max_list[i], k9max=k9max_list[i], k2b=k2b,
                                                                   Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state_nf(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=200000, print_pwd_it=2000, save_data_t=True, 
            save_data_t_it=10000, save_dist=True, term_check=False)]
        
        file_loc = Common.FileLocations("Unit_mesh", "2d_model_nf_amount_weak_beta100/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d_nf.run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)        
        
    if experiment == "2d_model_nf_amount_weak_beta75":
        r_list = np.concatenate([np.linspace(2.5, 4.0, 20)])
        k8max_list = [0.0063*0.5]
        k9max_list = [0.0044*0.5]
        k2b_list = [0.35]
        r_ref = 2.5
        beta=0.75
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        mod_param_list = []
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref *beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            for i in range(len(k9max_list)):
                for k2b in k2b_list:
                    mod_param_list.append(Common.ModelParametersNf(r=r, L=r, k8max=k8max_list[i], k9max=k9max_list[i], k2b=k2b,
                                                                   Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state_nf(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=200000, print_pwd_it=2000, save_data_t=True, 
            save_data_t_it=10000, save_dist=True, term_check=False)]
        
        file_loc = Common.FileLocations("Unit_mesh", "2d_model_nf_amount_weak_beta75/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d_nf.run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)        
        
    if experiment == "2d_model_nf_amount_weak_beta50":
        r_list = np.concatenate([np.linspace(2.5, 4.0, 20)])
        k8max_list = [0.0063*0.5]
        k9max_list = [0.0044*0.5]
        k2b_list = [0.35]
        r_ref = 2.5
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        beta=0.50
        mod_param_list = []
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref *beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            for i in range(len(k9max_list)):
                for k2b in k2b_list:
                    mod_param_list.append(Common.ModelParametersNf(r=r, L=r, k8max=k8max_list[i], k9max=k9max_list[i], k2b=k2b,
                                                                   Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state_nf(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=200000, print_pwd_it=2000, save_data_t=True, 
            save_data_t_it=10000, save_dist=True, term_check=False)]
        
        file_loc = Common.FileLocations("Unit_mesh", "2d_model_nf_amount_weak_beta50/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d_nf.run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)

    if experiment == "2d_model_nf_amount_weak_beta25":
        r_list = np.concatenate([np.linspace(2.5, 4.0, 10)])
        k8max_list = [0.0063*0.5]
        k9max_list = [0.0044*0.5]
        k2b_list = [0.35]
        r_ref = 2.5
        n_cdc42_ref = 1.0 * 4*np.pi/3*r_ref**3 * 1.0
        n_bemgef_ref = 0.017 * 4*np.pi/3*r_ref**3 * 1.0
        beta=0.25
        mod_param_list = []
        for r in r_list:
            ratio_new = ((r - r_ref)/r_ref *beta + 1)**3
            V_new = np.pi*4/3*r**3
            cdc42 = n_cdc42_ref*ratio_new / V_new
            bemgef = n_bemgef_ref*ratio_new / V_new
            for i in range(len(k9max_list)):
                for k2b in k2b_list:
                    mod_param_list.append(Common.ModelParametersNf(r=r, L=r, k8max=k8max_list[i], k9max=k9max_list[i], k2b=k2b,
                                                                   Cdc42_tot=cdc42, BemGEF_tot=bemgef))

        for mod_param in mod_param_list:
            Solve_steday_state.solve_for_steady_state_nf(mod_param)
        sol_opt_list = [Common.SolverOption(term_eps_clust=0.1, term_max_it=200000, print_pwd_it=2000, save_data_t=True, 
            save_data_t_it=10000, save_dist=True, term_check=False)]
        
        file_loc = Common.FileLocations("Unit_mesh", "2d_model_nf_amount_weak_beta25/")

        print("N_times_run = {}, n_threads = {}".format(n_times_run, n_threads))
        Model_2d_nf.run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, n_times_run, n_threads)
        sys.exit(0)                
        
    print("Provided experiment argv[1] = {} does not exist".format(sys.argv[1]))
    sys.exit(0)
