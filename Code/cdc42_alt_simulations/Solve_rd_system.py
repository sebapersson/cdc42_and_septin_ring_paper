"""
This file simulates, that is solves, the RD-model. The file will solve 
the reaction-diffusion model using the user provided parameters which 
are converted to correct format in Model_parameters. By the default options 
in the FEM-options (FEM-file) the solver terminates when a single pole has been 
formed (number of poles are calculated via the algorithm in Calc_poles). The simulation 
is launched from the steady-state (Find_steady_state) and the mesh is read from disk 
into a correct format using Read_process_mesh-file. Lastly, the result if written to disk 
using the file-structure created by the File_system module. 

This file also contain function that will from provided parameters launch multiple RD-simulations, 
where each simulation will run in parallel as its own process. 
"""


import numpy as np
import pandas as pd 
import os
import time 
from dolfin import *
import multiprocessing as mp
import copy
import shutil

# Project related files
import FEM 
import Model_parameters
import File_system 
import Read_process_mesh
import Calc_poles
import Config 


# Function which creates a lock to avoid multiple processes accessing
# the same file. 
def init(l):
    Config.lock = l


# Function that will simulate the bulk-surface model for a specific
# experimental condition. The function will output the pvd-file that
# is saved every tenth second. If a part of the mesh, the user can
# also provide information regarding the bud-scars. Upon finding a pole
# the final pole-finding time, max-concentration and ratio
# of surface covered. If a file with max-concentration etc. already
# exists that file will be appended. 
# Args:
#     file_loc, the file-location object for the experiment
#     mod_param, the model parameters for the experiment
#     fem_opt, the options for the fem-solver
#     bud_scar_param, a ParamBudScars object (false by default by assuming no scars)
def solve_rd_system(file_loc, mod_param, fem_opt, param_bud_scar=False):
    
    # Process the mesh for the experiment, note multiple proc can access same file
    Config.lock.acquire()
    Read_process_mesh.process_msh_file_to_xdmf(file_loc)
    mesh, subdomains, surfaces = Read_process_mesh.read_mesh(file_loc)
    Config.lock.release()
    
    # Different measures, note that:
    # dx -> volumes (in this case)
    # ds -> exterior boundaries
    # dS -> interior boundaries
    dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
    ds = Measure("ds", domain=mesh, subdomain_data=surfaces)
    dS = Measure("dS", domain=mesh, subdomain_data=surfaces)
    measures = [dx, ds, dS]
    
    # Define the finite-element space used  
    P1 = FiniteElement('P', tetrahedron, 1)
    element = MixedElement([P1, P1, P1])
    H = FunctionSpace(mesh, element)
    
    # Trial and test-functions for variational formulation 
    phi_1, phi_2, phi_3 = TestFunctions(H)
    u, v, V = TrialFunctions(H) # Current time-step
    U_prev = Function(H) # Previous time-step 
    
    # Membrane excluding bud-scars (for pole-finding algorithm)
    bmesh = BoundaryMesh(mesh, "exterior")
    bmesh = SubMesh(bmesh, surfaces, 1)
    bV = FunctionSpace(bmesh, "P", 1)
    
    # The case we want to scale the absolute number of molecules
    if fem_opt.scale_conc != False:
        mod_param.V0 *= fem_opt.scale_conc * (5 / 2.5)
    
    # Interpolate the initial values for the problem, by using time 
    # in nano seconds none of the processes uses the same seed.
    if fem_opt.seed == False:
        my_seed = int(time.time_ns() % (2**32 - 1))
    else:
        my_seed = fem_opt.seed
    np.random.seed(my_seed)
    u0_exp = Model_parameters.ic(ic_u=mod_param.u0, ic_v=mod_param.v0,
                                 ic_V=mod_param.V0, sigma=0.05,
                                 find_mesh=file_loc.find_mesh,
                                 scale_sigma=mod_param.scale_sigma,
                                 element=H.ufl_element())
    U_prev.interpolate(u0_exp)
    u_prev, v_prev, V_prev = split(U_prev)
    
    # Calculate components for the variational problem
    mass_form, stiffness_form, reaction_form = FEM.formulate_components_fem(
        mod_param, [u, v, V], [u_prev, v_prev, V_prev], [phi_1, phi_2, phi_3],
        [dx, ds, dS], mesh, fem_opt, param_bud_scar)
    
    # Calculate LHS (these will not change) for the variational problem. Note
    # RHS will change with the dynamics of the solution
    mass_matrix, stiffness_reaction_matrix, mass_form_rhs, stiffness_reaction_form_rhs = (
        FEM.compute_components_fem_system(mass_form, stiffness_form, reaction_form))
    
    # Redefine u, v and V as functions to make further calculations possible 
    U = Function(H)
    u, v, V = split(U)
    U.assign(U_prev)
    _u, _v, _V = U.split()
    
    # Residual form for calculations
    residual_form = FEM.formulate_residual_form(
        mod_param, measures, [u, v, V], [u_prev, v_prev, V_prev],
        [phi_1, phi_2, phi_3], fem_opt, param_bud_scar)
    
    # Mass, used to check stability of solution 
    mass_v = v*ds(1)
    mass_V = V*dx(1)
    mass_u = FEM.calc_mass_u(u, ds, param_bud_scar)
    u_tot_prev = assemble(mass_u)
    
    # Pole finding parameters 
    cell_volume = assemble(Constant(1.0)*dx(1))
    u_mem = interpolate(_u, bV)
    u_vec = u_mem.vector()
    t_it = 0
    # Membrane area accounts for bud-scars
    membrane_area = assemble(Constant(1.0)*ds(1))
    print("membrane_area = {:.3f}".format(membrane_area))
    if param_bud_scar != False:
        membrane_area += assemble(Constant(1.0)*ds(param_bud_scar.i_new))
        for i in param_bud_scar.i_ring:
            membrane_area += assemble(Constant(1.0)*ds(i))
        for i in param_bud_scar.i_old:
            membrane_area += assemble(Constant(1.0)*ds(i))
        for i in param_bud_scar.i_barrier:
            membrane_area += assemble(Constant(1.0)*ds(i))
    
    # Create and write initial step to pvd-file 
    t = 0.0
    dt = fem_opt.dt
    pvd_folder = File_system.create_pwd_folder(file_loc)
    vtkfile_u = File(pvd_folder + "u" + ".pvd")
    if fem_opt.print_v:
        vtkfile_v = File(pvd_folder + "v" + ".pvd")
        vtkfile_v << (_v, t)
    if fem_opt.print_V:
        vtkfile_V = File(pvd_folder + "V_nuc" + ".pvd")
        vtkfile_V << (_V, t)
    vtkfile_u << (_u, t)
    
    if FEM.negative_conc(mass_u, mass_v, mass_V):
        print("Negative simulation before starting")
        print("Terminating simulation")
        print("File_result = {}".format(file_loc.result_folder))
        return 1 
    
    start_time = time.time()
    # Termination criteria
    u0_threshold = mod_param.u0 * fem_opt.u0_factor
    check_u_max = False
    u_max_best, u_max_old = 0, 0
    times_smaller_diff, times_smaller_max = 0, 0
    pole_ratio_check, times_pole_not_change = 0, 0
    # Adaptive time stepping with a termination cirteria
    while t < fem_opt.T:
        # Update current time and iteration number
        t_it += 1
        t += dt
        k = Constant(dt)
        
        # If transient dynamics are applied FEM-system must be assembled again
        # Solve linear variation problem for time-step 
        A = mass_matrix + dt*stiffness_reaction_matrix
        b = assemble(mass_form_rhs + k*stiffness_reaction_form_rhs)
        solve(A,  U.vector(), b, fem_opt.solver, fem_opt.pre_cond)
        
        if FEM.negative_conc(mass_u, mass_v, mass_V):
            print("ERROR: Negative amounts detected ==> System unstable")
            print("TERMINATING SIMULATION")
            print("File_result = {}".format(file_loc.result_folder))
            break
        
        # Check pole characteristics by performing 3 tests. A pole should pass all 3.
        # 1st test: Is the total amount of u constant?
        u_tot = assemble(mass_u)
        u_tot_slope = (u_tot - u_tot_prev)/dt
        u_tot_prev = u_tot
        
        # Update user about progress at specific iterations 
        if t_it % fem_opt.print_it == 0 or t_it == 1:
            real_time = t * mod_param.time_conversion_factor
            print("Step = {}, dt = {:.2e}, t = {:.2e}, T = {}, real_time = {:.2e}".format(t_it, dt, t, fem_opt.T, real_time))
            print("u_tot = {:.3f}, u_tot_slope = {:.3f}".format(u_tot, u_tot_slope))
            print("Solver = {}, pre_cond = {}".format(fem_opt.solver, fem_opt.pre_cond))
            # If user provides iterations as a breaking point 
            if fem_opt.max_it < t_it:
                print("Will break on maximum number of iterations")
                _u, _v, _V = U.split()
                vtkfile_u << (_u, t)
                break 
        
        # Check if pole-formation has occurred by checking if maximum concentration has been
        # been reached for u. Note, criteria not checked until u passes a certain level.
        # Furthermore, the code is also made to terminate if we reach a flack region.
        # If more than one pole is present upon termination, the code continues until
        # there is only one pole present. Furthermore, if time is larger than 25 minutes,
        # and there are two strong poles the code terminates. Note, termination is only
        # checked if provided as an option by the user. 
        u_max = FEM.calc_u_max(U, bV)
        time_real = t * mod_param.time_conversion_factor
        if u_max > u0_threshold and fem_opt.check_termination:
            check_u_max = True
            __pole_data = Calc_poles.calc_pole_data(U, bV, fem_opt, t)
            pole_ratio_ = __pole_data.pole_surface_ratio
        if check_u_max:
            if u_max > u_max_best:
                pole_param_best_diff = Calc_poles.calc_pole_data(U, bV, fem_opt, t)
                u_max_best = u_max
                times_smaller_max = 0
            else:
                times_smaller_max += 1
            
            if abs(u_max - u_max_old) > (fem_opt.tol_diff * u_max_old):
                pole_param_best_max = Calc_poles.calc_pole_data(U, bV, fem_opt, t)
                times_smaller_diff = 0
            else:
                times_smaller_diff += 1

            if abs(pole_ratio_ - pole_ratio_check) > fem_opt.tol_ratio*pole_ratio_check:
                times_pole_not_change = 0
                pole_ratio_check = pole_ratio_
            else:
                times_pole_not_change += 1
                
            if times_smaller_diff > fem_opt.tol_step_diff and times_pole_not_change > fem_opt.ratio_not_change:
                # Check the number of poles
                u_surface_data = FEM.get_u_surface(bV, U, mesh)
                n_poles = Calc_poles.calc_n_poles_her(u_surface_data, fem_opt.eps_clust,
                                           fem_opt.min_samples_clust)
                if n_poles == 1:
                    criteria_term = "difference"
                    pole_param_best = pole_param_best_diff
                    u_surf_best = u_surface_data
                    _u, _v, _V = U.split()
                    vtkfile_u << (_u, t)
                    print("Pole found: Terminates on u_max being flat")
                    break
                else:
                    times_smaller_diff = 0
            if times_smaller_max > fem_opt.termination_no_change and times_pole_not_change > fem_opt.ratio_not_change:
                u_surface_data = FEM.get_u_surface(bV, U, mesh)
                n_poles = Calc_poles.calc_n_poles_her(u_surface_data, fem_opt.eps_clust, fem_opt.min_samples_clust)
                if n_poles == 1:
                    criteria_term = "maximum"
                    pole_param_best = pole_param_best_max
                    _u, _v, _V = U.split()
                    u_surf_best = u_surface_data
                    vtkfile_u << (_u, t)
                    print("Pole found: Terminates on u_max decreasing")
                    break
                else:
                    times_smaller_max = 0
            
        if time_real > fem_opt.max_time:
            print("Terminates on to long time")
            print("Time_real = {:.3e}".format(time_real))
            # Calculate number of poles with large u_max
            # to see if there are two strong poles
            u_surface_data = FEM.get_u_surface(bV, U, mesh)
            n_poles = Calc_poles.calc_n_poles(u_surface_data, u_max_filt=0.9)
            if n_poles > 1:
                criteria_term = "multiple_poles"
                pole_param_best = Calc_poles.calc_pole_data(U, bV, fem_opt, t)
                u_surf_best = u_surface_data
                _u, _v, _V = U.split()
                vtkfile_u << (_u, t)
                print("Terminates on multiple poles")
                break
            else:
                criteria_term = "one_pole_long_time"
                pole_param_best = Calc_poles.calc_pole_data(U, bV, fem_opt, t)
                u_surf_best = u_surface_data
                _u, _v, _V = U.split()
                vtkfile_u << (_u, t)
                print("Terminates on one pole but long time")
                break        
        
        # Compute next time-step adaptively using residuals 
        dt_old = dt
        R = assemble(residual_form)
        dt = fem_opt.TOL_tadapt / norm(R, 'l2')
        dt = min(2.0*dt_old*dt/(dt_old + dt), fem_opt.dt_max)
        
        # Only save data sparsely  
        if t_it % fem_opt.write_pvd == 0 or t + dt >= fem_opt.T:
            _u, _v, _V = U.split()
            vtkfile_u << (_u, t)
            if fem_opt.print_v:
                vtkfile_v << (_v, t)
            if fem_opt.print_V:
                vtkfile_V << (_V, t)
        # The case pool data is written
        if (fem_opt.write_pole_data != False and t_it % fem_opt.write_pole_data == 0):
            pole_data = Calc_poles.calc_pole_data(U, bV, fem_opt, t)
            should_terminate = Calc_poles.write_pole_to_file(pole_data, file_loc)
            if should_terminate:
                print("Simulations crashed")
                return 0 
        
        # Update previous solution
        U_prev.assign(U)
        u_max_old = u_max
    end_time = time.time()
    run_time = (end_time - start_time) / 60
    
    # Write u_max, pole surface ratio, end-time and run-time to file
    # If provided, also write u_surface distribution
    # Exception handling important to avoid code getting stuck which
    # happens upon crash inside a lock. 
    Config.lock.acquire()
    try:
        if fem_opt.write_pole_data_end == True:
            print("Terminating on time")
            pole_param_best = Calc_poles.calc_pole_data(U, bV, fem_opt, t)
            u_surface_data = FEM.get_u_surface(bV, U, mesh)
            t = pole_param_best.t
            u_max = pole_param_best.u_max
            ratio = pole_param_best.pole_surface_ratio
            ratio_hard = pole_param_best.pole_surface_ratio_hard
            binary = pole_param_best.pole_percent_binary
            criteria_term = "end_time"
            n_poles = Calc_poles.calc_n_poles_her(u_surface_data, fem_opt.eps_clust,
                                                  fem_opt.min_samples_clust)
        
        elif t > fem_opt.T or t_it > fem_opt.max_it or should_terminate:
            print("Pole has likely not been found")
            ratio = 0.0
            ratio_hard = 0.0
            u_max = 0.0
            binary = 0.0
            n_poles = 0.0
            if should_terminate:
                criteria_term = "system_crashed"
            else:
                criteria_term = "fail"
            u_surf_best = Calc_poles.calc_pole_data(U, bV, fem_opt, t)
        else:
            t = pole_param_best.t
            u_max = pole_param_best.u_max
            ratio = pole_param_best.pole_surface_ratio
            ratio_hard = pole_param_best.pole_surface_ratio_hard
            binary = pole_param_best.pole_percent_binary
        # Get mass-values 
        v_tot = assemble(mass_v)
        u_tot = assemble(mass_u)
        V_tot = assemble(mass_V)
        data_save = pd.DataFrame({"t_end": [t],
                                  "u_max": [u_max],
                                  "ratio_hard": [ratio_hard],
                                  "ratio": [ratio],
                                  "binary" : [binary],
                                  "run_time": [run_time],
                                  "term_criteria": [criteria_term],
                                  "v_tot": [v_tot],
                                  "u_tot": [u_tot],
                                  "V_tot": [V_tot],
                                  "n_poles": [n_poles]})
        
        if not os.path.isfile(file_loc.intermediate_file):
            data_save.to_csv(file_loc.intermediate_file)
        else:
            data_save.to_csv(file_loc.intermediate_file, header=False, mode='a')
        
        # Save the surface coordinates if provided
        if fem_opt.save_u_surface == True:
            u_surf_best["index"] = file_loc.value_write_pole
            if not os.path.isfile(file_loc.u_surface):
                u_surf_best.to_csv(file_loc.u_surface)
            else:
                u_surf_best.to_csv(file_loc.u_surface, header=False, mode='a')
    except:
        data_save = pd.DataFrame({"t_end": [0.0],
                                  "u_max": [0.0],
                                  "ratio_hard": [0.0],
                                  "ratio": [0.0],
                                  "binary" : [0.0],
                                  "run_time": [0.0],
                                  "term_criteria": ["crashed"],
                                  "v_tot": [0.0],
                                  "u_tot": [0.0],
                                  "V_tot": [0.0],
                                  "n_poles": [0.0]})
        
        if not os.path.isfile(file_loc.intermediate_file):
            data_save.to_csv(file_loc.intermediate_file)
        else:
            data_save.to_csv(file_loc.intermediate_file, header=False, mode='a')
    
    Config.lock.release()
    
    return 0


# Function that will simulate a provided experimental setup, note that for each
# parameter setting, or mesh, the experiment will be run a certain number of
# times, and the result will be written to file.
# Args:
#     param_list_dim, a list with parameters to use (each entry is a model
#         parameter object with dimensions)
#     bud_scar_list, a list with the number of bud-scars to investigate
#     fem_opt, the FEM-solver options
#     times_run, the number of times to run the simulations
#     nucleus, bool stating whether or not the mesh has a nucleus 
#     delete_old_experiment, whether or not intermediate old-files should be deleted
#     param_bud_scar, a param_bud_scar parameter object that is provided if bud-scars are used
#     tag_exp, a string tagging the experiment when saving the result
# Returns:
#    void 
def run_rd_simulations(param_list_dim, bud_scar_list, fem_opt_list=[FEM.FemOpt()], times_run=5, nucleus=True,
                       delete_old_experiment=False, delete_old_pvd=False, n_threads=4, check_mesh_size=False,
                       param_bud_scar_list=False, tag_exp=""):
        
    # Convert the parameters to non-dimensional parameters
    param_list = []
    for param_dim in param_list_dim:
        param_to_append = Model_parameters.make_param_non_dimension(param_dim)
        if param_to_append == False:
            continue
        else:
            param_list.append(param_to_append)
    
    # If the check mesh size is an option
    if check_mesh_size != False:
        dens_mesh_list = check_mesh_size
    else:
        dens_mesh_list = [False]
    
    # Delete experiment if provided by user, note avoid deleting all results 
    if tag_exp != "":
        dir_intermeidate = "../../Intermediate/Experiments/" + tag_exp
        if delete_old_experiment and os.path.isdir(dir_intermeidate):
            shutil.rmtree(dir_intermeidate)
    
    # Create the list which to run 
    comb_run_parts = []
    for fem_opt in fem_opt_list:
        for dens_mesh in dens_mesh_list:
            for n_bud_scars in bud_scar_list:
                if not Read_process_mesh.msh_file_exists(n_bud_scars, dens_mesh, nucleus):
                    return 1
                
                for param_dim in param_list_dim:
                    mod_param = Model_parameters.make_param_non_dimension(param_dim)
                    # Add bud-scars parameters if provided
                    fem_opt_deep = copy.deepcopy(fem_opt)
                    if param_bud_scar_list == False:
                        param_bud_scar = False
                        file_loc = File_system.FileLocations(nucleus, param_dim, n_bud_scars,
                                                 dens_mesh, param_bud_scar, tag_exp, fem_opt_deep)
                        comb_run_parts.append([file_loc, mod_param, fem_opt_deep])
                    else:
                        for param_bud_scar in param_bud_scar_list:
                            param_bud_scar_deep = copy.deepcopy(param_bud_scar)
                            file_loc = File_system.FileLocations(nucleus, param_dim, n_bud_scars,
                                                     dens_mesh, param_bud_scar, tag_exp, fem_opt_deep)
                            comb_run_parts.append([file_loc, mod_param, fem_opt_deep, param_bud_scar_deep])
                        
                    if delete_old_pvd and os.path.isdir(file_loc.pvd_folder):
                        shutil.rmtree(file_loc.pvd_folder)
                        os.mkdir(file_loc.pvd_folder)
    
    # Set the input combinations to the number to run 
    combinations_to_run = []
    for comb in comb_run_parts:
        for i in range(times_run):
            combinations_to_run.append(copy.deepcopy(comb))
    # Ensure correct pole index file if saving pole-data
    for comb in combinations_to_run:
        comb[0].calc_pole_index()
    
    # Run the code in parallel
    print("Starting simulations")
    l = mp.Lock()
    pool = mp.Pool(initializer=init, initargs=(l,), processes=n_threads)
    result = pool.starmap(solve_rd_system, combinations_to_run)
    pool.close()
    pool.join()



