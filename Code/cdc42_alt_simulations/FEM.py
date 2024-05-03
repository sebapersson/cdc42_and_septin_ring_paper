"""
This file contains functions and classes related to the FEM-solver which solves the RD-system. 

Overall this file contains:
  Fem-option class (stores all the options of the FEM-solver)
  Functions for formulating the weak-formula, residual-form (adaptive time-stepping), formulating FEM-system, 
  calculate maximum concentration of u, calculating mass of u and extracting u-values on the membrane. 
"""


from dolfin import *
import numpy as np
import pandas as pd 
import Model_parameters


# Class holding the options for the FEM-solver,
# such as tolerances for adaptive step and
# pole finding parameters
# Args:
#    TOL_tadapt, tolerance for adaptive time-step
#    dt_max, maximum time-step length
#    initial_data_percent_pole,
#    tol_diff, the tolerance when terminating based on difference in u_max
#    termination_no_change, threshold for the number of steps that should not change before termination 
#    binary_percent_pole_min, 
#    binary_percent_pole_max,
#    boundary_percent_pole,
#    T, end time when solving RD-system
#    dt, time-step size when solving RD-system
#    steps_diff, how many steps taken before terminating on difference 
#    write_pvd, at which iterations to write pwd-file
#    solver, the solver to use for the linear FEM problem
#    seed, the seed to use in the simulations
#    write_pole_data, whether or not pole-data should be written to file as often as print_pvd
#    u0_factor, the factor that u0 should be scaled with for starting to look for termination
#    save_u_surface, whether or not the u-surface value should be saved upon termination
#    eps_clust, the epsilon parameter for the db-scan clustering
#    min_samples_clust, the minimum amount of neighbours for clustering
#    check_termination, whether or not the termination criteria are checked
#    print_v, option if v (inhibitor should be printed) or not
#    bud_scar_empty, whether or not the bud-scar are modeled by zero-flux Neuman-instead
class FemOpt:
    def __init__(self, TOL_tadapt=5e-4, dt_max = 0.05, initial_percent_pole = 0.1, tol_diff=1e-6,
                 termination_no_change=150, binary_percent_pole_min = 0.45, binary_percent_pole_max = 0.45,
                 boundary_percent_pole = 0.9, print_it=100, T=100, dt=1e-12, tol_step_diff=150,
                 write_pvd=100, max_it=np.inf, solver='gmres', pre_cond='ilu', seed=False,
                 write_pole_data=50, u0_factor=5, save_u_surface=False, write_pole_data_end=False,
                 eps_clust=0.1, min_samples_clust=5, check_termination=True, print_v=False, print_V=False,
                 bud_scar_empty=False, scale_conc=False, max_time=50, tol_ratio=1e-6, ratio_not_change=150):

        self.ratio_not_change = ratio_not_change
        self.tol_ratio = tol_ratio
        self.TOL_tadapt = TOL_tadapt
        self.dt_max = dt_max
        self.tol_diff = tol_diff
        self.initial_percent_pole = initial_percent_pole
        self.binary_percent_pole_min = binary_percent_pole_min
        self.binary_percent_pole_max = binary_percent_pole_max
        self.boundary_percent_pole = boundary_percent_pole
        self.print_it = print_it
        self.T = T
        self.dt = dt
        self.write_pvd = write_pvd
        self.max_it = max_it
        self.solver = solver
        self.pre_cond = pre_cond
        self.seed = seed
        self.write_pole_data = write_pole_data
        self.u0_factor = u0_factor
        self.termination_no_change = termination_no_change
        self.tol_step_diff = tol_step_diff
        self.save_u_surface = save_u_surface
        self.eps_clust = eps_clust
        self.min_samples_clust = min_samples_clust
        self.check_termination = check_termination
        self.print_v = print_v
        self.print_V = print_V
        self.bud_scar_empty = bud_scar_empty
        self.scale_conc = False
        self.max_time = max_time
        self.write_pole_data_end = write_pole_data_end


# Function that defines the components in the variational
# finite-element form.
# Args:
#     mod_param, the model parameters
#     measures, the different measures dx, ds, dS
#     trial_func, the trial functions u, v, V
#     trial_func_prev, the previous time-step of the trial functions 
#     mesh, the finite element mesh 
# Returns:
#     mass_form, mass-form for assembling problem 
#     stiffnes_form, stiffnes_form for assembling problem 
#     reaction_form, reactions form for assembling problem 
def formulate_components_fem(mod_param, trial_func, trial_func_prev,
                             test_func, measures, mesh, fem_opt, param_bud_scar=False):
    
    # Process the input arguments 
    u, v, V = trial_func
    u_prev, v_prev, V_prev = trial_func_prev
    phi_1, phi_2, phi_3 = test_func
    dx, ds, dS = measures
    
    # For calculating gradients easier 
    n = FacetNormal(mesh)
    def grad_G(w):
        return grad(w) - dot(grad(w), n)*n
    
    # Loop over the variational formulation to include subdomains (bud-scars)
    # Note: i_ds[1] should always be the main membrane 
    i_ds_list = [1]
    if param_bud_scar != False:
        if fem_opt.bud_scar_empty == False:
            # Barriers and old-bud-scars indices are lists 
            for i in param_bud_scar.i_old:
                i_ds_list.append(i)
            for i in param_bud_scar.i_barrier:
                i_ds_list.append(i)
            
            # Newest bud-scar and activation ring are scalars 
            i_ds_list.append(param_bud_scar.i_new)
        
        # Add ring
        for i in param_bud_scar.i_ring:
            i_ds_list.append(i)
    
    membrane_active_mass, membrane_active_stiffness = 0.0, 0.0
    membrane_inactive_mass, membrane_inactive_stiffness = 0.0, 0.0
    membrane_reaction, interface_reaction = 0.0, 0.0
    for i_ds in i_ds_list:
        # Use the correct parameters
        if i_ds == 1:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(mod_param)
        elif i_ds in param_bud_scar.i_old:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(param_bud_scar.param_old)
        elif i_ds in param_bud_scar.i_barrier:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(param_bud_scar.param_barrier)
        elif i_ds == param_bud_scar.i_new:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(param_bud_scar.param_new)
        elif i_ds in param_bud_scar.i_ring:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(param_bud_scar.param_ring)
        else:
            print("Error: ds-index {} is not marked".format(i_ds))
            sys.exit(1)
        
        # Reaction terms 
        f = ( (c2 * v) - u + (u_prev * u_prev * v_prev) ) 
        q = (c1 * V_prev * (c_max - (u_prev+v_prev))) - (c_1 * v)
        
        # Mass and stiffness matrices for membrane
        membrane_active_mass += (u-u_prev)*phi_1*ds(i_ds)
        membrane_active_stiffness += dot(grad_G(u), grad_G(phi_1))*ds(i_ds)
        membrane_inactive_mass += (v-v_prev)*phi_2*ds(i_ds)
        membrane_inactive_stiffness += d*dot(grad_G(v), grad_G(phi_2))*ds(i_ds)
        membrane_reaction += gamma*f*(phi_2 - phi_1)*ds(i_ds)
        
        # Mass and stiffnes matrix for membrane and cytosol 
        interface_reaction += gamma*q*(phi_3 - phi_2)*ds(i_ds)
    
    # Mass and stiffnes matrix for membrane and cytosol 
    cytosol_mass = (V - V_prev)*phi_3*dx(1)
    cytosol_stiffness = D*dot(grad(V), grad(phi_3))*dx(1)
    
    # Variational formulation 
    mass_form = membrane_active_mass + membrane_inactive_mass + cytosol_mass
    stiffness_form = (membrane_active_stiffness
                      + membrane_inactive_stiffness
                      + cytosol_stiffness)
    reaction_form = membrane_reaction + interface_reaction
    
    return mass_form, stiffness_form, reaction_form


# Function that defines the residual components used for terminating
# the simulations. Note, here u, v and V must have been redefined.
# Args:
#     mod_param, the model parameters
#     measures, the different measures dx, ds, dS
#     conc_func, the concentration functions 
#     trial_func_prev, the previous time-step of the trial functions
#     test_func, the test-functions 
#     mesh, the finite element mesh 
def formulate_residual_form(mod_param, measures, conc_func,
                            trial_func_prev, test_func, fem_opt, param_bud_scar=False):
    
    # The different functions used 
    u, v, V = conc_func
    u_prev, v_prev, V_prev = trial_func_prev
    phi_1, phi_2, phi_3 = test_func
    dx, ds, dS = measures
    
    # Loop over the variational formulation to include subdomains (bud-scars)
    # Note: i_ds[1] should always be the main membrane
    i_ds_list = [1]
    if param_bud_scar != False:
        if fem_opt.bud_scar_empty == False:
            # Barriers and old-bud-scars indices are lists 
            for i in param_bud_scar.i_old:
                i_ds_list.append(i)
            for i in param_bud_scar.i_barrier:
                i_ds_list.append(i)
            
            # Newest bud-scar and activation ring are scalars 
            i_ds_list.append(param_bud_scar.i_new)
        
        # Add ring
        for i in param_bud_scar.i_ring:
            i_ds_list.append(i)
    
    membrane_active_res, membrane_inactive_res = 0.0, 0.0
    membrane_reaction_res, membrane_res = 0.0, 0.0
    interface_res = 0.0
    for i_ds in i_ds_list:
        # Use the correct parameters
        if i_ds == 1:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(mod_param)
        elif i_ds in param_bud_scar.i_old:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(param_bud_scar.param_old)
        elif i_ds in param_bud_scar.i_barrier:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(param_bud_scar.param_barrier)
        elif i_ds == param_bud_scar.i_new:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(param_bud_scar.param_new)
        elif i_ds in param_bud_scar.i_ring:
            c1, c_1, c2, gamma, d, D, c_max = Model_parameters.convert_to_fenics_constants(param_bud_scar.param_ring)
        else:
            print("Error: ds-index {} is not marked".format(i_ds))
            sys.exit(1)
        
        f_res = c2*v - u + u*u*v
        q_res = c1*V*(c_max - (u+v)) - c_1*v
        
        membrane_active_res += (u-u_prev)*phi_1*ds(i_ds) 
        membrane_inactive_res += (v-v_prev)*phi_2*ds(i_ds) 
        membrane_reaction_res += gamma*f_res*(phi_2 - phi_1)*ds(i_ds)
        membrane_res += (membrane_active_res
                        + membrane_inactive_res
                        + membrane_reaction_res)
        
        interface_res += gamma*q_res*(phi_3 - phi_2)*ds(i_ds)
    
    cytosol_res = (V - V_prev)*phi_3*dx(1) + D*dot(grad(V), grad(phi_3))*dx(1)
    residual_form = membrane_res + interface_res + cytosol_res
    
    return residual_form


# Function that compute the components of the FEM equation-system. It will
# assemble the left-hand side, and formulates the right hand side (which is
# dynamic over time).
# Args:
#     mass_form, stiffnes_form, reaction_form of the problem
# Returns:
#     mass_matrix
#     stiffnes_reaction_matrix
#     mass_form_rhs
#     stiffnes_form_rhs 
def compute_components_fem_system(mass_form, stiffness_form, reaction_form):
    # Assemble lhs-side 
    mass_matrix = assemble(lhs(mass_form), keep_diagonal = True)
    mass_matrix.ident_zeros()
    stiffness_reaction_matrix = assemble(lhs(stiffness_form + reaction_form),
                                         keep_diagonal = True)
    stiffness_reaction_matrix.ident_zeros()
    
    # Forms for vectors on RHS
    mass_form_rhs = rhs(mass_form)
    stiffness_reaction_form_rhs = rhs(stiffness_form + reaction_form)
    
    return mass_matrix, stiffness_reaction_matrix, mass_form_rhs, stiffness_reaction_form_rhs


# Function that will check if iteration should be
# terminated based on negative concentrations
# Args:
#     mass_u/v/V, the mass of u, v and V
# Returns
#     true if should terminate
def negative_conc(mass_u, mass_v, mass_V):
    # Check mass conservation
    u_tot = assemble(mass_u)
    v_tot = assemble(mass_v)
    V_tot = assemble(mass_V)
    
    if u_tot < 0.0 or v_tot < 0.0 or V_tot < 0.0:
        print("Negative concentration")
        print("u_tot = {:.3f}, v_tot = {:.3f}, V_tot = {:.3f}".format(u_tot, v_tot, V_tot))
        return True
    else:
        return False


# Function that calculates u_max on the membrane.
# Args:
#    U, the values of the current solution
#    bV, the boundary mesh for the cell surface
# Returns:
#    u_max on the membrane
def calc_u_max(U, bV):
    _u, _v, _V = U.split()
    u_mem = interpolate(_u, bV)
    u_vec = u_mem.vector()
    return u_vec.max()


# Function that calculates the mass of the u-component in the model
# by including all the measures on the outer surface.
# Args:
#    u, the u-component
#    ds, the surface measures
#    param_bud_scar, the parameters for the bud-scars 
def calc_mass_u(u, ds, param_bud_scar=False):
    mass_u = u * ds(1)
    if param_bud_scar != 0:
        for i_ds in param_bud_scar.i_old:
            mass_u += u * ds(i_ds)
        for i_ds in param_bud_scar.i_barrier:
            mass_u += u * ds(i_ds)
        mass_u += u * ds(param_bud_scar.i_new)
        for i_ds in param_bud_scar.i_ring:
            mass_u += u * ds(i_ds)
    
    return mass_u


# Function that extracts the surface coordinates and the values of
# u on the surface (or rather the boundary mesh) and puts it into
# a pandas data frame.
# Args:
#    bV, the boundary mesh object
#    U, the current values of U
#    mesh, the mesh object for the function 
# Returns:
#   data_u_surface, the u-value on the surface 
def get_u_surface(bV, U, mesh):
    _u, _v, _V = U.split()
    n = bV.dim()  
    d = mesh.geometry().dim()
    dof_coord = bV.tabulate_dof_coordinates().reshape(n, d)
    dof_x = dof_coord[:, 0]
    dof_y = dof_coord[:, 1]
    dof_z = dof_coord[:, 2]
    u_val = (interpolate(_u, bV)).vector().get_local()
    data_u_surface = pd.DataFrame({"x": dof_x,
                                   "y": dof_y,
                                   "z": dof_z,
                                   "u": u_val})
    
    return data_u_surface
