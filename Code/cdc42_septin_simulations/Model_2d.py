#!/usr/bin/env python

from dolfin import *
from tqdm import tqdm
import numpy as np
import pandas as pd
import scipy as sp 
import os
import sys 
import copy 
import time 
import multiprocessing as mp
import Common
import Config 

def init(l):
    Config.lock = l


class IC_2d(UserExpression):
    def __init__(self, *args, **kwargs):
        self.ic_Cdc42T = kwargs.pop('ic_Cdc42T')
        self.ic_Cdc42D = kwargs.pop('ic_Cdc42D')
        self.ic_BemGEF42 = kwargs.pop('ic_BemGEF42')
        self.ic_BemGEFm = kwargs.pop('ic_BemGEFm')
        self.sigma = kwargs.pop('sigma')
        super(IC_2d, self).__init__(*args, **kwargs)
    
    def eval(self, values, x):
        
        values[0] = self.ic_Cdc42T + np.random.normal(scale=0.001)
        values[1] = self.ic_Cdc42D + np.random.normal(scale=0.001)
        values[2] = self.ic_BemGEF42 + np.random.normal(scale=0.001)
        values[3] = self.ic_BemGEFm + np.random.normal(scale=0.001)
    
    def value_shape(self):
        return(4,)        


def formulate_components_fem_2d(mod_param, trial_f, trial_f_prev, test_f, measures, mesh, Cdc42c_exp, BemGEFc_exp, bud_scar_list=None):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm = trial_f
    Cdc42T_prev, Cdc42D_prev, BemGEF42_prev, BemGEFm_prev = trial_f_prev
    phi1, phi2, phi3, phi4 = test_f
    dx = measures
    
    mod_param_list, dx_list = calc_param_list(mod_param, bud_scar_list)
    mass_form, stiff_form, reaction_form = 0, 0, 0

    for i in range(len(dx_list)):
        k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, Dm_gdp = mod_param.make_param_constants()

        # Check if Cdc42-GTP should be capped or not, in case it is capped limit Cdc42T production 
        if mod_param.cap_CDC42T != None:
            k2a = conditional(lt(Cdc42T_prev, mod_param.cap_CDC42T), k2a, 0.0)
            k3 = conditional(lt(Cdc42T_prev, mod_param.cap_CDC42T), k3, 0.0)
            k4b = conditional(lt(Cdc42T_prev, mod_param.cap_CDC42T), k4b, 0.0)

        Cdc42T_mass = (Cdc42T - Cdc42T_prev)*phi1 * dx(dx_list[i])
        Cdc42D_mass = (Cdc42D - Cdc42D_prev)*phi2 * dx(dx_list[i])
        BemGEF42_mass = (BemGEF42 - BemGEF42_prev)*phi3 * dx(dx_list[i])
        BemGEFm_mass = (BemGEFm - BemGEFm_prev)*phi4 * dx(dx_list[i])
    
        Cdc42T_stiff = Dm * dot(grad(Cdc42T), grad(phi1)) * dx(dx_list[i])
        Cdc42D_stiff = Dm_gdp * dot(grad(Cdc42D), grad(phi2)) * dx(dx_list[i])
        BemGEF42_stiff = Dm * dot(grad(BemGEF42), grad(phi3)) * dx(dx_list[i])
        BemGEFm_stiff = Dm * dot(grad(BemGEFm), grad(phi4)) * dx(dx_list[i])

        # Cdc42T
        f1 = (k2a*BemGEFm_prev + k3*BemGEF42_prev) * Cdc42D_prev - (k2b + k4a*BemGEFm_prev + k7*BemGEFc_exp)*Cdc42T_prev + k4b*BemGEF42
        
        # Cdc42D
        f2 = k2b*Cdc42T - (k2a*BemGEFm_prev + k3*BemGEF42_prev)*Cdc42D_prev - k5b*Cdc42D + k5a*Cdc42c_exp
    
        # BemGEF42
        f3 = (k4a*BemGEFm_prev + k7*BemGEFc_exp)*Cdc42T_prev - k4b*BemGEF42

        # BemGEF
        f4 = k1a*BemGEFc_exp - k1b*BemGEFm + k4b*BemGEF42 - k4a*BemGEFm_prev*Cdc42T_prev
    
        reaction_form += (-f1*phi1 - f2*phi2 - f3*phi3 - f4*phi4) * dx(dx_list[i])
        mass_form += Cdc42T_mass + Cdc42D_mass + BemGEF42_mass + BemGEFm_mass
        stiff_form += Cdc42T_stiff + Cdc42D_stiff + BemGEF42_stiff + BemGEFm_stiff 
    
    return mass_form, stiff_form, reaction_form



def compute_fem_system_2d(mass_form, stiff_form, reaction_form):
    
    # Forms for vectors on RHS
    mass_form_rhs = rhs(mass_form)
    stiffness_reaction_form_rhs = rhs(stiff_form + reaction_form)
    
    mass_form_lhs = lhs(mass_form)
    stiffness_reaction_form_lhs = lhs(stiff_form + reaction_form)
    
    return mass_form_lhs, stiffness_reaction_form_lhs, mass_form_rhs, stiffness_reaction_form_rhs


def calc_ode_vals(U, mod_param, exp_list, eta_div_cell_area, dx):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm = U.split()
    Cdc42c_exp, BemGEFc_exp = exp_list
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm = U.split()
    
    x_old = np.zeros(2)
    # Cdc42 
    x_old[0] = exp_list[0].Cdc42c
    # BemGEFc
    x_old[1] = exp_list[1].BemGEFc
    Cdc42c_exp, BemGEFc_exp = exp_list

    ode_form1, ode_form2 = 0, 0
    
    k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, Dm_gdp = mod_param.return_param()
    eta_div_cell_area = Constant(eta_div_cell_area)
    ode_form1 = (eta_div_cell_area * (k5b*Cdc42D - k5a*Cdc42c_exp))*dx(1)
    ode_form2 = (eta_div_cell_area * (k1b*BemGEFm - k1a*BemGEFc_exp - k7*BemGEFc_exp*Cdc42T))*dx(1)

    ode_list = np.zeros(2)
    ode_list[0] = assemble(ode_form1)
    ode_list[1] = assemble(ode_form2)

    return ode_list


def compute_save_pole_data(U_high, exp_list_high, mesh, mod_param, eta_div_cell_area, dx):
    
    ode1, ode2 = calc_ode_vals(U_high, mod_param, exp_list_high, eta_div_cell_area, dx)
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm = U_high.split()
    Cdc42T_surface = assemble(Cdc42T * dx(1))
    Cdc42D_surface = assemble(Cdc42D * dx(1))
    BemGEF42_surface = assemble(BemGEF42 * dx(1))
    BemGEFm_surface = assemble(BemGEFm * dx(1))

    Cdc42T, Cdc42D, BemGEF42, BemGEFm = U_high.split(True)
    Cdc42c_exp, BemGEFc_exp = exp_list_high

    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    V = FunctionSpace(mesh, MixedElement([P1]))
    n = V.dim()
    d = mesh.geometry().dim()
    dof_coord = V.tabulate_dof_coordinates().reshape(n, d)
    dof_x = dof_coord[:, 0]
    dof_y = dof_coord[:, 1]
    dof_z = dof_coord[:, 2]
    Cdc42T_val = Cdc42T.vector().get_local()
    Cdc42D_val = Cdc42D.vector().get_local()
    BemGEF42_val = BemGEF42.vector().get_local()
    BemGEFm_val = BemGEFm.vector().get_local()

    # For the cytosolic values 
    cdc42_max = Cdc42T_val.max()
    cdc42_min = Cdc42T_val.min()
    max_TOL = 0.45 * (cdc42_max - cdc42_min)
    index_vec = abs(Cdc42T_val - cdc42_max) < max_TOL

    data_pole = pd.DataFrame({"x": dof_x,
                              "y": dof_y,
                              "z": dof_z,
                              "Cdc42T": Cdc42T_val, 
                              "Cdc42D" : Cdc42D_val, 
                              "BemGEF42" : BemGEF42_val, 
                              "BemGEFm" : BemGEFm_val})
    # Filter to keep the pole 
    data_pole = data_pole[index_vec]

    data_pole["Cdc42c"] = [Cdc42c_exp.Cdc42c for j in range(len(data_pole["x"]))]
    data_pole["BemGEFc"] = [BemGEFc_exp.BemGEFc for j in range(len(data_pole["x"]))]

    data_pole["Cdc42T_surface"] = [Cdc42T_surface for j in range(len(data_pole["x"]))]
    data_pole["Cdc42D_surface"] = [Cdc42D_surface for j in range(len(data_pole["x"]))]
    data_pole["BemGEF42_surface"] = [BemGEF42_surface for j in range(len(data_pole["x"]))]
    data_pole["BemGEFm_surface"] = [BemGEFm_surface for j in range(len(data_pole["x"]))]

    data_pole["Ode1"] = [ode1 for j in range(len(data_pole["x"]))]
    data_pole["Ode2"] = [ode2 for j in range(len(data_pole["x"]))]
    
    return data_pole 


def calc_find_pole_data(mesh, Cdc42T):
    
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    V = FunctionSpace(mesh, MixedElement([P1]))
    n = V.dim()
    d = mesh.geometry().dim()
    dof_coord = V.tabulate_dof_coordinates().reshape(n, d)
    dof_x = dof_coord[:, 0]
    dof_y = dof_coord[:, 1]
    dof_z = dof_coord[:, 2]
    cdc42_val = Cdc42T.vector().get_local()
    find_pole_data = pd.DataFrame({"x": dof_x,
                                   "y": dof_y,
                                   "z": dof_z,
                                   "u": cdc42_val})
    
    return find_pole_data


def calc_n_poles_2d(mesh, U_high, eps_clust, min_samples_clust):
    
    Cdc42T_tmp, Cdc42D_tmp, BemGEF42_tmp, BemGEFm_tmp = U_high.split(True)
    find_pole_data = calc_find_pole_data(mesh, Cdc42T_tmp)
    n_poles = Common.calc_n_poles_her(find_pole_data, eps_clust, min_samples_clust)
    
    return n_poles


def calc_cdc42_max(U_high):
    Cdc42T_tmp, Cdc42D_tmp, BemGEF42_tmp, BemGEFm_tmp = U_high.split(True)
    return Cdc42T_tmp.vector().max()


def calc_param_list(mod_param, bud_scar_param):

    mod_param_list = [mod_param]
    dx_list = [1]

    if bud_scar_param != None:
        dx_list.append(bud_scar_param.i_ring)
        dx_list.append(bud_scar_param.i_bud_scar)

        mod_param_list.append(bud_scar_param.param_ring)
        mod_param_list.append(bud_scar_param.param_bud_scar)

    return mod_param_list, dx_list


def solve_ode_system(U, mod_param, exp_list, dt, eta_div_cell_area, dx, param_bud_scar=None):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm = U.split()
    
    x_old = np.zeros(2)
    # Cdc42 
    x_old[0] = exp_list[0].Cdc42c
    # BemGEFc
    x_old[1] = exp_list[1].BemGEFc
    Cdc42c_exp, BemGEFc_exp = exp_list

    ode_form1, ode_form2 = 0, 0
    mod_param_list, dx_list = calc_param_list(mod_param, param_bud_scar)
    for i in range(len(dx_list)):
        k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, Dm_gdp = mod_param_list[i].return_param()
        eta_div_cell_area = Constant(eta_div_cell_area)
        ode_form1 += (eta_div_cell_area * (k5b*Cdc42D - k5a*Cdc42c_exp))*dx(dx_list[i])
        ode_form2 += (eta_div_cell_area * (k1b*BemGEFm - k1a*BemGEFc_exp - k7*BemGEFc_exp*Cdc42T))*dx(dx_list[i])

    ode_list_old = np.zeros(2)
    # Cdc42
    ode_list_old[0] = assemble(ode_form1)
    # BemGEFc
    ode_list_old[1] = assemble(ode_form2)
    
    # k1 RK2 method 
    k1 = dt * ode_list_old

    # k2 RK2 method 
    exp_list[0].Cdc42c += k1[0] * 0.5
    exp_list[1].BemGEFc += k1[1] * 0.5
    ode_list_old[0] = assemble(ode_form1)
    ode_list_old[1] = assemble(ode_form2)
    k2 = dt * ode_list_old
    
    x_new = x_old + 0.5*(k1 + k2)

    # Cdc42c
    exp_list[0].Cdc42c = x_new[0]
    # BemGEFc
    exp_list[1].BemGEFc = x_new[1]


def print_mass(U_high, dx, Cdc42_exp, BemGEF_exp):
    Cdc42T, Cdc42D, BemGEF42, BemGEFm = U_high.split()
    print("Mass Cdc42T = {:.3e}".format(assemble(Cdc42T * dx(1))))
    print("Mass Cdc42D = {:.3e}".format(assemble(Cdc42D * dx(1))))
    print("Mass BemGEF42 = {:.3e}".format(assemble(BemGEF42 * dx(1))))
    print("Mass BemGEFm = {:.3e}".format(assemble(BemGEFm * dx(1))))
    
    print("Mass Cdc42c = {:.3e}".format(Cdc42_exp.Cdc42c))
    print("Mass BemGEFc = {:.3e}".format(BemGEF_exp.BemGEFc))


def check_for_neg_conc(U, exp_list):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm = U.split()
    
    print("Checking negative concentrations")
    if Cdc42T.vector().min() < 0.0:
        print("Warning : Negative concentration in Cdc42T")
    if Cdc42D.vector().min() < 0.0:
        print("Warning : Negative concentration in Cdc42D")
    if BemGEF42.vector().min() < 0.0:
        print("Warning : Negative concentration BemGEF42")
    if BemGEFm.vector().min() < 0.0:
        print("Warning : Negative concentration BemGEFm")

    if exp_list[0].Cdc42c < 0.0:
        print("Warning : Negative concentration in Cdc42c")
    if exp_list[1].BemGEFc < 0.0:
        print("Warning : Negative concentration in BemGEFm")


def perform_time_step_dU(A_list, b_list, exp_list, dt_list, U_high, U_prev, U_low, U_norm, 
    U_old, exp_val_old, mod_param, solver, sol_opt, eta_div_cell_area, dx, delta_U):

    A, A_mass, A_stiff = A_list
    mass_rhs, stiffness_reaction_rhs = b_list

    dt_low, dt_high = dt_list

    while True:
        dt_use = dt_low

        A.zero()
        A.axpy(1.0, A_mass, True)
        A.axpy(dt_low, A_stiff, True)

        b_low = assemble(mass_rhs + Constant(dt_low) * stiffness_reaction_rhs)
        solver.solve(A,  U_low.vector(), b_low)
        delta_U.vector()[:] = U_low.vector() - U_prev.vector()

        # Update high accuracy system
        A.zero()
        A.axpy(1.0, A_mass, True)
        A.axpy(dt_high, A_stiff, True)
                
        # High accuracy update t -> t + dt_high
        U_high.vector()[:] = U_prev.vector() + dt_high / dt_low * delta_U.vector()
        U_prev.assign(U_high)
        solve_ode_system(U_high, mod_param, exp_list, dt_high, eta_div_cell_area, dx)

        # High accuracy update t + dt_high -> t + dt_high*2
        b_high = assemble(mass_rhs + Constant(dt_high)*stiffness_reaction_rhs)
        solver.solve(A,  U_high.vector(), b_high)
        U_prev.assign(U_high)

        U_norm.vector()[:] = U_low.vector() - U_high.vector()
            
        dt_use = dt_low
        dt_low, dt_high, progress = Common.update_step_length(U_norm, dt_low, sol_opt.adaptive_rel_tol, sol_opt.adaptive_rho)

        if progress == False:
            exp_list[0].Cdc42c = exp_val_old[0]
            exp_list[1].BemGEFc = exp_val_old[0]
            U_prev.vector()[:] = U_old.vector()
        else:
            # Only update ODE-system if time-step is accepted 
            solve_ode_system(U_high, mod_param, exp_list, dt_use/2.0, eta_div_cell_area, dx)
            break 

    return dt_low, dt_high, dt_use


def perform_time_step(A_list, b_list, exp_list, dt_list, U_high, U_prev, U_low, 
    U_norm, U_old, exp_val_old, mod_param, solver, sol_opt, eta_div_cell_area, dx):

    A, A_mass, A_stiff = A_list
    mass_rhs, stiffness_reaction_rhs = b_list

    dt_low, dt_high = dt_list

    while True:
        dt_use = dt_low

        A.zero()
        A.axpy(1.0, A_mass, True)
        A.axpy(dt_low, A_stiff, True)

        b_low = assemble(mass_rhs + Constant(dt_low) * stiffness_reaction_rhs)
        solver.solve(A,  U_low.vector(), b_low)

        # Update high accuracy system
        A.zero()
        A.axpy(1.0, A_mass, True)
        A.axpy(dt_high, A_stiff, True)
            
        # High accuracy update t -> t + dt_high
        b_high = assemble(mass_rhs + Constant(dt_high)*stiffness_reaction_rhs)
        solver.solve(A,  U_high.vector(), b_high)
        U_prev.assign(U_high)
        solve_ode_system(U_high, mod_param, exp_list, dt_high, eta_div_cell_area, dx)

        # High accuracy update t + dt_high -> t + dt_high*2
        b_high = assemble(mass_rhs + Constant(dt_high)*stiffness_reaction_rhs)
        solver.solve(A,  U_high.vector(), b_high)
        U_prev.assign(U_high)

        U_norm.vector()[:] = U_low.vector() - U_high.vector()
        dt_low, dt_high, progress = Common.update_step_length(U_norm, dt_low, sol_opt.adaptive_rel_tol, sol_opt.adaptive_rho)
        if progress == False:
            exp_list[0].Cdc42c = exp_val_old[0]
            exp_list[1].BemGEFc = exp_val_old[1]
            U_prev.vector()[:] = U_old.vector()
        else:
            # Only update ODE-system if time-step is accepted 
            solve_ode_system(U_high, mod_param, exp_list, dt_use / 2.0, eta_div_cell_area, dx)
            break 

    return dt_low, dt_high, dt_use


def perform_time_step_new(A_list, b_list, exp_list, dt_list, U_high, U_prev, U_low, 
    U_norm, U_old, exp_val_old, mod_param, solver, sol_opt, eta_div_cell_area, dx, k1, k2, q1, q2):

    A = A_list
    stiffness_reaction_rhs = b_list

    dt_low, dt_high = dt_list

    while True:
        dt_use = dt_low

        # Compute k1 
        b1 = assemble(Constant(dt_low) * stiffness_reaction_rhs)
        solver.solve(A,  k1.vector(), b1)
        q1.vector()[:] = k1.vector() + U_prev.vector()
        U_prev.vector()[:] = q1.vector()[:]

        # Compute second step and get the low solution 
        b2 = assemble(Constant(dt_low) * stiffness_reaction_rhs)
        solver.solve(A,  k2.vector(), b2)
        q2.vector()[:] = 0.75 * U_old.vector() + 0.25 * (q1.vector() + k2.vector())
        U_low.vector()[:] = q2.vector()
        U_prev.vector()[:] = q2.vector()

        # Update high accuracy system
        b2 = assemble(Constant(dt_low) * stiffness_reaction_rhs)
        solver.solve(A,  k2.vector(), b2)
        U_high.vector()[:] = 1.0/3.0 * U_old.vector() + 2.0/3.0*(q2.vector() + k2.vector())
            
        U_norm.vector()[:] = U_low.vector() - U_high.vector()
        dt_low, dt_high, progress = Common.update_step_length(U_norm, dt_low, sol_opt.adaptive_rel_tol, sol_opt.adaptive_rho)
        if progress == False:
            exp_list[0].Cdc42c = exp_val_old[0]
            exp_list[1].BemGEFc = exp_val_old[1]
            U_prev.vector()[:] = U_old.vector()
        else:
            # Only update ODE-system if time-step is accepted 
            solve_ode_system(U_high, mod_param, exp_list, dt_low, eta_div_cell_area, dx)
            break 

    return dt_low, dt_high, dt_use


def solve_2d(file_loc, mod_param, sol_opt, bud_scar_param, exocyte_data):
    
    Config.lock.acquire()
    file_loc.process_msh_file_to_xdmf_2d()
    file_loc.calc_pole_index()
    mesh, subdomains = file_loc.read_mesh_2d()
    Config.lock.release()
    
    dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
    
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement([P1, P1, P1, P1])
    H = FunctionSpace(mesh, element)

    if exocyte_data != None:
        exocyte_data.calc_coords(H)
        compute_time_exocyte = True
        t_hit = 0.0
    else:
        compute_time_exocyte = False
        t_hit = np.inf
    # If endocytosis is modelled to occur 
    if exocyte_data != None and exocyte_data.endocytosis == True:
        compute_time_endocyte = True
        t_hit_endo = 0.0
    else:
        compute_time_endocyte = False
        t_hit_endo = np.inf
    
    # Initial values
    if sol_opt.seed != False:
        np.random.seed(sol_opt.seed)
    u0_exp = IC_2d(ic_Cdc42T=mod_param.ic_Cdc42T,
                   ic_Cdc42D=mod_param.ic_Cdc42D,
                   ic_BemGEF42=mod_param.ic_BemGEF42,
                   ic_BemGEFm=mod_param.ic_BemGEFm,
                   sigma=0.1,
                   element=H.ufl_element())

    # Components for low accuracy soluations 
    U_prev = Function(H) 
    U_prev.interpolate(u0_exp)
    exp_list = [Expression('Cdc42c', Cdc42c = mod_param.ic_Cdc42I, element=P1), Expression('BemGEFc', BemGEFc = mod_param.ic_BemGEFc, element=P1)]
    U_low = Function(H)
    U_low.assign(U_prev)
    
    dt_low, dt_high = 1e-5, 1e-3 / 2

    # Weak formulation high accuracy Euler solution 
    test_functions = TestFunctions(H)
    trial_functions = TrialFunctions(H)
    exp_val_old = np.array([exp_list[0].Cdc42c, exp_list[1].BemGEFc])
    # The explicit-implicit scheme -> linear system. Here extracting components of said linear system
    if sol_opt.solver != "new_solver":
        mass, stiff, reaction = formulate_components_fem_2d(mod_param, trial_functions, U_prev, test_functions,
                dx, mesh, exp_list[0], exp_list[1], bud_scar_param)
        mass_lhs, stiffness_reaction_lhs, mass_rhs, stiffness_reaction_rhs = compute_fem_system_2d(mass, stiff, reaction)
        b_list = [mass_rhs, stiffness_reaction_rhs]

        mass_lhs_high = assemble(mass_lhs)
        stiffness_reaction_lhs_high = assemble(stiffness_reaction_lhs)
        A = mass_lhs_high + dt_low * stiffness_reaction_lhs_high
        A_mass = 1.0 * mass_lhs_high + 0.0 * stiffness_reaction_lhs_high
        A_stiff = 0.0 * mass_lhs_high + 1.0 * stiffness_reaction_lhs_high
        A_list = [A, A_mass, A_stiff]
    elif sol_opt.solver == "new_solver":
        mass_form_lhs, stiffness_reaction_form_rhs = formulate_components_fem_2d_exp(mod_param, trial_functions, U_prev, test_functions, dx,
             mesh, exp_list[0], exp_list[1])
        A_list = assemble(mass_form_lhs)
        b_list = stiffness_reaction_form_rhs
        k1, k2 = Function(H), Function(H)
        q1, q2 = Function(H), Function(H)


    if sol_opt.read_old_run == True:
        if mod_param.read_run != None:
            dir_read = "../../Intermediate/Experiments/Saved_runs/" + mod_param.read_run + "/"
        else:
            dir_read = file_loc.dir_read_run
        print("Reading old run from {}".format(dir_read))
        U_vector, exp_list_vector = Common.read_sim_result(dir_read)
        U_prev.vector()[:] = U_vector
        exp_list[0].Cdc42c = exp_list_vector[0]
        exp_list[1].BemGEFc = exp_list_vector[1]
        
    # Assigning initial values 
    U_high = Function(H)
    delta_U = Function(H)
    U_high.assign(U_prev)

    # Writing states to pvd-files
    t_curr = 0.0
    _Cdc42T_high, _Cdc42D_high, _BemGEF42_high, _BemGEFm_high = U_high.split()
    Cdc42T_high, Cdc42D_high, BemGEF42_high, BemGEFm_high = split(U_high)
    vtkfile_Cdc42T = File(file_loc.dir_pvd + "Cdc42T" + ".pvd")
    vtkfile_Cdc42T << (_Cdc42T_high, t_curr)

    # Actual solver to use 
    solver = KrylovSolver("gmres", "ilu")
    parameters['form_compiler']['optimize'] = True
    parameters['form_compiler']['cpp_optimize'] = True
    
    # Parameters for efficient time-stepping 
    eta_div_cell_area = mod_param.eta / (4*np.pi*(mod_param.r / mod_param.L)**2) # For ODE:s
    t_it = 1
    # Parameters and variables for adaptive time solver
    U_norm = Function(H)
    U_old = Function(H)
    # Old values for adaptive solver 
    U_old.vector()[:] = U_prev.vector()
    U_high.set_allow_extrapolation(True)
    
    # Termination criteria parameters 
    cdc42_limit_check_term = mod_param.ic_Cdc42T * sol_opt.term_inc_cdc42T
    cdc42_max_old = 0.0
    cdc42_max_best = 0.0
    term_crit = "Max_it"

    if sol_opt.term_time != False:
        term_time = sol_opt.term_time
    else:
        term_time = np.inf

    # Add actin cables if used 
    if exocyte_data != None and exocyte_data.use_cables == True and exocyte_data.cables_windows == False:
        Cdc42T_tmp, tmp1, tmp2, tmp3 = U_high.split(True)
        exocyte_data.calc_cable_attach(Cdc42T_tmp, 0, first=True)
        for i in range(1, exocyte_data.n_cables):
            exocyte_data.calc_cable_attach(Cdc42T_tmp, i, n_cables_use=i)

    if exocyte_data != None and exocyte_data.use_cables == True and exocyte_data.cables_windows == True:
        Cdc42T_tmp, tmp1, tmp2, tmp3 = U_high.split(True)
        for i in range(exocyte_data.n_cables):
            exocyte_data.calc_cable_attach_window(Cdc42T_tmp, i)

    # If location of cables are plotted 
    if exocyte_data != None and exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
        exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

    dt_use = 0.0
    for i in tqdm(range(sol_opt.term_max_it)):

        if t_curr > term_time:
            break

        if exocyte_data != None and compute_time_exocyte == True:
            t_hit = exocyte_data.calc_time_hit()
            compute_time_exocyte = False
            t_hit = t_curr + t_hit
            print("Time-hit = {}".format(t_hit))

        if exocyte_data != None and compute_time_endocyte == True:
            compute_time_endocyte = False
            t_hit_endo = exocyte_data.calc_time_endo()
            t_hit_endo += t_curr
            print("Time hit endo = {}".format(t_hit_endo))
        
        if sol_opt.solver == "dU":
            dt_low, dt_high, dt_use = perform_time_step_dU(A_list, b_list, exp_list, [dt_low, dt_high], U_high, U_prev, U_low, 
                U_norm, U_old, exp_val_old, mod_param, solver, sol_opt, eta_div_cell_area, dx, delta_U)
        elif sol_opt.solver == "two_step":
            dt_low, dt_high, dt_use = perform_time_step(A_list, b_list, exp_list, [dt_low, dt_high], 
                U_high, U_prev, U_low, U_norm, U_old, exp_val_old, mod_param, solver, sol_opt, 
                eta_div_cell_area, dx)
        elif sol_opt.solver == "new_solver":
            dt_low, dt_high, dt_use = perform_time_step_new(A_list, b_list, exp_list, [dt_low, dt_high], U_high, 
                    U_prev, U_low, U_norm, U_old, exp_val_old, mod_param, solver, sol_opt, eta_div_cell_area, 
                    dx, k1, k2, q1, q2)

        t_curr += dt_use

        # Perform endocytosis if required 
        if t_curr > t_hit_endo:
            print("Endocytosis occuring, t_curr = {:.3f}".format(t_curr))
            Cdc42T_tmp, tmp1, tmp2, tmp3 = U_high.split(True)
            i_hit = exocyte_data.calc_site_hit_endo(Cdc42T_tmp)
            coord_hit_sphere = exocyte_data.get_coord_sphere_hit(i_hit)
            for j in range(exocyte_data.n_coord):
                exocyte_data.calc_coord_interpolate(j, coord_hit_sphere, exocyte=False)
            exocyte_data.set_values_post_endocyte(U_high)
            U_high.vector()[:] = exocyte_data.U_new.vector()
            U_prev.assign(U_high)
            compute_time_endocyte = True
    
        # Check if exocytosis time should be computed 
        if t_curr > t_hit:
            print("Exocytosis occuring")
            Cdc42T_tmp, tmp1, tmp2, tmp3 = U_high.split(True)
            if exocyte_data.use_cables == False:
                i_hit = exocyte_data.calc_site_hit(Cdc42T_tmp)
            else:
                i_hit = exocyte_data.calc_site_hit_cable()

            coord_hit_sphere = exocyte_data.get_coord_sphere_hit(i_hit)
            for j in range(exocyte_data.n_coord):
                exocyte_data.calc_coord_interpolate(j, coord_hit_sphere)
            
            exocyte_data.set_values_post_exocyte(U_high)
            U_high.vector()[:] = exocyte_data.U_new.vector()
            U_prev.assign(U_high)
            compute_time_exocyte = True

            vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
            if exocyte_data != None and exocyte_data.plot_cables == True:
                exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

            # Correct the cytosolic diffusion using the practical total concentration 
            # obtained from solving for the steady state. 
            correct_cdc42 = mod_param.Cdc42_tot_pract - assemble((Cdc42T_high + Cdc42D_high + BemGEF42_high)*dx(1)) *eta_div_cell_area
            correct_BemGEF = mod_param.BemGEF_tot_pract - assemble((BemGEFm_high + BemGEF42_high)*dx(1)) *eta_div_cell_area

            exp_list[0].Cdc42c = correct_cdc42
            exp_list[1].BemGEFc = correct_BemGEF

        # Calculate new locations for cables if required 
        if exocyte_data != None and exocyte_data.use_cables == True:
            Cdc42T_tmp, tmp1, tmp2, tmp3 = U_high.split(True)
            exocyte_data.update_cable_pos(t_curr, Cdc42T_tmp)


        # Set low accuracy solution at high accuracy solution
                
        exp_val_old[0] = exp_list[0].Cdc42c
        exp_val_old[1] = exp_list[1].BemGEFc
        U_prev.vector()[:] = U_high.vector()
        U_old.vector()[:] = U_high.vector()
        
        # Export solution to disk 
        if t_it % sol_opt.print_pwd_it == 0:
            vtkfile_Cdc42T << (_Cdc42T_high, t_curr)            
            # If cable data is plotted 
            if exocyte_data != None and exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
                exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

        if t_it % sol_opt.print_pwd_it == 0 and t_curr > 0.0 and sol_opt.save_dist:
            data_pole = compute_save_pole_data(U_high, exp_list, mesh, mod_param, eta_div_cell_area, dx)
            file_name = file_loc.dir_intermediate + "Pole_dist_t.csv"
            n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
            Common.save_pole_dist_data(data_pole, t_curr, file_loc, file_name, n_poles)
        
        # Termination algorithm
        cdc42_max = calc_cdc42_max(U_high)
        if cdc42_max > cdc42_limit_check_term and sol_opt.term_check and t_curr > 50:
            
            # See if concentration flattens out or if concentration has peaked
            if cdc42_max > cdc42_max_best:
                cdc42_max_best = cdc42_max
                n_times_decrease = 0
            else:
                n_times_decrease += 1
            if abs(cdc42_max - cdc42_max_old) > (5e-5 * cdc42_max_old):
                n_times_flat = 0
            else:
                n_times_flat += 1
            
            if n_times_flat > sol_opt.term_times_decrease or n_times_decrease > sol_opt.term_times_flat:
                n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
                print("n_poles = {}".format(n_poles))
                if n_poles == 1:
                    print("Terminates upon a single pole")
                    print("n_poles = {}".format(n_poles))
                    term_crit = "One_pole"
                    break
            
            cdc42_max_old = cdc42_max
        
        t_it += 1
    
    # Write end coordinates
    vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
    n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
    
    # Calculate pole data and export to disk 
    Config.lock.acquire()
    Cdc42T_high, Cdc42D_high, BemGEF42_high, BemGEFm_high = U_high.split(True)
    pole_data = Common.PoleData(Cdc42T_high, mesh, t_curr, n_poles, term_crit=term_crit)
    pole_data.save_data(file_loc.dir_intermediate + "Pole_data.csv", file_loc)
    n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
    pole_data.print_f()
    print("n_poles = {:.0f}".format(n_poles))
    Config.lock.release()

    # If saving run-result 
    if sol_opt.save_run:
        if mod_param.save_run != None:
            dir_save = "../Intermediate/Experiments/Saved_runs/" + mod_param.save_run + "/"
        else:
            dir_save = file_loc.dir_save_run
        print("Saving results into {}".format(dir_save))
        exp_list_vector = np.array([exp_list[0].Cdc42c, exp_list[1].BemGEFc])
        Common.save_sim_result(U_high, exp_list_vector, dir_save)

    if sol_opt.save_dist:
        Config.lock.acquire()
        try:
            data_pole = compute_save_pole_data(U_high, exp_list, mesh, mod_param, eta_div_cell_area, dx)
            file_name = file_loc.dir_intermediate + "Pole_dist_end.csv"
            n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
            Common.save_pole_dist_data(data_pole, t_curr, file_loc, file_name, n_poles)
        except:
            print("Simulation failed saving pole data")
        Config.lock.release()
    
    if sol_opt.save_train == True:
        Cdc42T_tmp, Cdc42D_tmp, BemGEF42_tmp, BemGEFm_tmp = U_high.split(True)
        find_pole_data = calc_find_pole_data(mesh, Cdc42T_tmp)
        find_pole_data["index"] = [file_loc.run_index for j in range(len(find_pole_data["x"]))]
        file_path = file_loc.dir_intermediate + "Train_data.csv"
        if not os.path.isfile(file_path):
            find_pole_data.to_csv(file_path)
        else:
            find_pole_data.to_csv(file_path, header=False, mode='a')
        


def run_2d_model(file_loc, mod_param_list, sol_opt_list, times_run, n_threads, bud_scar_param=None, exocyte_data=None):

    try:
        len_exo_list = len(exocyte_data)
    except:
        len_exo_list = 0

    print("len_exo_list = ", len_exo_list)

    input_list = []

    if len_exo_list == 0:
        for mod_param in mod_param_list:
            for sol_opt in sol_opt_list:
                file_loc_tmp = Common.FileLocations(file_loc.mesh_name, file_loc.save_name, mod_param.calc_save_tag(bud_scar_param), name_store=file_loc.name_store, name_save=file_loc.name_save)
                input_list.append([copy.deepcopy(file_loc_tmp), copy.deepcopy(mod_param), copy.deepcopy(sol_opt), 
                    copy.deepcopy(bud_scar_param), copy.deepcopy(exocyte_data)])
    else:
        for mod_param in mod_param_list:
            for sol_opt in sol_opt_list:
                for exo_data in exocyte_data:
                    file_loc_tmp = Common.FileLocations(file_loc.mesh_name, file_loc.save_name, mod_param.calc_save_tag(bud_scar_param, exo_data), name_store=file_loc.name_store, name_save=file_loc.name_save)
                    input_list.append([copy.deepcopy(file_loc_tmp), copy.deepcopy(mod_param), copy.deepcopy(sol_opt), 
                        copy.deepcopy(bud_scar_param), copy.deepcopy(exo_data)])
    
    len_input_list = len(input_list)
    for i in range(times_run-1):
        for j in range(len_input_list):
            input_list.append(copy.deepcopy(input_list[j]))
    
    print("Starting simulations for 2d system")
    l = mp.Lock()
    pool = mp.Pool(initializer=init, initargs=(l,), processes=n_threads)
    result = pool.starmap(solve_2d, input_list)
    pool.close()
    pool.join()
    print("Done with simulations for 2d system")
    
