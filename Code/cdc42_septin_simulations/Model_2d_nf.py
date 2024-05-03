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

# TODO : Add common function for calculating pole-data 


class IC_2d_nf(UserExpression):
    def __init__(self, *args, **kwargs):
        self.ic_Cdc42T = kwargs.pop('ic_Cdc42T')
        self.ic_Cdc42D = kwargs.pop('ic_Cdc42D')
        self.ic_BemGEF42 = kwargs.pop('ic_BemGEF42')
        self.ic_BemGEFm = kwargs.pop('ic_BemGEFm')
        self.ic_BemGEF42_star = kwargs.pop('ic_BemGEF42_star')
        self.ic_BemGEFm_star = kwargs.pop('ic_BemGEFm_star')
        self.sigma = kwargs.pop('sigma')
        super(IC_2d_nf, self).__init__(*args, **kwargs)
    
    def eval(self, values, x):
        
        values[0] = self.ic_Cdc42T + np.random.normal(scale=0.001)
        values[1] = self.ic_Cdc42D + np.random.normal(scale=0.001)
        values[2] = self.ic_BemGEF42 + np.random.normal(scale=0.001)
        values[3] = self.ic_BemGEFm + np.random.normal(scale=0.001)
        values[4] = self.ic_BemGEF42_star + np.random.normal(scale=0.001)
        values[5] = self.ic_BemGEFm_star + np.random.normal(scale=0.001)
    
    def value_shape(self):
        return(6,)        


def formulate_components_fem_2d_nf(mod_param, trial_f, trial_f_prev, test_f, measures, exp_list):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star = trial_f
    Cdc42T_prev, Cdc42D_prev, BemGEF42_prev, BemGEFm_prev, BemGEF42_star_prev, BemGEFm_star_prev = trial_f_prev
    phi1, phi2, phi3, phi4, phi5, phi6 = test_f
    dx = measures

    Cdc42c_exp, BemGEFc_exp, BemGEFc_star_exp = exp_list
    
    mass_form, stiff_form, reaction_form = 0, 0, 0
    
    k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, k8max, k8h, k8n, k9max, k9h, k9n, Dm_gdp = mod_param.make_param_constants()

    Cdc42T_mass = (Cdc42T - Cdc42T_prev)*phi1 * dx(1)
    Cdc42D_mass = (Cdc42D - Cdc42D_prev)*phi2 * dx(1)
    BemGEF42_mass = (BemGEF42 - BemGEF42_prev)*phi3 * dx(1)
    BemGEFm_mass = (BemGEFm - BemGEFm_prev)*phi4 * dx(1)
    BemGEF42_star_mass = (BemGEF42_star - BemGEF42_star_prev)*phi5 * dx(1)
    BemGEFm_star_mass = (BemGEFm_star - BemGEFm_star_prev)*phi6 * dx(1)
    
    Cdc42T_stiff = Dm * dot(grad(Cdc42T), grad(phi1)) * dx(1)
    Cdc42D_stiff = Dm_gdp * dot(grad(Cdc42D), grad(phi2)) * dx(1)
    BemGEF42_stiff = Dm * dot(grad(BemGEF42), grad(phi3)) * dx(1)
    BemGEFm_stiff = Dm * dot(grad(BemGEFm), grad(phi4)) * dx(1)
    BemGEF42_star_stiff = Dm * dot(grad(BemGEF42_star), grad(phi5)) * dx(1)
    BemGEFm_star_stiff = Dm * dot(grad(BemGEFm_star), grad(phi6)) * dx(1)

    # Regulate strength of nf by altering the maximum phosphorylation rate 
    strength_feedback = mod_param.strength_feedback

    # Helper terms 
    k8 = k8max * (BemGEF42_prev + BemGEF42_star_prev)**k8n / (k8h**k8n + (BemGEF42_prev + BemGEF42_star_prev)**k8n)
    k8 *= strength_feedback
    if mod_param.second_feedback == True:
        k8max2, k8h2 = Constant(mod_param.k8max2), Constant(mod_param.k8h2)
        k8_p2 = k8max2 * 1 / (1 + exp(-k8h2*(Cdc42T_prev - 150)))
        k8_p1 = k8max * (BemGEF42_prev + BemGEF42_star_prev)**k8n / (k8h**k8n + (BemGEF42_prev + BemGEF42_star_prev)**k8n)
        k8 = (k8_p1 + k8_p2)
    
    k7_term = k7 * (BemGEFc_star_exp  + BemGEFc_exp)
    k4a_term = k4a * (BemGEFm_prev + BemGEFm_star_prev)
    k4b_term = k4b * (BemGEF42 + BemGEF42_star)
    k8_term = k8 * (BemGEF42_prev + BemGEF42_star_prev) * BemGEF42_prev

    # Cdc42T
    f1 = (k2a*BemGEFm_prev + k3*BemGEF42_prev) * Cdc42D_prev - (k2b + k4a_term + k7_term)*Cdc42T_prev + k4b_term
        
    # Cdc42D
    f2 = k2b*Cdc42T - (k2a*BemGEFm_prev + k3*BemGEF42_prev)*Cdc42D_prev - k5b*Cdc42D + k5a*Cdc42c_exp
    
    # BemGEF42
    f3 = (k4a*BemGEFm_prev + k7*BemGEFc_exp)*Cdc42T_prev - k4b*BemGEF42 - k8_term

    # BemGEF
    f4 = k1a*BemGEFc_exp - k1b*BemGEFm + k4b*BemGEF42 - k4a*BemGEFm_prev * Cdc42T_prev

    # BemGEF42*
    f5 = (k4a * BemGEFm_star_prev + k7*BemGEFc_star_exp)*Cdc42T_prev - k4b*BemGEF42_star + k8_term

    # BemGEF*
    f6 = k1a*BemGEFc_star_exp - k1b*BemGEFm_star + k4b*BemGEF42_star - k4a*BemGEFm_star_prev * Cdc42T_prev
    
    reaction_form += (-f1*phi1 - f2*phi2 - f3*phi3 - f4*phi4 - f5*phi5 - f6*phi6) * dx(1)
    mass_form += Cdc42T_mass + Cdc42D_mass + BemGEF42_mass + BemGEFm_mass + BemGEF42_star_mass + BemGEFm_star_mass
    stiff_form += Cdc42T_stiff + Cdc42D_stiff + BemGEF42_stiff + BemGEFm_stiff + BemGEF42_star_stiff + BemGEFm_star_stiff
    
    return mass_form, stiff_form, reaction_form


def compute_fem_system_2d_nf(mass_form, stiff_form, reaction_form):
    
    # Forms for vectors on RHS
    mass_form_rhs = rhs(mass_form)
    stiffness_reaction_form_rhs = rhs(stiff_form + reaction_form)
    
    mass_form_lhs = lhs(mass_form)
    stiffness_reaction_form_lhs = lhs(stiff_form + reaction_form)
    
    return mass_form_lhs, stiffness_reaction_form_lhs, mass_form_rhs, stiffness_reaction_form_rhs


def solve_ode_system(U, mod_param, exp_list, dt, eta_div_cell_area, dx):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star = U.split()
    Cdc42c_exp, BemGEFc_exp, BemGEFc_star_exp = exp_list
    
    x_old = np.zeros(3)
    # Cdc42 
    x_old[0] = exp_list[0].Cdc42c
    # BemGEFc
    x_old[1] = exp_list[1].BemGEFc
    # BemGEFc_start
    x_old[2] = exp_list[2].BemGEFc_star

    Cdc42c_exp, BemGEFc_exp, BemGEFc_star_exp = exp_list

    ode_form1, ode_form2, ode_form3 = 0.0, 0.0, 0.0

    k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, k8max, k8h, k8n, k9max, k9h, k9n, Dm_gdp = mod_param.return_param()    
    eta_div_cell_area = Constant(eta_div_cell_area)

    # Regulate strength of the negative feedback 
    strength_feedback = mod_param.strength_feedback
    k9 = k9max * (exp_list[2].BemGEFc_star**k9n) / (k9h**k9n + exp_list[2].BemGEFc_star**k9n)
    k9 *= strength_feedback
    k9_term = Constant(k9 * exp_list[2].BemGEFc_star)

    ode_form1 += eta_div_cell_area * (k5b*Cdc42D - k5a*Cdc42c_exp) * dx(1)
    ode_form2 += eta_div_cell_area * (k1b * BemGEFm - k1a * BemGEFc_exp - k7*BemGEFc_exp*Cdc42T) * dx(1) 
    ode_form3 += eta_div_cell_area * (k1b * BemGEFm_star - k1a * BemGEFc_star_exp - k7*BemGEFc_star_exp*Cdc42T) *dx(1)

    ode_list_old = np.zeros(3)
    # Cdc42
    ode_list_old[0] = assemble(ode_form1)
    # BemGEFc
    ode_list_old[1] = assemble(ode_form2) + k9_term
    # BemGEFs_star 
    ode_list_old[2] = assemble(ode_form3) - k9_term
    
    # k1 RK2 method 
    k1 = dt * ode_list_old

    # k2 RK2 method 
    exp_list[0].Cdc42c += k1[0] * 0.5
    exp_list[1].BemGEFc += k1[1] * 0.5
    exp_list[2].BemGEFc_star += k1[2] * 0.5

    k9 = k9max * (exp_list[2].BemGEFc_star**k9n) / (k9h**k9n + exp_list[2].BemGEFc_star**k9n)
    k9_term = Constant(k9 * exp_list[2].BemGEFc_star)

    ode_list_old[0] = assemble(ode_form1)
    ode_list_old[1] = assemble(ode_form2) + k9_term
    ode_list_old[2] = assemble(ode_form3) - k9_term
    k2 = dt * ode_list_old
    
    x_new = x_old + 0.5*(k1 + k2)

    # Cdc42c
    exp_list[0].Cdc42c = x_new[0] 
    # BemGEFc
    exp_list[1].BemGEFc = x_new[1]
    # BemGEFc_star
    exp_list[2].BemGEFc_star = x_new[2]


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
    
    Cdc42T_tmp, Cdc42D_tmp, BemGEF42_tmp, BemGEFm_tmp, BemGEF42_star_tmp, BemGEFm_star_tmp = U_high.split(True)
    find_pole_data = calc_find_pole_data(mesh, Cdc42T_tmp)
    n_poles = Common.calc_n_poles_her(find_pole_data, eps_clust, min_samples_clust)
    
    return n_poles


def calc_cdc42_max(U_high):
    Cdc42T_tmp, Cdc42D_tmp, BemGEF42_tmp, BemGEFm_tmp, BemGEF42_star_tmp, BemGEFm_star_tmp = U_high.split(True)
    return Cdc42T_tmp.vector().max()



def calc_ode_vals(U, mod_param, exp_list, eta_div_cell_area, dx):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star = U.split()
    Cdc42c_exp, BemGEFc_exp, BemGEFc_star_exp = exp_list
    
    x_old = np.zeros(3)
    # Cdc42 
    x_old[0] = exp_list[0].Cdc42c
    # BemGEFc
    x_old[1] = exp_list[1].BemGEFc
    # BemGEFc_start
    x_old[2] = exp_list[2].BemGEFc_star

    Cdc42c_exp, BemGEFc_exp, BemGEFc_star_exp = exp_list

    ode_form1, ode_form2, ode_form3 = 0.0, 0.0, 0.0

    k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, k8max, k8h, k8n, k9max, k9h, k9n, Dm_gdp = mod_param.return_param()    
    eta_div_cell_area = Constant(eta_div_cell_area)

    k9 = k9max * (exp_list[2].BemGEFc_star**k9n) / (k9h**k9n + exp_list[2].BemGEFc_star**k9n)
    k9_term = Constant(k9 * exp_list[2].BemGEFc_star)

    ode_form1 += eta_div_cell_area * (k5b*Cdc42D - k5a*Cdc42c_exp) * dx(1)
    ode_form2 += eta_div_cell_area * (k1b * BemGEFm - k1a * BemGEFc_exp - k7*BemGEFc_exp*Cdc42T) * dx(1) 
    ode_form3 += eta_div_cell_area * (k1b * BemGEFm_star - k1a * BemGEFc_star_exp - k7*BemGEFc_star_exp*Cdc42T) *dx(1)

    ode_list_old = np.zeros(3)
    # Cdc42
    ode_list_old[0] = assemble(ode_form1)
    # BemGEFc
    ode_list_old[1] = assemble(ode_form2) + k9_term
    # BemGEFs_star 
    ode_list_old[2] = assemble(ode_form3) - k9_term

    return ode_list_old


def compute_save_pole_data(U_high, exp_list_high, mesh, mod_param, eta_div_cell_area, dx):
    
    ode1, ode2, ode3 = calc_ode_vals(U_high, mod_param, exp_list_high, eta_div_cell_area, dx)
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star = U_high.split()
    Cdc42T_surface = assemble(Cdc42T * dx(1))
    Cdc42D_surface = assemble(Cdc42D * dx(1))
    BemGEF42_surface = assemble(BemGEF42 * dx(1))
    BemGEFm_surface = assemble(BemGEFm * dx(1))
    BemGEF42_star_surface = assemble(BemGEF42_star * dx(1))
    BemGEFm_star_surface = assemble(BemGEFm_star * dx(1))

    Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star = U_high.split(True)
    Cdc42c_exp, BemGEFc_exp, BemGEFc_star_exp = exp_list_high

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
    BemGEF42_star_val = BemGEF42_star.vector().get_local()
    BemGEFm_star_val = BemGEFm_star.vector().get_local()

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
                              "BemGEFm" : BemGEFm_val, 
                              "BemGEF42_star" : BemGEF42_star_val, 
                              "BemGEFm_star" : BemGEFm_star_val})
    # Filter to keep the pole 
    data_pole = data_pole[index_vec]

    data_pole["Cdc42c"] = [Cdc42c_exp.Cdc42c for j in range(len(data_pole["x"]))]
    data_pole["BemGEFc"] = [BemGEFc_exp.BemGEFc for j in range(len(data_pole["x"]))]
    data_pole["BemGEFc_star"] = [BemGEFc_star_exp.BemGEFc_star for j in range(len(data_pole["x"]))]

    data_pole["Cdc42T_surface"] = [Cdc42T_surface for j in range(len(data_pole["x"]))]
    data_pole["Cdc42D_surface"] = [Cdc42D_surface for j in range(len(data_pole["x"]))]
    data_pole["BemGEF42_surface"] = [BemGEF42_surface for j in range(len(data_pole["x"]))]
    data_pole["BemGEFm_surface"] = [BemGEFm_surface for j in range(len(data_pole["x"]))]
    data_pole["BemGEF42_star_surface"] = [BemGEF42_star_surface for j in range(len(data_pole["x"]))]
    data_pole["BemGEFm_star_surface"] = [BemGEFm_star_surface for j in range(len(data_pole["x"]))]

    data_pole["Ode1"] = [ode1 for j in range(len(data_pole["x"]))]
    data_pole["Ode2"] = [ode2 for j in range(len(data_pole["x"]))]
    data_pole["Ode3"] = [ode3 for j in range(len(data_pole["x"]))]
    
    return data_pole 


def compute_residual_form(mod_param, trial_f, trial_f_prev, test_f, measures, dt_use, exp_list):

    Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star = trial_f
    Cdc42T_prev, Cdc42D_prev, BemGEF42_prev, BemGEFm_prev, BemGEF42_star_prev, BemGEFm_star_prev = trial_f_prev
    phi1, phi2, phi3, phi4, phi5, phi6 = test_f
    dx = measures

    Cdc42c_exp, BemGEFc_exp, BemGEFc_star_exp = exp_list
    
    k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, k8max, k8h, k8n, k9max, k9h, k9n, Dm_gdp = mod_param.make_param_constants()

    Cdc42T_mass = (Cdc42T - Cdc42T_prev)*phi1 * dx(1)
    Cdc42D_mass = (Cdc42D - Cdc42D_prev)*phi2 * dx(1)
    BemGEF42_mass = (BemGEF42 - BemGEF42_prev)*phi3 * dx(1)
    BemGEFm_mass = (BemGEFm - BemGEFm_prev)*phi4 * dx(1)
    BemGEF42_star_mass = (BemGEF42_star - BemGEF42_star_prev)*phi5 * dx(1)
    BemGEFm_star_mass = (BemGEFm_star - BemGEFm_star_prev)*phi6 * dx(1)
    
    Cdc42T_stiff = Dm * dt_use * dot(grad(Cdc42T), grad(phi1)) * dx(1)
    Cdc42D_stiff = Dm_gdp * dt_use * dot(grad(Cdc42D), grad(phi2)) * dx(1)
    BemGEF42_stiff = Dm * dt_use * dot(grad(BemGEF42), grad(phi3)) * dx(1)
    BemGEFm_stiff = Dm * dt_use * dot(grad(BemGEFm), grad(phi4)) * dx(1)
    BemGEF42_star_stiff = Dm * dt_use * dot(grad(BemGEF42_star), grad(phi5)) * dx(1)
    BemGEFm_star_stiff = Dm * dt_use * dot(grad(BemGEFm_star), grad(phi6)) * dx(1)

    # Regulate strength of nf by altering the maximum phosphorylation rate 
    strength_feedback = mod_param.strength_feedback

    # Helper termks 
    k8 = k8max * (BemGEF42 + BemGEF42_star)**k8n / (k8h**k8n + (BemGEF42 + BemGEF42_star)**k8n)
    k8 *= strength_feedback
    k7_term = k7 * (BemGEFc_star_exp  + BemGEFc_exp)
    k4a_term = k4a * (BemGEFm + BemGEFm_star)
    k4b_term = k4b * (BemGEF42 + BemGEF42_star)
    k8_term = k8 * (BemGEF42 + BemGEF42_star) * BemGEF42
    # Cdc42T
    f1 = (k2a*BemGEFm + k3*BemGEF42) * Cdc42D - (k2b + k4a_term + k7_term)*Cdc42T + k4b_term        
    # Cdc42D
    f2 = k2b*Cdc42T - (k2a*BemGEFm + k3*BemGEF42)*Cdc42D - k5b*Cdc42D + k5a*Cdc42c_exp
    # BemGEF42
    f3 = (k4a*BemGEFm + k7*BemGEFc_exp)*Cdc42T - k4b*BemGEF42 - k8_term
    # BemGEF
    f4 = k1a*BemGEFc_exp - k1b*BemGEFm + k4b*BemGEF42 - k4a*BemGEFm * Cdc42T
    # BemGEF42*
    f5 = (k4a * BemGEFm_star + k7*BemGEFc_star_exp)*Cdc42T - k4b*BemGEF42_star + k8_term
    # BemGEF*
    f6 = k1a*BemGEFc_star_exp - k1b*BemGEFm_star + k4b*BemGEF42_star - k4a*BemGEFm_star * Cdc42T
    
    residual_form = (-f1*phi1 - f2*phi2 - f3*phi3 - f4*phi4 - f5*phi5 - f6*phi6) * dt_use * dx(1)
    residual_form += Cdc42T_mass + Cdc42D_mass + BemGEF42_mass + BemGEFm_mass + BemGEF42_star_mass + BemGEFm_star_mass
    
    return residual_form
    
# dt_list = [dt_low, dt_high]
def perform_time_step(A_low_list, b_low_list, A_high_list, b_high_list, exp_list_high, dt_list,
                      U_high, U_prev_high, U_low, U_norm, U_old, exp_val_old, mod_param, solver, sol_opt, 
                      eta_div_cell_area, dx):

    A_low, A_low_mass, A_low_stiff = A_low_list
    mass_rhs_low, stiffness_reaction_rhs_low = b_low_list

    A_high, A_high_mass, A_high_stiff = A_high_list
    mass_rhs_high, stiffness_reaction_rhs_high = b_high_list

    dt_low, dt_high = dt_list

    while True:
            dt_use = dt_low

            A_low.zero()
            A_low.axpy(1.0, A_low_mass, True)
            A_low.axpy(dt_low, A_low_stiff, True)

            b_low = assemble(mass_rhs_low + Constant(dt_low) * stiffness_reaction_rhs_low)
            solver.solve(A_low,  U_low.vector(), b_low)

            # Update high accuracy system
            A_high.zero()
            A_high.axpy(1.0, A_high_mass, True)
            A_high.axpy(dt_high, A_high_stiff, True)
            
            # High accuracy update t -> t + dt_high
            b_high = assemble(mass_rhs_high + Constant(dt_high)*stiffness_reaction_rhs_high)
            solver.solve(A_high,  U_high.vector(), b_high)
            U_prev_high.assign(U_high)
            solve_ode_system(U_high, mod_param, exp_list_high, dt_high, eta_div_cell_area, dx)

            # High accuracy update t + dt_high -> t + dt_high*2
            b_high = assemble(mass_rhs_high + Constant(dt_high)*stiffness_reaction_rhs_high)
            solver.solve(A_high,  U_high.vector(), b_high)
            U_prev_high.assign(U_high)

            U_norm.vector()[:] = U_low.vector() - U_high.vector()
            dt_low, dt_high, progress = Common.update_step_length(U_norm, dt_low, sol_opt.adaptive_rel_tol, sol_opt.adaptive_rho)
            if progress == False:
                exp_list_high[0].Cdc42c = exp_val_old[0]
                exp_list_high[1].BemGEFc = exp_val_old[1]
                exp_list_high[2].BemGEFc_star = exp_val_old[2]
                U_prev_high.vector()[:] = U_old.vector()
            else:
                # Only update ODE-system if time-step is accepted 
                solve_ode_system(U_high, mod_param, exp_list_high, dt_use / 2.0, eta_div_cell_area, dx)
                break 

    return dt_low, dt_high, dt_use


def perform_time_step_residual(A_high_list, b_high_list, exp_list_high, dt_high, U_high, U_prev_high, 
        mod_param, solver, sol_opt, eta_div_cell_area, dx, residual_form, dt_res):


    A_high, A_high_mass, A_high_stiff = A_high_list
    mass_rhs_high, stiffness_reaction_rhs_high = b_high_list

    while True:
            dt_use = dt_high

            # Update high accuracy system
            A_high.zero()
            A_high.axpy(1.0, A_high_mass, True)
            A_high.axpy(dt_high, A_high_stiff, True)
            
            # High accuracy update t -> t + dt_high
            b_high = assemble(mass_rhs_high + Constant(dt_high)*stiffness_reaction_rhs_high)
            solver.solve(A_high,  U_high.vector(), b_high)

            dt_res.assign(dt_use)
            R = assemble(residual_form)
            norm_error = norm(R, 'l2')
            dt_high, progress = Common.update_step_length_new(norm_error, dt_use, sol_opt.adaptive_rel_tol, sol_opt.adaptive_rho)

            if progress == True:
                # Only update ODE-system if time-step is accepted 
                U_prev_high.assign(U_high)
                solve_ode_system(U_high, mod_param, exp_list_high, dt_use / 2.0, eta_div_cell_area, dx)
                break 

    return dt_high, dt_use


def solve_2d_nf(file_loc, mod_param, sol_opt, exocyte_data):
    
    Config.lock.acquire()
    file_loc.process_msh_file_to_xdmf_2d()
    file_loc.calc_pole_index()
    mesh, subdomains = file_loc.read_mesh_2d()
    Config.lock.release()

    print("r = {:.3e}, k2b = {:.3e}, k2a = {:.3e}".format(mod_param.r, mod_param.k2b, mod_param.k2a))
    
    dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
    
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement([P1, P1, P1, P1, P1, P1])
    H = FunctionSpace(mesh, element)

    if exocyte_data != None:
        exocyte_data.calc_coords(H)
        compute_time_exocyte = True
        t_hit = 0.0
    else:
        compute_time_exocyte = False
        t_hit = np.inf
    
    # Initial values
    if sol_opt.seed != False:
        np.random.seed(sol_opt.seed)
    u0_exp = IC_2d_nf(ic_Cdc42T=mod_param.ic_Cdc42T,
                      ic_Cdc42D=mod_param.ic_Cdc42D,
                      ic_BemGEF42=mod_param.ic_BemGEF42,
                      ic_BemGEFm=mod_param.ic_BemGEFm,
                      ic_BemGEF42_star=mod_param.ic_BemGEF42_star,
                      ic_BemGEFm_star=mod_param.ic_BemGEFm_star,
                      sigma=0.1,
                      element=H.ufl_element())

    # Weak formulation low accuracy Euler solution
    test_functions_low = TestFunctions(H)
    trial_functions_low = TrialFunctions(H)
    U_prev_low = Function(H) 
    U_prev_low.interpolate(u0_exp)
    # Expressions for ODE-system components 
    Cdc42c_exp_low = Expression('Cdc42c', Cdc42c = mod_param.ic_Cdc42I, element=P1)
    BemGEFc_exp_low = Expression('BemGEFc', BemGEFc = mod_param.ic_BemGEFc, element=P1)
    BemGEFc_star_exp_low = Expression('BemGEFc_star', BemGEFc_star = mod_param.ic_BemGEFc_star, element=P1)
    exp_list_low = [Cdc42c_exp_low, BemGEFc_exp_low, BemGEFc_star_exp_low]

    # If saved results are used 
    if sol_opt.read_old_run == True:
        U_vector, exp_list_vector = Common.read_sim_result(file_loc.dir_save_run)
        U_prev_low.vector()[:] = U_vector
        Cdc42c_exp_low.Cdc42c = exp_list_vector[0]
        BemGEFc_exp_low.BemGEFc = exp_list_vector[1]
        BemGEFc_star_exp_low.BemGEFc_star = exp_list_vector[2]

    # The explicit-implicit scheme -> linear system. Here extracting components of said linear system
    mass_low, stiff_low, reaction_low = formulate_components_fem_2d_nf(mod_param, trial_functions_low, U_prev_low, test_functions_low, dx, exp_list_low)
    mass_lhs_low, stiffness_reaction_lhs_low, mass_rhs_low, stiffness_reaction_rhs_low = compute_fem_system_2d_nf(mass_low, stiff_low, reaction_low)
    b_low_list = [mass_rhs_low, stiffness_reaction_rhs_low]
    # Assigning initial values 
    U_low = Function(H)
    U_low.assign(U_prev_low)

    # Weak formulation high accuracy Euler solution
    test_functions_high = TestFunctions(H)
    trial_functions_high = TrialFunctions(H)
    U_prev_high = Function(H) 
    U_prev_high.vector()[:] = U_prev_low.vector()
    # Expressions for ODE-system components 
    Cdc42c_exp_high = Expression('Cdc42c', Cdc42c = Cdc42c_exp_low.Cdc42c, element=P1)
    BemGEFc_exp_high = Expression('BemGEFc', BemGEFc = BemGEFc_exp_low.BemGEFc, element=P1)
    BemGEFc_star_exp_high = Expression('BemGEFc_star', BemGEFc_star = BemGEFc_star_exp_low.BemGEFc_star, element=P1)
    exp_list_high = [Cdc42c_exp_high, BemGEFc_exp_high, BemGEFc_star_exp_high]
    # The explicit-implicit scheme -> linear system. Here extracting components of said linear system
    mass_high, stiff_high, reaction_high = formulate_components_fem_2d_nf(mod_param, trial_functions_high, U_prev_high, test_functions_high, dx, exp_list_high)
    mass_lhs_high, stiffness_reaction_lhs_high, mass_rhs_high, stiffness_reaction_rhs_high = compute_fem_system_2d_nf(mass_high, stiff_high, reaction_high)
    b_high_list = [mass_rhs_high, stiffness_reaction_rhs_high]
    # Assigning initial values 
    U_high = Function(H)
    Cdc42T_high, Cdc42D_high, BemGEF42_high, BemGEFm_high, BemGEF42_star_high, BemGEFm_star_high = split(U_high)
    U_high.assign(U_prev_high)
    _Cdc42T_high, _Cdc42D_high, _BemGEF42_high, _BemGEFm_high, _BemGEF42_star_high, _BemGEFm_star_high = U_high.split()
    # For second time-stepper 
    dt_res = Constant(1.0)
    residual_form = compute_residual_form(mod_param, U_high, U_prev_high, test_functions_high, dx, dt_res, exp_list_high)

    # Writing states to pvd-files
    t_curr = 0.0
    vtkfile_Cdc42T = File(file_loc.dir_pvd + "Cdc42T" + ".pvd")
    vtkfile_BemGEFm_star = File(file_loc.dir_pvd + "BemGEFm_star" + ".pvd")
    vtkfile_BemGEFm = File(file_loc.dir_pvd + "BemGEFm" + ".pvd")
    vtkfile_BemGEF42 = File(file_loc.dir_pvd + "BemGEF42" + ".pvd")    
    vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
    vtkfile_BemGEFm << (_BemGEFm_high, t_curr)
    vtkfile_BemGEFm_star << (_BemGEFm_star_high, t_curr)
    vtkfile_BemGEF42 << (_BemGEF42_high, t_curr)
    
    # Parameters for efficient time-stepping 
    eta_div_cell_area = mod_param.eta / (4*np.pi*(mod_param.r / mod_param.L)**2) # For ODE:s
    dt_low = 1e-3
    dt_high = 1e-3 / 2
    t_it = 1
    # Parameters and variables for adaptive time solver
    U_norm = Function(H)
    U_old = Function(H)
    U_old.vector()[:] = U_prev_high.vector()

    solver = KrylovSolver("gmres", "ilu")
    #solver.parameters["nonzero_initial_guess"] = True
    exp_val_old = np.array([Cdc42c_exp_high.Cdc42c, BemGEFc_exp_high.BemGEFc, BemGEFc_star_exp_high.BemGEFc_star])

    parameters['form_compiler']['optimize'] = True
    parameters['form_compiler']['cpp_optimize'] = True
    
    # Only assemble left-hand size matrices once (speed-up)
    mass_lhs_low = assemble(mass_lhs_low)
    stiffness_reaction_lhs_low = assemble(stiffness_reaction_lhs_low)
    mass_lhs_high = assemble(mass_lhs_high)
    stiffness_reaction_lhs_high = assemble(stiffness_reaction_lhs_high)

    # For more efficient matrix calculations 
    A_low = mass_lhs_low + dt_low * stiffness_reaction_lhs_low
    A_low_mass = 1.0 * mass_lhs_low + 0.0 * stiffness_reaction_lhs_low
    A_low_stiff = 0.0 * mass_lhs_low + 1.0 * stiffness_reaction_lhs_low
    A_low_list = [A_low, A_low_mass, A_low_stiff]
    A_high = mass_lhs_high + dt_low * stiffness_reaction_lhs_high
    A_high_mass = 1.0 * mass_lhs_high + 0.0 * stiffness_reaction_lhs_high
    A_high_stiff = 0.0 * mass_lhs_high + 1.0 * stiffness_reaction_lhs_high
    A_high_list = [A_high, A_high_mass, A_high_stiff]

    dt_use = 0.0

    # For termination critera 
    cdc42_limit_check_term = mod_param.ic_Cdc42T * sol_opt.term_inc_cdc42T
    pole_ratio_strict_old = 0.0
    cdc42_max_best = 0.0
    cdc42_max_old = 0.0
    term_crit = "Max_it"
    n_times_flat_max = 0
    n_times_flat_area = 0    

    U_high.set_allow_extrapolation(True)

    # Add actin cables if used 
    if exocyte_data != None and exocyte_data.use_cables == True:
        Cdc42T_tmp, tmp1, tmp2, tmp3, tmp4, tmp5 = U_high.split(True)
        for i in range(exocyte_data.n_cables):
            exocyte_data.calc_cable_attach(Cdc42T_tmp, i, n_cables_use=i)

    # If location of cables are plotted 
    if exocyte_data != None:
        if exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
            exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

    for i in tqdm(range(sol_opt.term_max_it)):

        if sol_opt.solver == "two_step":
            dt_low, dt_high, dt_use = perform_time_step(A_low_list, b_low_list, A_high_list, b_high_list, exp_list_high, 
                    [dt_low, dt_high], U_high, U_prev_high, U_low, U_norm, U_old, exp_val_old, mod_param, solver, sol_opt, 
                    eta_div_cell_area, dx)
        else:
            dt_high, dt_use = perform_time_step_residual(A_high_list, b_high_list, exp_list_high, dt_high, U_high, U_prev_high, 
                mod_param, solver, sol_opt, eta_div_cell_area, dx, residual_form, dt_res)

        if exocyte_data != None and compute_time_exocyte == True:
            t_hit = exocyte_data.calc_time_hit()
            compute_time_exocyte = False
            t_hit = t_curr + t_hit
            print("Time-hit = {}".format(t_hit))
            
        t_curr += dt_use
        t_it += 1
        
        # Set low accuracy solution at high accuracy solution

        # Check if exocytosis time should be computed 
        if t_curr > t_hit:
            print("Exocytosis occuring")
            Cdc42T_tmp, tmp1, tmp2, tmp3, tmp4, tmp5 = U_high.split(True)
            if exocyte_data.use_cables == False:
                i_hit = exocyte_data.calc_site_hit(Cdc42T_tmp)
            else:
                i_hit = exocyte_data.calc_site_hit_cable()

            coord_hit_sphere = exocyte_data.get_coord_sphere_hit(i_hit)
            for j in range(exocyte_data.n_coord):
                exocyte_data.calc_coord_interpolate(j, coord_hit_sphere)
            
            exocyte_data.set_values_post_exocyte(U_high)
            U_high.vector()[:] = exocyte_data.U_new.vector()
            U_prev_high.assign(U_high)
            compute_time_exocyte = True

            vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
            # If location of cables are plotted 
            if exocyte_data != None:
                if exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
                    exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

            # Correct the cytosolic diffusion 
            correct_cdc42c = mod_param.Cdc42_tot_pract - assemble((Cdc42T_high + Cdc42D_high + BemGEF42_high + BemGEF42_star_high)*dx(1)) *eta_div_cell_area
            correct_BemGEFc = mod_param.BemGEF_tot_pract - assemble((BemGEFm_high + BemGEF42_high)*dx(1)) *eta_div_cell_area
            correct_BemGEFc_star = mod_param.BemGEF_star_tot_pract - assemble((BemGEFm_star_high + BemGEF42_star_high)*dx(1)) *eta_div_cell_area

            exp_list_high[0].Cdc42c = correct_cdc42c
            exp_list_high[1].BemGEFc = correct_BemGEFc
            exp_list_high[2].BemGEFc_star = correct_BemGEFc_star

        # Calculate new locations for cables if required 
        if exocyte_data != None and exocyte_data.use_cables == True:
            Cdc42T_tmp, tmp1, tmp2, tmp3, tmp4, tmp5 = U_high.split(True)
            exocyte_data.update_cable_pos(t_curr, Cdc42T_tmp)

        # Set low accuracy solution at high accuracy solution
        exp_val_old[0] = exp_list_high[0].Cdc42c
        exp_val_old[1] = exp_list_high[1].BemGEFc
        exp_val_old[2] = exp_list_high[2].BemGEFc_star
        Cdc42c_exp_low.Cdc42c = exp_list_high[0].Cdc42c
        BemGEFc_exp_low.BemGEFc = exp_list_high[1].BemGEFc
        BemGEFc_star_exp_low.BemGEFc_star = exp_list_high[2].BemGEFc_star
        U_prev_low.vector()[:] = U_high.vector()
        U_old.vector()[:] = U_high.vector()
        
        U_prev_low.vector()[:] = U_high.vector()
        U_old.vector()[:] = U_high.vector()
        
        # Export solution to disk 
        if t_it % sol_opt.print_pwd_it == 0:
            vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
            vtkfile_BemGEFm << (_BemGEFm_high, t_curr)
            vtkfile_BemGEFm_star << (_BemGEFm_star_high, t_curr)
            vtkfile_BemGEF42 << (_BemGEF42_high, t_curr)
            if exocyte_data != None:
                if exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
                    exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

        # If saving run-result 
        if sol_opt.save_run:    
            exp_list_vector = np.array([exp_list_high[0].Cdc42c, exp_list_high[1].BemGEFc, exp_list_high[2].BemGEFc_star])
            Common.save_sim_result(U_high, exp_list_vector, file_loc.dir_save_run)

        if t_it % sol_opt.print_pwd_it == 0 and t_curr > 0.0 and sol_opt.save_dist:
            data_pole = compute_save_pole_data(U_high, exp_list_high, mesh, mod_param, eta_div_cell_area, dx)
            file_name = file_loc.dir_intermediate + "Pole_dist_t.csv"
            n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
            Common.save_pole_dist_data(data_pole, t_curr, file_loc, file_name, n_poles)

        # If we want to export data during simulations 
        if sol_opt.save_data_t and t_it % sol_opt.save_data_t_it == 0:
            Cdc42T_high, Cdc42D_high, BemGEF42_high, BemGEFm_high, BemGEF42_star_high, BemGEFm_star_high = U_high.split(True)
            n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
            pole_data = Common.PoleData(Cdc42T_high, mesh, t_curr, n_poles, term_crit="Not")
            pole_data.save_data(file_loc.dir_intermediate + "Pole_data_t.csv", file_loc)

        # Check if exocyte hit should occur, and update actin cables (if present)

        # Termination algorithm
        cdc42_max = calc_cdc42_max(U_high)
        if cdc42_max > cdc42_limit_check_term and sol_opt.term_check and t_curr > 50.0:
            
            # See if concentration flattens out or if concentration has peaked
            if cdc42_max > cdc42_max_best:
                cdc42_max_best = cdc42_max
                n_times_smaller = 0
            else:
                n_times_smaller += 1

            # Check if we are 0.8 below the max 
            below_max = False
            if cdc42_max < cdc42_max_best * 0.8:
                below_max = True
            
            # For flattnes, better to look at pole-ratio 
            Cdc42T_high_tmp, x1, x2, x3, x4, x5 = U_high.split(True)
            pole_ratio, pole_ratio_strict, cdc42_max = Common.calc_pole_ratios(Cdc42T_high_tmp)

            if abs(pole_ratio_strict - pole_ratio_strict_old) > (sol_opt.term_tol_flat):
                n_times_flat_area = 0
            else:
                n_times_flat_area += 1

            if abs(cdc42_max - cdc42_max_old) > (sol_opt.term_tol_flat_max):
                n_times_flat_max = 0
            else:
                n_times_flat_max += 1
            
            print("n_times_flat = {}".format(n_times_flat_area))
            
            if n_times_flat_max > sol_opt.term_times_flat and n_times_flat_area > sol_opt.term_times_flat and below_max:
                n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
                print("n_poles = {}".format(n_poles))
                if n_poles == 1:
                    print("Terminates upon a single pole")
                    print("n_poles = {}".format(n_poles))
                    term_crit = "One_pole"
                    break
            
            pole_ratio_strict_old = pole_ratio_strict
            cdc42_max_old = cdc42_max

    # Write end coordinates
    vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
    vtkfile_BemGEFm << (_BemGEFm_high, t_curr)
    vtkfile_BemGEFm_star << (_BemGEFm_star_high, t_curr)
    vtkfile_BemGEF42 << (_BemGEF42_high, t_curr)
    n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)

    if sol_opt.save_dist:
        Config.lock.acquire()
        try:
            data_pole = compute_save_pole_data(U_high, exp_list_high, mesh, mod_param, eta_div_cell_area, dx)
            file_name = file_loc.dir_intermediate + "Pole_dist_end.csv"
            n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
            Common.save_pole_dist_data(data_pole, t_curr, file_loc, file_name, n_poles)
        except:
            print("Simulation failed saving pole data")
        Config.lock.release()
    
    # Calculate pole data and export to disk 
    Config.lock.acquire()
    try:
        Cdc42T_high, Cdc42D_high, BemGEF42_high, BemGEFm_high, BemGEF42_star_high, BemGEFm_star_high = U_high.split(True)
        pole_data = Common.PoleData(Cdc42T_high, mesh, t_curr, n_poles, term_crit="Max_it")
        pole_data.save_data(file_loc.dir_intermediate + "Pole_data.csv", file_loc)
        n_poles = calc_n_poles_2d(mesh, U_high, sol_opt.term_eps_clust, sol_opt.term_min_samp_clust)
        pole_data.print_f()
        print("n_poles = {:.0f}".format(n_poles))
    except:
        print("Simulation failed")
    Config.lock.release()


def run_2d_model_nf(file_loc, mod_param_list, sol_opt_list, times_run, n_threads, exocyte_data=None):
    
    input_list = []
    for mod_param in mod_param_list:
        for sol_opt in sol_opt_list:
            file_loc_tmp = Common.FileLocations(file_loc.mesh_name, file_loc.save_name, mod_param.calc_save_tag(), name_store=file_loc.name_store)
            input_list.append([copy.deepcopy(file_loc_tmp), copy.deepcopy(mod_param), copy.deepcopy(sol_opt), 
            copy.deepcopy(exocyte_data)])
    
    len_input_list = len(input_list)
    for i in range(times_run-1):
        for j in range(len_input_list):
            input_list.append(copy.deepcopy(input_list[j]))
    
    print("Length of input argument list = {}".format(len(input_list)))
    
    print("Starting simulations for 2d system nf")
    l = mp.Lock()
    pool = mp.Pool(initializer=init, initargs=(l,), processes=n_threads)
    result = pool.starmap(solve_2d_nf, input_list)
    pool.close()
    pool.join()
    print("Done with simulations for 2d system nf")
