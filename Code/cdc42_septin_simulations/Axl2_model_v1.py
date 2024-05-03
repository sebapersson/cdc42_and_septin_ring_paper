#!/usr/bin/env python

from dolfin import *
from tqdm import tqdm
import numpy as np
import pandas as pd
import os
import copy 
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
        self.ic_P = kwargs.pop('ic_P')
        self.ic_GapS = kwargs.pop('ic_GapS')
        self.ic_S = kwargs.pop('ic_S')
        self.ic_Axl2 = kwargs.pop('ic_Axl2')
        self.ic_Axl2_S = kwargs.pop('ic_Axl2_S')
        super(IC_2d, self).__init__(*args, **kwargs)
    
    def eval(self, values, x):
        
        values[0] = self.ic_Cdc42T + np.random.normal(scale=0.001)
        values[1] = self.ic_Cdc42D + np.random.normal(scale=0.001)
        values[2] = self.ic_BemGEF42 + np.random.normal(scale=0.001)
        values[3] = self.ic_BemGEFm + np.random.normal(scale=0.001)
        values[4] = self.ic_S + 0.0
        values[5] = self.ic_P + 0.0
        values[6] = self.ic_GapS + 0.0    
        values[7] = self.ic_Axl2 + 0.0
        values[8] = self.ic_Axl2_S + 0.0
    
    def value_shape(self):
        return(9,)        


def formulate_components_fem_2d(mod_param, trial_f, trial_f_prev, test_f, measures, mesh, exp_list):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm, S, P, GapS, Axl2, Axl2_S = trial_f
    Cdc42T_prev, Cdc42D_prev, BemGEF42_prev, BemGEFm_prev, S_prev, P_prev, GapS_prev, Axl2_prev, Axl2_S_prev = trial_f_prev
    phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9 = test_f
    Cdc42c_exp, BemGEFc_exp, Sc_exp, GapSc_exp, Axl2c_exp = exp_list

    dx = measures
    
    mass_form, stiff_form, reaction_form = 0, 0, 0
    k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, k12a, k12b, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, Dm, Dms, Dmss_val, Dc, eta, Dm_gdp = mod_param.make_param_constants()

    Cdc42T_mass = (Cdc42T - Cdc42T_prev)*phi1 * dx
    Cdc42D_mass = (Cdc42D - Cdc42D_prev)*phi2 * dx
    BemGEF42_mass = (BemGEF42 - BemGEF42_prev)*phi3 * dx
    BemGEFm_mass = (BemGEFm - BemGEFm_prev)*phi4 * dx
    S_mass = (S - S_prev) * phi5 * dx
    P_mass = (P - P_prev) * phi6 * dx
    GapS_mass = (GapS - GapS_prev) * phi7 * dx
    Axl2_mass = (Axl2 - Axl2_prev) * phi8 * dx
    Axl2_S_mass = (Axl2_S - Axl2_S_prev) * phi9 * dx

    if mod_param.use_gamma == True:
        gamma_d = mod_param.gamma_d
        gamma_k = mod_param.gamma_k 
        Dmss = Dmss_val / (1.0 + P_prev / gamma_d)
        k18_use = k18 / (1.0 + P_prev / gamma_k)
    else:
        k18_use = k18
        Dmss = Dmss_val

    if mod_param.crowding_p == True:
        k16 = conditional(lt(P_prev, 200.0), k16, 0.0)
        k17 = conditional(lt(P_prev, 200.0), k17, 0.0)
    
    Cdc42T_stiff = Dm * dot(grad(Cdc42T), grad(phi1)) * dx
    Cdc42D_stiff = Dm_gdp * dot(grad(Cdc42D), grad(phi2)) * dx
    BemGEF42_stiff = Dm * dot(grad(BemGEF42), grad(phi3)) * dx
    BemGEFm_stiff = Dm * dot(grad(BemGEFm), grad(phi4)) * dx
    S_stiff = Dms * dot(grad(S), grad(phi5)) * dx
    P_stiff = Dmss * dot(grad(P), grad(phi6)) * dx
    GapS_stiff = Dms * dot(grad(GapS), grad(phi7)) * dx
    Axl2_stiff = Dm * dot(grad(Axl2), grad(phi8)) * dx
    Axl2_S_stiff = Dms * dot(grad(Axl2_S), grad(phi9)) * dx

    # Cdc42T
    f1 = (k2a*BemGEFm_prev + k3*BemGEF42_prev) * Cdc42D_prev - (k2b + k4a*BemGEFm_prev + k7*BemGEFc_exp + k13*GapS_prev + k14*GapSc_exp)*Cdc42T_prev + k4b*BemGEF42
        
    # Cdc42D
    f2 = (k2b + k13*GapS_prev + k14*GapSc_exp)*Cdc42T_prev - (k2a*BemGEFm_prev + k3*BemGEF42_prev)*Cdc42D_prev - k5b*Cdc42D + k5a*Cdc42c_exp
    
    # BemGEF42
    f3 = (k4a*BemGEFm_prev + k7*BemGEFc_exp)*Cdc42T_prev - k4b*BemGEF42

    # BemGEF
    f4 = k1a*BemGEFc_exp - k1b*BemGEFm + k4b*BemGEF42 - k4a*BemGEFm_prev*Cdc42T_prev

    # Septin S
    if mod_param.k17_alt == False:
        k17_term = conditional(lt(P_prev, 25.0), 0.0, 1.0) * S_prev * k17
    else:
        k17_term = k17 * P_prev * S_prev
    
    f5 = k15*Axl2_S - (2*k16*S_prev*S_prev + k17_term) + k18_use*P - k21*S - k19*Axl2_prev*S_prev

    # P 
    f6 = 2*k16*S_prev*S_prev + k17_term - k18_use*P

    # GapS
    f7 = k12a*GapSc_exp*(P_prev) - k12b*GapS

    # Axl2
    f8 = k15*Axl2_S - k19*Axl2_prev*S_prev - k20*Axl2_prev*Sc_exp - k22*Axl2
    if mod_param.recruited == True:
        if mod_param.use_k23_exp == True:
            k23 = mod_param.k23_exp
        else:
            k23 = Constant(mod_param.k23)
        f8 = k15*Axl2_S - k19*Axl2_prev*S_prev - k20*Axl2_prev*Sc_exp - k22*Axl2 + k23*BemGEF42_prev*Axl2c_exp

    # Axl2_S
    f9 = k19*Axl2_prev*S_prev + k20*Axl2_prev*Sc_exp - k15*Axl2_S
    
    reaction_form += (-f1*phi1 - f2*phi2 - f3*phi3 - f4*phi4 - f5*phi5 - f6*phi6 - f7*phi7 - f8*phi8 - f9*phi9) * dx
    mass_form += Cdc42T_mass + Cdc42D_mass + BemGEF42_mass + BemGEFm_mass + S_mass + P_mass + GapS_mass + Axl2_mass + Axl2_S_mass
    stiff_form += Cdc42T_stiff + Cdc42D_stiff + BemGEF42_stiff + BemGEFm_stiff + S_stiff + P_stiff + GapS_stiff + Axl2_stiff + Axl2_S_stiff
    
    return mass_form, stiff_form, reaction_form

def compute_fem_system_2d(mass_form, stiff_form, reaction_form):
    
    # Forms for vectors on RHS
    mass_form_rhs = rhs(mass_form)
    stiffness_reaction_form_rhs = rhs(stiff_form + reaction_form)
    
    mass_form_lhs = lhs(mass_form)
    stiffness_reaction_form_lhs = lhs(stiff_form + reaction_form)
    
    return mass_form_lhs, stiffness_reaction_form_lhs, mass_form_rhs, stiffness_reaction_form_rhs


def solve_ode_system(U, mod_param, exp_list, dt, eta_div_cell_area, dx):
    
    Cdc42T, Cdc42D, BemGEF42, BemGEFm, S, P, GapS, Axl2, Axl2_S = U.split()
    Cdc42c_exp, BemGEFc_exp, Sc_exp, GapSc_exp, Axl2c_exp = exp_list
    x_old = np.zeros(5)
    
    # Cdc42 
    x_old[0] = exp_list[0].Cdc42c
    # BemGEFc
    x_old[1] = exp_list[1].BemGEFc
    # Sc
    x_old[2] = exp_list[2].Sc
    # GapS
    x_old[3] = exp_list[3].GapSc
    # Axl2
    x_old[4] = exp_list[4].Axl2c

    k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, k12a, k12b, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, Dm, Dms, Dmss, Dc, eta, Dm_gdp = mod_param.return_param()
    eta_div_cell_area = Constant(eta_div_cell_area)

    ode_form1 = (eta_div_cell_area * (k5b*Cdc42D - k5a*Cdc42c_exp))*dx
    ode_form2 = (eta_div_cell_area * (k1b*BemGEFm - k1a*BemGEFc_exp - k7*BemGEFc_exp*Cdc42T))*dx
    ode_form3 = (eta_div_cell_area * (k21*S - k20*Axl2*Sc_exp))*dx
    ode_form4 = (eta_div_cell_area * (-k12a * GapSc_exp * (P) + k12b * GapS)) * dx
    ode_form5 = (eta_div_cell_area * (k22*Axl2))*dx

    if mod_param.recruited == True:
        if mod_param.use_k23_exp == True:
            k23 = mod_param.k23_exp 
        else:
            k23 = mod_param.k23
        ode_form5 = (eta_div_cell_area * (k22*Axl2 - k23*BemGEF42*Axl2c_exp))*dx
    
    ode_list_old = np.zeros(5)
    ode_list_old[0] = assemble(ode_form1)
    ode_list_old[1] = assemble(ode_form2)
    ode_list_old[2] = assemble(ode_form3)
    ode_list_old[3] = assemble(ode_form4)
    ode_list_old[4] = assemble(ode_form5)
    
    # k1 RK2 method 
    k1 = dt * ode_list_old

    # k2 RK2 method 
    exp_list[0].Cdc42c += k1[0] * 0.5
    exp_list[1].BemGEFc += k1[1] * 0.5
    exp_list[2].Sc += k1[2] * 0.5
    exp_list[3].GapSc += k1[3] * 0.5
    exp_list[4].Axl2c += k1[4] * 0.5
    ode_list_old[0] = assemble(ode_form1)
    ode_list_old[1] = assemble(ode_form2)
    ode_list_old[2] = assemble(ode_form3)
    ode_list_old[3] = assemble(ode_form4)
    ode_list_old[4] = assemble(ode_form5)
    
    k2 = dt * ode_list_old
    
    x_new = x_old + 0.5*(k1 + k2)

    # Ensure positivity for ODE-system 
    x_new[x_new < 0] = 0.0
    exp_list[0].Cdc42c = x_new[0]    
    exp_list[1].BemGEFc = x_new[1]
    exp_list[2].Sc = x_new[2]
    exp_list[3].GapSc = x_new[3]
    exp_list[4].Axl2c = x_new[4]


def compute_mass(U, exp_list, eta_div_cell_area, dx):
    Cdc42T, Cdc42D, BemGEF42, BemGEFm, S, P, GapS, Axl2, Axl2_S = U.split(True)
    Cdc42c_exp, BemGEFc_exp, Sc_exp, GapSc_exp, Axl2c_exp = exp_list

    mass_Cdc42 = Cdc42c_exp.Cdc42c + assemble((Cdc42T + Cdc42D + BemGEF42)*dx) *eta_div_cell_area
    mass_bem = BemGEFc_exp.BemGEFc + assemble((BemGEFm + BemGEF42)*dx) *eta_div_cell_area
    mass_S = Sc_exp.Sc + assemble((Axl2_S + S + P)*dx) *eta_div_cell_area
    mass_Axl2 = Axl2c_exp.Axl2c + assemble((Axl2_S + Axl2)*dx) *eta_div_cell_area
    mass_GapS = GapSc_exp.GapSc + assemble((GapS)*dx) *eta_div_cell_area

    print()
    print("mass_S = {:.3e}".format(mass_S))
    print("mass_Cdc42 = {:.3e}".format(mass_Cdc42))
    print("mass_bem = {:.3e}".format(mass_bem))
    print("mass_axl2 = {:.3e}".format(mass_Axl2))
    print("mass_gapS = {:.3e}".format(mass_GapS))
    

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
            exp_list[2].Sc = exp_val_old[2]
            exp_list[3].GapS = exp_val_old[3]
            exp_list[4].Axl2c = exp_val_old[4]
            U_prev.vector()[:] = U_old.vector()
        else:
            # Only update ODE-system if time-step is accepted
            # First ensure positive PDE-states 
            Uvector = as_backend_type(U_high.vector()).get_local()
            Uvector[Uvector <= DOLFIN_EPS] = 0.0
            U_high.vector().set_local(Uvector)

            solve_ode_system(U_high, mod_param, exp_list, dt_use / 2.0, eta_div_cell_area, dx)
            break 


    return dt_low, dt_high, dt_use


def calc_max_coord_cdc42(U_high, H):

    Cdc42T_tmp, tmp1, tmp2, tmp3, tmp4, P_tmp, tmp5, tmp6, tmp7 = U_high.split(True)
    n_states = 9
    Cdc42T_vec = Cdc42T_tmp.vector().get_local()
    i_max = np.argmax(Cdc42T_vec)
    coord_eu = H.tabulate_dof_coordinates()
    coord_eu = coord_eu[range(0, len(coord_eu), n_states)] # Want to avoid duplicating 
    x_c, y_c, z_c = coord_eu[i_max]
    
    return x_c, y_c, z_c  


def compute_save_dist_data(U_high, mesh, t, file_loc):
    
    print("Saving distribution data")
    Cdc42T, Cdc42D, BemGEF42, BemGEFm, S, P, GapS, Axl2, Axl2_S = U_high.split(True)

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
    S_val = S.vector().get_local()
    P_val = P.vector().get_local()
    GapS_val = GapS.vector().get_local()
    Axl2_val =Axl2.vector().get_local()
    Axl2_S_val = Axl2_S.vector().get_local()

    data_dist = pd.DataFrame({"x": dof_x,
                              "y": dof_y,
                              "z": dof_z,
                              "Cdc42T": Cdc42T_val, 
                              "Cdc42D" : Cdc42D_val, 
                              "BemGEF42" : BemGEF42_val, 
                              "BemGEFm" : BemGEFm_val,
                              "S" : S_val, 
                              "P" : P_val, 
                              "GapS" : GapS_val, 
                              "Axl2" : Axl2_val, 
                              "Axl2_S_val" : Axl2_S_val,
                              "t" : Axl2_S_val*0.0 + t,
                              "index" : Axl2_S_val*0.0 + file_loc.run_index})

    return data_dist                              


def solve_2d(file_loc, mod_param, sol_opt, bud_scar_param, exocyte_data):
    
    Config.lock.acquire()
    file_loc.process_msh_file_to_xdmf_2d()
    file_loc.calc_pole_index()
    mesh, subdomains = file_loc.read_mesh_2d()
    Config.lock.release()
    
    dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
    
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement([P1, P1, P1, P1, P1, P1, P1, P1, P1])
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
    u0_exp = IC_2d(ic_Cdc42T=mod_param.ic_Cdc42T,
                   ic_Cdc42D=mod_param.ic_Cdc42D,
                   ic_BemGEF42=mod_param.ic_BemGEF42,
                   ic_BemGEFm=mod_param.ic_BemGEFm,
                   ic_P=mod_param.ic_P,
                   ic_GapS=mod_param.ic_GapS, 
                   ic_S=mod_param.ic_S,
                   ic_Axl2=mod_param.ic_Axl2,
                   ic_Axl2_S=mod_param.ic_Axl2_S,
                   element=H.ufl_element())

    # If k23 confined to local region 
    if mod_param.use_k23_exp == True:
        x_c, y_c, z_c = 0.0, 0.0, 0.0
        k23 = Expression("pow(x[0]-x_c, 2) + pow(x[1]-y_c, 2) + pow(x[2]-z_c, 2) < tol ? k23 : 0.0", x_c=x_c, y_c=y_c, z_c=z_c, tol=mod_param.tol_k23, k23=mod_param.k23, element=P1)
        mod_param.k23_exp = k23
        

    # Components for low accuracy soluations 
    U_prev = Function(H) 
    U_prev.interpolate(u0_exp)
    exp_list = [Expression('Cdc42c', Cdc42c = mod_param.ic_Cdc42I, element=P1), 
                Expression('BemGEFc', BemGEFc = mod_param.ic_BemGEFc, element=P1), 
                Expression('Sc', Sc = mod_param.ic_Sc, element=P1), 
                Expression('GapSc', GapSc = mod_param.ic_GapSc, element=P1), 
                Expression('Axl2c', Axl2c = mod_param.ic_Axl2c, element=P1)]
    U_low = Function(H)
    U_low.assign(U_prev)
    
    dt_low, dt_high = 1e-5, 1e-3 / 2

    # Weak formulation high accuracy Euler solution 
    test_functions = TestFunctions(H)
    trial_functions = TrialFunctions(H)
    exp_val_old = np.array([exp_list[0].Cdc42c, exp_list[1].BemGEFc, exp_list[2].Sc, exp_list[3].GapSc, exp_list[4].Axl2c])
    # The explicit-implicit scheme -> linear system. Here extracting components of said linear system
    mass, stiff, reaction = formulate_components_fem_2d(mod_param, trial_functions, U_prev, test_functions, dx, mesh, exp_list)
    mass_lhs, stiffness_reaction_lhs, mass_rhs, stiffness_reaction_rhs = compute_fem_system_2d(mass, stiff, reaction)
    b_list = [mass_rhs, stiffness_reaction_rhs]
    mass_lhs_high = assemble(mass_lhs)
    stiffness_reaction_lhs_high = assemble(stiffness_reaction_lhs)
    A = mass_lhs_high + dt_low * stiffness_reaction_lhs_high
    A_mass = 1.0 * mass_lhs_high + 0.0 * stiffness_reaction_lhs_high
    A_stiff = 0.0 * mass_lhs_high + 1.0 * stiffness_reaction_lhs_high
    A_list = [A, A_mass, A_stiff]

    if sol_opt.read_old_run == True:
        if mod_param.read_run != None:
            dir_read = "../../Intermediate/Experiments/Saved_runs/" + mod_param.read_run + "/"
        else:
            dir_read = file_loc.dir_read_run
        U_vector, exp_list_vector = Common.read_sim_result(dir_read)
        print("Dir_read = {}".format(dir_read))
        for i in range(int(len(U_vector)/4)):
            i_start_new = i * 9
            i_start_old = i * 4
            U_prev.vector()[i_start_new] = U_vector[i_start_old]
            U_prev.vector()[i_start_new+1] = U_vector[i_start_old+1]
            U_prev.vector()[i_start_new+2] = U_vector[i_start_old+2]
            U_prev.vector()[i_start_new+3] = U_vector[i_start_old+3]

        exp_list[0].Cdc42c = exp_list_vector[0]
        exp_list[1].BemGEFc = exp_list_vector[1]
        
    # Assigning initial values 
    U_high = Function(H)
    U_high.assign(U_prev)

    # Writing states to pvd-files
    t_curr = 0.0
    _Cdc42T_high, _Cdc42D_high, _BemGEF42_high, _BemGEFm_high, _S_high, _P_high, _GapS_high, _Axl2_high, _Axl2_S_high = U_high.split()
    Cdc42T_high, Cdc42D_high, BemGEF42_high, BemGEFm_high, S_high, P_high, GapS_high, Axl2_high, Axl2_S_high = split(U_high)
    vtkfile_Cdc42T = File(file_loc.dir_pvd + "Cdc42T" + ".pvd")
    vtkfile_P = File(file_loc.dir_pvd + "P" + ".pvd")
    vtkfile_GapS = File(file_loc.dir_pvd + "GapS" + ".pvd")
    vtkfile_S = File(file_loc.dir_pvd + "S" + ".pvd")
    vtkfile_Axl2 = File(file_loc.dir_pvd + "Axl2" + ".pvd")
    vtkfile_Axl2_S = File(file_loc.dir_pvd + "Axl2_S" + ".pvd")
    
    vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
    vtkfile_P << (_P_high, t_curr)    
    vtkfile_GapS << (_GapS_high, t_curr)
    vtkfile_S << (_S_high, t_curr)
    vtkfile_Axl2 << (_Axl2_high, t_curr)
    vtkfile_Axl2_S << (_Axl2_S_high, t_curr)

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

    if sol_opt.term_time != False:
        term_time = sol_opt.term_time
    else:
        term_time = np.inf

    # Add actin cables if used 
    if exocyte_data != None and exocyte_data.use_cables == True and exocyte_data.cables_windows == False:
        Cdc42T_tmp, tmp1, tmp2, tmp3, tmp4, P_tmp, tmp5, tmp6, tmp7 = U_high.split(True)
        exocyte_data.calc_cable_attach_septin(Cdc42T_tmp, P_tmp, i, n_cables_use=i, first=True)
        for i in range(1, exocyte_data.n_cables):
            exocyte_data.calc_cable_attach_septin(Cdc42T_tmp, P_tmp, i, n_cables_use=i)
    if exocyte_data != None and exocyte_data.use_cables == True and exocyte_data.cables_windows == True:
        Cdc42T_tmp, tmp1, tmp2, tmp3, tmp4, P_tmp, tmp5, tmp6, tmp7 = U_high.split(True)
        for i in range(exocyte_data.n_cables):
            exocyte_data.calc_cable_attach_window(Cdc42T_tmp, i)
            
    # If location of cables are plotted 
    if exocyte_data != None and exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
        exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

    # Ensure correct mass if reading old run 
    mod_param.Cdc42_tot_pract = exp_list[0].Cdc42c + assemble((Cdc42T_high + Cdc42D_high + BemGEF42_high)*dx(1)) *eta_div_cell_area
    mod_param.BemGEF_tot_pract = exp_list[1].BemGEFc + assemble((BemGEFm_high + BemGEF42_high)*dx(1)) *eta_div_cell_area

    if mod_param.use_k23_exp == True:
        x_c, y_c, z_c = calc_max_coord_cdc42(U_high, H)
        print("(x_c, y_c, z_c) = ({:.3f}, {:.3f}, {:.3f})".format(x_c, y_c, z_c))
        mod_param.k23_exp.x_c = x_c
        mod_param.k23_exp.y_c = y_c
        mod_param.k23_exp.z_c = z_c

    dt_use = 0.0
    for i in tqdm(range(sol_opt.term_max_it)):

        if t_curr > term_time:
            break

        if exocyte_data != None and compute_time_exocyte == True:
            t_hit = exocyte_data.calc_time_hit()
            compute_time_exocyte = False
            t_hit = t_curr + t_hit
            print("Time-hit = {}".format(t_hit))
        
        dt_low, dt_high, dt_use = perform_time_step(A_list, b_list, exp_list, [dt_low, dt_high], 
                U_high, U_prev, U_low, U_norm, U_old, exp_val_old, mod_param, solver, sol_opt, 
                eta_div_cell_area, dx)
        t_curr += dt_use
        t_it += 1
        
        # Check if exocytosis time should be computed 
        if t_curr > t_hit:
            print("Exocytosis occuring")
            Cdc42T_tmp, tmp1, tmp2, tmp3, tmp4, P_tmp, tmp5, tmp6, tmp7 = U_high.split(True)
            if exocyte_data.use_cables == False:
                i_hit = exocyte_data.calc_site_hit(Cdc42T_tmp)
            else:
                i_hit = exocyte_data.calc_site_hit_cable()

            coord_hit_sphere = exocyte_data.get_coord_sphere_hit(i_hit)
            for j in range(exocyte_data.n_coord):
                exocyte_data.calc_coord_interpolate(j, coord_hit_sphere, exocyte=True, P_vec=P_tmp.vector().get_local(), mod_param=mod_param)
            
            exocyte_data.set_values_post_exocyte(U_high)
            U_high.vector()[:] = exocyte_data.U_new.vector()
            U_prev.assign(U_high)
            compute_time_exocyte = True
            vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
            vtkfile_P << (_P_high, t_curr)    
            vtkfile_GapS << (_GapS_high, t_curr)
            vtkfile_S << (_S_high, t_curr)
            vtkfile_Axl2 << (_Axl2_high, t_curr)
            vtkfile_Axl2_S << (_Axl2_S_high, t_curr)
            if exocyte_data != None and exocyte_data.plot_cables == True:
                exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

            # Correct the cytosolic diffusion using the practical total concentration 
            # obtained from solving for the steady state. 
            Cdc42T_high, Cdc42D_high, BemGEF42_high, BemGEFm_high, S_high, P_high, GapS_high, Axl2_high, Axl2_S_high = split(U_high)
            correct_cdc42 = mod_param.Cdc42_tot_pract - assemble((Cdc42T_high + Cdc42D_high + BemGEF42_high)*dx(1)) *eta_div_cell_area
            correct_BemGEF = mod_param.BemGEF_tot_pract - assemble((BemGEFm_high + BemGEF42_high)*dx(1)) *eta_div_cell_area
            correct_S = mod_param.S_tot_pract - assemble((S_high + P_high + Axl2_S_high)*dx) *eta_div_cell_area
            correct_GapSc = mod_param.GapS_tot_pract - assemble(GapS_high*dx) * eta_div_cell_area
            correct_Axl2 = mod_param.Axl2_tot_pract - assemble((Axl2_high + Axl2_S_high)*dx) *eta_div_cell_area
        
            exp_list[0].Cdc42c = correct_cdc42
            exp_list[1].BemGEFc = correct_BemGEF
            exp_list[2].Sc = correct_S
            exp_list[3].GapSc = correct_GapSc
            exp_list[4].Axl2 = correct_Axl2

        # Calculate new locations for cables if required 
        if exocyte_data != None and exocyte_data.use_cables == True:
            Cdc42T_tmp, tmp1, tmp2, tmp3, tmp4, P_tmp, tmp5, tmp6, tmp7 = U_high.split(True)
            exocyte_data.update_cable_pos_septin(t_curr, Cdc42T_tmp, P_tmp)

        # Set low accuracy solution at high accuracy solution        
        exp_val_old[0] = exp_list[0].Cdc42c
        exp_val_old[1] = exp_list[1].BemGEFc
        exp_val_old[2] = exp_list[2].Sc
        exp_val_old[3] = exp_list[3].GapSc
        exp_val_old[4] = exp_list[4].Axl2c
        U_prev.vector()[:] = U_high.vector()
        U_old.vector()[:] = U_high.vector()
        
        # Export solution to disk 
        if t_it % sol_opt.print_pwd_it == 0:
            vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
            vtkfile_P << (_P_high, t_curr)    
            vtkfile_GapS << (_GapS_high, t_curr)
            vtkfile_S << (_S_high, t_curr)
            vtkfile_Axl2 << (_Axl2_high, t_curr)
            vtkfile_Axl2_S << (_Axl2_S_high, t_curr)
            # If cable data is plotted 
            if exocyte_data != None and exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
                exocyte_data.plot_location_cables(U_high, _Cdc42T_high, file_loc, t_curr)

        if sol_opt.save_dist != False and t_it % sol_opt.save_dist == 0:
            data_dist = compute_save_dist_data(U_high, mesh, t_curr, file_loc)
            file_path = file_loc.dir_intermediate + "Distribution_data.csv"
            if not os.path.isfile(file_path):
                data_dist.to_csv(file_path)
            else:
                data_dist.to_csv(file_path, header=False, mode='a')
    
    # Write end coordinates
    vtkfile_Cdc42T << (_Cdc42T_high, t_curr)
    vtkfile_P << (_P_high, t_curr)    
    vtkfile_GapS << (_GapS_high, t_curr)
    vtkfile_S << (_S_high, t_curr)
    vtkfile_Axl2 << (_Axl2_high, t_curr)
    vtkfile_Axl2_S << (_Axl2_S_high, t_curr)

    if sol_opt.save_dist != False:
        data_dist = compute_save_dist_data(U_high, mesh, t_curr, file_loc)
        file_path = file_loc.dir_intermediate + "Distribution_data.csv"
        if not os.path.isfile(file_path):
            data_dist.to_csv(file_path)
        else:
            data_dist.to_csv(file_path, header=False, mode='a')
    
    # If saving run-result 
    if sol_opt.save_run:
        exp_list_vector = np.array([exp_list[0].Cdc42c, exp_list[1].BemGEFc])
        Common.save_sim_result(U_high, exp_list_vector, file_loc.dir_save_run)


def run_axl2_v1(file_loc, mod_param_list, sol_opt_list, times_run, n_threads, bud_scar_param=None, exocyte_data=None, couple_mod_param_and_exo=False):

    try:
        len_exo_list = len(exocyte_data)
    except:
        len_exo_list = 0

    print("len_exo_list = ", len_exo_list)

    input_list = []

    if len_exo_list == 0:
        for mod_param in mod_param_list:
            for sol_opt in sol_opt_list:
                file_loc_tmp = Common.FileLocations(file_loc.mesh_name, file_loc.save_name, mod_param.calc_save_tag(bud_scar_param), name_store=file_loc.name_store)
                input_list.append([copy.deepcopy(file_loc_tmp), copy.deepcopy(mod_param), copy.deepcopy(sol_opt),
                    copy.deepcopy(bud_scar_param), copy.deepcopy(exocyte_data)])
    elif couple_mod_param_and_exo == True:
        for i in range(len(mod_param_list)):
            for sol_opt in sol_opt_list:
                mod_param = mod_param_list[i]
                exo_data = exocyte_data[i]
                file_loc_tmp = Common.FileLocations(file_loc.mesh_name, file_loc.save_name, mod_param.calc_save_tag(bud_scar_param), name_store=file_loc.name_store)
                input_list.append([copy.deepcopy(file_loc_tmp), copy.deepcopy(mod_param), copy.deepcopy(sol_opt),
                        copy.deepcopy(bud_scar_param), copy.deepcopy(exo_data)])
    else:
        for mod_param in mod_param_list:
            for sol_opt in sol_opt_list:
                for exo_data in exocyte_data:
                    file_loc_tmp = Common.FileLocations(file_loc.mesh_name, file_loc.save_name, mod_param.calc_save_tag(bud_scar_param, exo_data),
                                                        name_store=file_loc.name_store)
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
