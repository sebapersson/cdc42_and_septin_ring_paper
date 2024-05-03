from dolfin import *
from tqdm import tqdm
import numpy as np
import copy 
import multiprocessing as mp

import Common
import Config 

def init(l):
    Config.lock = l

class IC_septin(UserExpression):
    def __init__(self, *args, **kwargs):
        
        self.ic_Cdc42T = kwargs.pop('ic_Cdc42T')
        self.ic_Cdc42D = kwargs.pop('ic_Cdc42D')
        self.ic_Cdc42I = kwargs.pop('ic_Cdc42I')
        self.ic_BemGEF42 = kwargs.pop('ic_BemGEF42')
        self.ic_Cdc24 = kwargs.pop('ic_Cdc24')
        self.ic_GapS = kwargs.pop('ic_GapS')
        self.ic_F = kwargs.pop('ic_F')
        self.ic_Cdc42TF = kwargs.pop('ic_Cdc42TF')
        self.ic_FS = kwargs.pop('ic_FS')
        self.ic_S = kwargs.pop('ic_S')
        self.ic_P = kwargs.pop('ic_P')
        super(IC_septin, self).__init__(*args, **kwargs)

        print("self.ic_Cdc42T = {:.3e}".format(self.ic_Cdc42T))
        print("self.ic_Cdc42D = {:.3e}".format(self.ic_Cdc42D))
        print("self.ic_Cdc42I = {:.3e}".format(self.ic_Cdc42I))
        print("self.ic_BemGEF42 = {:.3e}".format(self.ic_BemGEF42))
        print("self.ic_Cdc24 = {:.3e}".format(self.ic_Cdc24))
        print("self.ic_GapS = {:.3e}".format(self.ic_GapS))
        print("self.ic_F = {:.3e}".format(self.ic_F))
        print("self.ic_Cdc42TF = {:.3e}".format(self.ic_Cdc42TF))
        print("self.ic_FS = {:.3e}".format(self.ic_FS))
        print("self.ic_S = {:.3e}".format(self.ic_S))
        print("self.ic_P = {:.3e}".format(self.ic_P))
        
    
    def eval(self, values, x):
        # Cdc42T, Cdc42D, BemGEF42, Cdc24, GapS, F, Cdc42TF, FS, S, P
        values[0] = self.ic_Cdc42T + np.abs(np.random.normal(scale=0.01))
        values[1] = self.ic_Cdc42D + np.abs(np.random.normal(scale=0.01))
        values[2] = self.ic_BemGEF42 + np.abs(np.random.normal(scale=0.01))
        values[3] = self.ic_Cdc24 + np.abs(np.random.normal(scale=0.01))
        values[4] = self.ic_GapS + np.abs(np.random.normal(scale=0.01))
        values[5] = self.ic_F + np.abs(np.random.normal(scale=0.01))
        values[6] = self.ic_Cdc42TF + np.abs(np.random.normal(scale=0.01))
        values[7] = self.ic_FS + np.abs(np.random.normal(scale=0.01))
        values[8] = self.ic_S + np.abs(np.random.normal(scale=0.01))
        values[9] = self.ic_P + np.abs(np.random.normal(scale=0.01))
        values[10] = self.ic_Cdc42I + np.abs(np.random.normal(scale=0.01))
    
    def value_shape(self):
        return(11,)        


def formulate_components_fem_2d(mod_param, trial_f, trial_f_prev, test_f, measures, mesh, exp_list):
    
    Cdc42T, Cdc42D, BemGEF42, Cdc24, GapS, F, Cdc42TF, FS, S, P, Cdc42I = trial_f
    Cdc42T_prev, Cdc42D_prev, BemGEF42_prev, Cdc24_prev, GapS_prev, F_prev, Cdc42TF_prev, FS_prev, S_prev, P_prev, Cdc42I_prev = trial_f_prev
    phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9, phi10, phi11 = test_f
    dx = measures
    Cdc42I_exp, Cdc24c_exp, GapSc_exp, Fc_exp, Sc_exp, FSc_exp, Ic_exp = exp_list
    
    (k1a, k1b, k2, k3, k4, k5, k6a, k6b, k7, k8a, k8b, k9a, k9b, k10, k11, k12a, k12b, k13a, k13b,
     k14, k15, k16, k17, k18, k19, k20, k23, k24a, k24b, k25, k26a, k26b, k27, Dm, Dms, Dmss, Dm_gdp) = mod_param.make_param_constants()
    
    Cdc42T_mass = (Cdc42T - Cdc42T_prev)*phi1 * dx(1)
    Cdc42D_mass = (Cdc42D - Cdc42D_prev)*phi2 * dx(1)
    BemGEF42_mass = (BemGEF42 - BemGEF42_prev)*phi3 * dx(1)
    Cdc24_mass = (Cdc24 - Cdc24_prev)*phi4 * dx(1)
    GapS_mass = (GapS - GapS_prev)*phi5 * dx(1)
    F_mass = (F - F_prev)*phi6 * dx(1)
    Cdc42TF_mass = (Cdc42TF - Cdc42TF_prev)*phi7 * dx(1)
    FS_mass = (FS - FS_prev)*phi8 * dx(1)
    S_mass = (S - S_prev)*phi9 * dx(1)
    P_mass = (P - P_prev)*phi10 * dx(1)
    Cdc42I_mass = (Cdc42I - Cdc42I_prev)*phi11 * dx(1)

    # Whether or not having a septin diffusion barrier 
    if mod_param.barrier == False:
        Cdc42T_stiff = Dm * dot(grad(Cdc42T), grad(phi1)) * dx(1)
        Cdc42D_stiff = Dm_gdp * dot(grad(Cdc42D), grad(phi2)) * dx(1)
        BemGEF42_stiff = Dm * dot(grad(BemGEF42), grad(phi3)) * dx(1)
        Cdc24_stiff = Dm * dot(grad(Cdc24), grad(phi4)) * dx(1)
        GapS_stiff = Dmss * dot(grad(GapS), grad(phi5)) * dx(1)
        F_stiff = Dm * dot(grad(F), grad(phi6)) * dx(1)
        Cdc42TF_stiff = Dm * dot(grad(Cdc42TF), grad(phi7)) * dx(1)
        FS_stiff = Dms * dot(grad(FS), grad(phi8)) * dx(1)
        S_stiff = Dms * dot(grad(S), grad(phi9)) * dx(1)
        P_stiff = Dmss * dot(grad(P), grad(phi10)) * dx(1)
        Cdc42I_stiff = Dm * dot(grad(Cdc42I), grad(phi11)) * dx(1)
    else:
        Kd = Constant(mod_param.Kd)
        f = Constant(50.0)
        # Dm
        Dm0 = Dm
        Dm1 = Constant(Dm0 / f)
        Dm = (Dm0 - Dm1) * exp(-P_prev / Kd) + Dm1
        # Dmss
        Dms0 = Dms
        Dms1 = Constant(Dms0 / f)
        Dms = (Dms0 - Dms1) * exp(-P_prev / Kd) + Dms1
        # Dmss
        Dmss0 = Dmss
        Dmss1 = Constant(Dmss0 / f)
        Dmss = (Dmss0 - Dmss1) * exp(-P_prev / Kd) + Dmss1
        # Dm_gpd
        Dm_gdp0 = Dm_gdp
        Dm_gdp1 = Constant(Dm_gdp0 / f)
        Dm_gdp = (Dm_gdp0 - Dm_gdp1) * exp(-P_prev / Kd) + Dm_gdp1
        
        Cdc42T_stiff = Dm * dot(grad(Cdc42T), grad(phi1)) * dx - (Dm1 - Dm)/Kd * dot(grad(Cdc42T_prev), grad(P_prev)) * phi1 * dx
        Cdc42D_stiff = Dm_gdp * dot(grad(Cdc42D), grad(phi2)) * dx - (Dm_gdp1 - Dm_gdp)/Kd * dot(grad(Cdc42D_prev), grad(P_prev)) * phi2 * dx
        BemGEF42_stiff = Dm * dot(grad(BemGEF42), grad(phi3)) * dx - (Dm1 - Dm)/Kd * dot(grad(BemGEF42_prev), grad(P_prev)) * phi3 * dx
        Cdc24_stiff = Dm * dot(grad(Cdc24), grad(phi4)) * dx - (Dm1 - Dm)/Kd * dot(grad(Cdc24_prev), grad(P_prev)) * phi4 * dx
        GapS_stiff = Dmss * dot(grad(GapS), grad(phi5)) * dx - (Dmss1 - Dmss)/Kd * dot(grad(GapS_prev), grad(P_prev)) * phi5 * dx
        F_stiff = Dm * dot(grad(F), grad(phi6)) * dx - (Dm1 - Dm)/Kd * dot(grad(F_prev), grad(P_prev)) * phi6 * dx
        Cdc42TF_stiff = Dm * dot(grad(Cdc42TF), grad(phi7)) * dx - (Dm1 - Dm)/Kd * dot(grad(Cdc42TF_prev), grad(P_prev)) * phi7 * dx
        FS_stiff = Dms * dot(grad(FS), grad(phi8)) * dx - (Dms1 - Dms)/Kd * dot(grad(FS_prev), grad(P_prev)) * phi8 * dx
        S_stiff = Dms * dot(grad(S), grad(phi9)) * dx - (Dms1 - Dms)/Kd * dot(grad(S_prev), grad(P_prev)) * phi9 * dx
        P_stiff = Dmss * dot(grad(P), grad(phi10)) * dx - (Dmss1 - Dmss)/Kd * dot(grad(P_prev), grad(P_prev)) * phi10 * dx
        Cdc42I_stiff = Dm * dot(grad(Cdc42I), grad(phi11)) * dx - (Dm1 - Dm)/Kd * dot(grad(Cdc42I_prev), grad(P_prev)) * phi11 * dx

    k16 = conditional(lt(P_prev, 25.0), 0.0, k16)

    if mod_param.crowding == True:
        X = 1 - (S_prev + P_prev + FS_prev) / mod_param.Kt
        k1a *= X 
        k9a *= X 
        k10 *= X 
        k13a *= X
        k18 *= X
        k25 *= X
        
    v1a = k1a * Cdc24c_exp
    v1b = k1b * Cdc24
    v2 = k2 * Cdc24_prev * Cdc42D_prev
    v3 = k3 * Cdc42T
    v4 = k4 * GapS_prev * Cdc42T_prev
    v5 = k5 * GapSc_exp * Cdc42T_prev
    v6a = k6a * Cdc42T_prev * Cdc24_prev
    v6b = k6b * BemGEF42
    v7 = k7 * Cdc42D_prev * BemGEF42_prev
    v8a = k8a * Cdc42I_exp
    v8b = k8b * Cdc42I
    v9a = k9a * Ic_exp * Cdc42D_prev
    v9b = k9b * Cdc42I
    v10 = k10 * Cdc24c_exp * Cdc42T_prev
    v11 = k11 * F
    v12a = k12a * Cdc42T_prev * F_prev
    v12b = k12b * Cdc42TF
    v13a = k13a * Cdc42T_prev * Fc_exp
    v13b = k13b * Cdc42TF
    v14 = k14 * S
    v15 = k15 * S_prev * S_prev
    v16 = k16 * S
    v17 = k17 * P
    v18 = k18 * GapSc_exp * P
    v19 = k19 * GapS
    v20 = k20 * 0.9
    v23 = k23 * FS
    v24a = k24a * Fc_exp * Sc_exp
    v24b = k24b * FSc_exp
    v25 = k25 * Cdc42T_prev * FSc_exp
    v26a = k26a * F_prev * S_prev
    v26b = k26b * FS
    v27 = k27 * Cdc42T_prev * FS_prev
    
    # Reaction terms
    # Cdc42T, Cdc42D, BemGEF42, Cdc24, GapS, F, Cdc42TF, FS, S, P
    # Cdc42T
    f1 = v2 - v3 - v4 - v5 - v6a + v6b + v7- v10 - v12a + v12b - v13a + v13b - v20 - v25 - v27
    # Cdc42D
    f2 = -v2 + v3 + v4 + v5 - v7 - v9a - v9b + v20
    # BemGEF42
    f3 = v6a - v6b + v10
    # Cdc24
    f4 = v1a - v1b - v6a + v6b
    # GapS
    f5 = v18 - v19
    # F
    f6 = -v11 - v12a + v12b - v26a + v26b
    # Cdc42TF
    f7 = v12a - v12b + v13a - v13b + v25 + v27
    # FS
    f8 = v26a - v26b - v27 - v23
    # S
    f9 = v17 - 2*v15 - v16 - v14 + v25 + v27 - v26a + v26b
    # P
    f10 = 2*v15 + v16 - v17
    # Cdc42I
    f11 = v8a - v8b + v9a - v9b
    
    reaction_form = (-f1*phi1*dx(1) - f2*phi2*dx(1) - f3*phi3*dx(1) - f4*phi4*dx(1) - f5*phi5*dx(1) -
                     f6*phi6*dx(1) - f7*phi7*dx(1) - f8*phi8*dx(1) - f9*phi9*dx(1) - f10*phi10*dx(1) - f11*phi11*dx(1))
    
    mass_form = Cdc42T_mass + Cdc42D_mass + BemGEF42_mass + Cdc24_mass + GapS_mass + F_mass + Cdc42TF_mass + FS_mass + S_mass + P_mass + Cdc42I_mass
    stiff_form = Cdc42T_stiff + Cdc42D_stiff + BemGEF42_stiff + Cdc24_stiff + GapS_stiff + F_stiff + Cdc42TF_stiff + FS_stiff + S_stiff + P_stiff + Cdc42I_stiff
    
    return mass_form, stiff_form, reaction_form


def compute_fem_system_2d(mass_form, stiff_form, reaction_form):
    
    # Forms for vectors on RHS
    mass_form_rhs = rhs(mass_form)
    stiffness_reaction_form_rhs = rhs(stiff_form + reaction_form)
    
    mass_form_lhs = lhs(mass_form)
    stiffness_reaction_form_lhs = lhs(stiff_form + reaction_form)
    
    return mass_form_lhs, stiffness_reaction_form_lhs, mass_form_rhs, stiffness_reaction_form_rhs


def calc_rates_ode(U, exp_list, mod_param):
    
    Cdc42c, Cdc24c, GapSc, Fc, Sc, FSc, Ic = exp_list
    Cdc42T, Cdc42D, BemGEF42, Cdc24, GapS, F, Cdc42TF, FS, S, P, Cdc42I = split(U)
    
    (k1a, k1b, k2, k3, k4, k5, k6a, k6b, k7, k8a, k8b, k9a, k9b, k10, k11, k12a, k12b, k13a, k13b,
     k14, k15, k16, k17, k18, k19, k20, k23, k24a, k24b, k25, k26a, k26b, k27, Dm, Dms, Dmss, Dm_gdp) = mod_param.make_param_constants()

    if mod_param.crowding == True:
        X = 1 - (S + P + FS) / mod_param.Kt
        k1a *= X 
        k9a *= X 
        k10 *= X 
        k13a *= X
        k18 *= X
        k25 *= X
    
    v8a = k8a * Cdc42c
    v8b = k8b * Cdc42I
    v9a = k9a * Ic * Cdc42D
    v9b = k9b * Cdc42I

    v1a = k1a * Cdc24c
    v10 = k10 * Cdc24c * Cdc42T
    v13a = k13a * Cdc42T * Fc
    v18 = k18 * GapSc * P
    v25 = k25 * Cdc42T * FSc

    v1b = k1b * Cdc24
    v11 = k11 * F
    v13b = k13b * Cdc42TF
    v14 = k14 * S
    v19 = k19 * GapS
    v23 = k23 * FS
      
    return v1a, v1b, v8a, v8b, v9a, v9b, v10, v11, v13a, v13b, v14, v18, v19, v23, v25


def solve_ode_system(U, mod_param, exp_list, dt, eta_div_cell_area, dx):
    
    # Solving the ODE
    x_old = np.zeros(7)
    x_old[0] = exp_list[0].Cdc42c
    x_old[1] = exp_list[1].Cdc24c
    x_old[2] = exp_list[2].GapSc
    x_old[3] = exp_list[3].Fc
    x_old[4] = exp_list[4].Sc
    x_old[5] = exp_list[5].FSc
    x_old[6] = exp_list[6].Ic
    
    v1a, v1b, v8a, v8b, v9a, v9b, v10, v11, v13a, v13b, v14, v18, v19, v23, v25 = calc_rates_ode(U, exp_list, mod_param)
    v24a = mod_param.k24a * exp_list[3].Fc * exp_list[4].Sc
    v24b = mod_param.k24b * exp_list[5].FSc
    ode_list_old = np.zeros(7)
    # Cdc42c
    eta_div_cell_area = Constant(eta_div_cell_area)
    ode_form0 = (eta_div_cell_area * (v8b - v8a))*dx(1)
    ode_list_old[0] = assemble(ode_form0)
    # Cdc24c
    ode_form1 = (eta_div_cell_area * (v1b - v1a - v10))*dx(1)
    ode_list_old[1] = assemble(ode_form1)
    # GapSc
    ode_form2 = (eta_div_cell_area * (v19 - v18))*dx(1)
    ode_list_old[2] = assemble(ode_form2)
    # Fc, surface integral integrates to the area times the concentrations 
    ode_form3 = (eta_div_cell_area * (v11 - v13a + v13b))*dx(1)
    ode_list_old[3] = assemble(ode_form3) - v24a + v24b
    # Sc
    ode_form4 = (eta_div_cell_area * (v14))*dx(1)
    ode_list_old[4] = assemble(ode_form4) - v24a + v24b
    # FSc
    ode_form5 = (eta_div_cell_area * (v23 - v25))*dx(1)
    ode_list_old[5] = assemble(ode_form5) + v24a - v24b
    # Ic 
    ode_form6 = (eta_div_cell_area * (v9b - v9a)) * dx(1)
    ode_list_old[6] = assemble(ode_form6)

    # k1 RK2 method 
    k1 = dt * ode_list_old

    # k2 RK2
    exp_list[0].Cdc42c += k1[0] * 0.5
    exp_list[1].Cdc24c += k1[1] * 0.5
    exp_list[2].GapSc += k1[2] * 0.5
    exp_list[3].Fc += k1[3] * 0.5
    exp_list[4].Sc += k1[4] * 0.5
    exp_list[5].FSc += k1[5] * 0.5
    exp_list[6].Ic += k1[6] * 0.5
    v24a = mod_param.k24a * exp_list[3].Fc * exp_list[4].Sc
    v24b = mod_param.k24b * exp_list[5].FSc
    ode_list_old[0] = assemble(ode_form0)
    ode_list_old[1] = assemble(ode_form1)
    ode_list_old[2] = assemble(ode_form2)
    ode_list_old[3] = assemble(ode_form3) - v24a + v24b
    ode_list_old[4] = assemble(ode_form4) - v24a + v24b
    ode_list_old[5] = assemble(ode_form5) + v24a - v24b
    ode_list_old[6] = assemble(ode_form6)
    k2 = dt * ode_list_old
    
    x_new = x_old + 0.5*(k1 + k2)

    exp_list[0].Cdc42c = x_new[0]
    exp_list[1].Cdc24c = x_new[1]
    exp_list[2].GapSc = x_new[2]
    exp_list[3].Fc = x_new[3]
    exp_list[4].Sc = x_new[4]
    exp_list[5].FSc = x_new[5]
    exp_list[6].Ic  = x_new[6]


def update_exp_list(old_list, new_list):
    
    old_list[0].Cdc42c = new_list[0].Cdc42c
    old_list[1].Cdc24c = new_list[1].Cdc24c
    old_list[2].GapSc = new_list[2].GapSc
    old_list[3].Fc = new_list[3].Fc
    old_list[4].Sc = new_list[4].Sc
    old_list[5].FSc = new_list[5].FSc
    old_list[6].Ic = new_list[6].Ic


def solve_okada(file_loc, mod_param, sol_opt, exocyte_data=None):

    if sol_opt.seed != False:
        np.random.seed(sol_opt.seed)
    
    Config.lock.acquire()
    file_loc.process_msh_file_to_xdmf_2d()
    file_loc.calc_pole_index()
    mesh, subdomains = file_loc.read_mesh_2d()
    Config.lock.release()
    
    dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
    
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement([P1, P1, P1, P1, P1, P1, P1, P1, P1, P1, P1])
    H = FunctionSpace(mesh, element)

    if exocyte_data != None:
        exocyte_data.calc_coords(H)
        compute_time_exocyte = True
        t_hit = 0.0
    else:
        compute_time_exocyte = False
        t_hit = np.inf
    
    u0_exp = IC_septin(ic_Cdc42T = mod_param.ic_Cdc42T,
                       ic_Cdc42D = mod_param.ic_Cdc42D,
                       ic_BemGEF42 = mod_param.ic_BemGEF42,
                       ic_Cdc24 = mod_param.ic_Cdc24,
                       ic_GapS = mod_param.ic_GapS,
                       ic_F = mod_param.ic_F,
                       ic_Cdc42TF = mod_param.ic_Cdc42TF,
                       ic_FS = mod_param.ic_FS,
                       ic_S = mod_param.ic_S,
                       ic_P = mod_param.ic_P, 
                       ic_Cdc42I = mod_param.ic_Cdc42I)
    
    # Low resolution solution 
    test_f_low = TestFunctions(H)
    trial_f_low = TrialFunctions(H)
    U_prev_low = Function(H) 
    U_prev_low.interpolate(u0_exp)
    Cdc42I_exp_low = Expression('Cdc42c', Cdc42c = mod_param.ic_Cdc42Ic, element=P1)
    Cdc24c_exp_low = Expression('Cdc24c', Cdc24c = mod_param.ic_Cdc24c, element=P1)
    GapSc_exp_low = Expression('GapSc', GapSc = mod_param.ic_GapSc, element=P1)
    Fc_exp_low = Expression('Fc', Fc = mod_param.ic_Fc, element=P1)
    Sc_exp_low = Expression('Sc', Sc = mod_param.ic_Sc, element=P1)
    FSc_exp_low = Expression('FSc', FSc = mod_param.ic_FSc, element=P1)
    Ic_exp_low = Expression('Ic', Ic = mod_param.ic_Ic, element=P1)
    exp_list_low = [Cdc42I_exp_low, Cdc24c_exp_low, GapSc_exp_low, Fc_exp_low, Sc_exp_low, FSc_exp_low, Ic_exp_low]

    mass_low, stiff_low, reaction_low = formulate_components_fem_2d(mod_param, trial_f_low, U_prev_low, test_f_low, dx, mesh, exp_list_low)
    mass_lhs_low, stiffness_reaction_lhs_low, mass_rhs_low, stiffness_reaction_rhs_low = compute_fem_system_2d(mass_low, stiff_low, reaction_low)
    U_low = Function(H)
    
    # High resolution solution 
    test_f_high = TestFunctions(H)
    trial_f_high = TrialFunctions(H)
    U_prev_high = Function(H) 
    U_prev_high.vector()[:] = U_prev_low.vector()
    Cdc42I_exp_high = Expression('Cdc42c', Cdc42c = mod_param.ic_Cdc42Ic, element=P1)
    Cdc24c_exp_high = Expression('Cdc24c', Cdc24c = mod_param.ic_Cdc24c, element=P1)
    GapSc_exp_high = Expression('GapSc', GapSc = mod_param.ic_GapSc, element=P1)
    Fc_exp_high = Expression('Fc', Fc = mod_param.ic_Fc, element=P1)
    Sc_exp_high = Expression('Sc', Sc = mod_param.ic_Sc, element=P1)
    FSc_exp_high = Expression('FSc', FSc = mod_param.ic_FSc, element=P1)
    Ic_exp_high = Expression('Ic', Ic = mod_param.ic_Ic, element=P1)
    exp_list_high = [Cdc42I_exp_high, Cdc24c_exp_high, GapSc_exp_high, Fc_exp_high, Sc_exp_high, FSc_exp_high, Ic_exp_high]
    update_exp_list(exp_list_high, exp_list_low)
    
    mass_high, stiff_high, reaction_high = formulate_components_fem_2d(mod_param, trial_f_high, U_prev_high, test_f_high, dx, mesh, exp_list_high)
    mass_lhs_high, stiffness_reaction_lhs_high, mass_rhs_high, stiffness_reaction_rhs_high = compute_fem_system_2d(mass_high, stiff_high, reaction_high)
    
    U_high = Function(H)
    U_high.assign(U_prev_high)
    _Cdc42T, _Cdc42D, _BemGEF42, _Cdc24, _GapS, _F, _Cdc42TF, _FS, _S, _P, _Cdc42I = U_high.split()
    
    eta_div_cell_area = mod_param.eta / (4*np.pi*(mod_param.r / mod_param.L)**2)
    dt_low = 1e-3
    dt_high = dt_low / 2.0
    t_curr = 0.0
    
    vtkfile_Cdc42T = File(file_loc.dir_pvd + "Cdc42T" + ".pvd")
    vtkfile_S = File(file_loc.dir_pvd + "S" + ".pvd")
    vtkfile_P = File(file_loc.dir_pvd + "P" + ".pvd")
    vtkfile_GapS = File(file_loc.dir_pvd + "GapS" + ".pvd")
    vtkfile_Cdc42T << (_Cdc42T, t_curr)
    vtkfile_S << (_S, t_curr)
    vtkfile_P << (_P, t_curr)
    vtkfile_GapS << (_GapS, t_curr)
    
    # For adaptive time-stepping 
    U_norm = Function(H)
    U_old = Function(H)
    U_old.vector()[:] = U_prev_high.vector()    
    solver = KrylovSolver("gmres", "ilu")
    #solver.parameters["nonzero_initial_guess"] = True

    stiffness_reaction_lhs_high = assemble(stiffness_reaction_lhs_high)
    mass_lhs_high = assemble(mass_lhs_high)    
    stiffness_reaction_lhs_low = assemble(stiffness_reaction_lhs_low)
    mass_lhs_low = assemble(mass_lhs_low)
    
    # For more efficient matrix calculations 
    A_low = mass_lhs_low + dt_low * stiffness_reaction_lhs_low
    A_low_mass = 1.0 * mass_lhs_low + 0.0 * stiffness_reaction_lhs_low
    A_low_stiff = 0.0 * mass_lhs_low + 1.0 * stiffness_reaction_lhs_low
    A_high = mass_lhs_high + dt_high * stiffness_reaction_lhs_high
    A_high_mass = 1.0 * mass_lhs_high + 0.0 * stiffness_reaction_lhs_high
    A_high_stiff = 0.0 * mass_lhs_high + 1.0 * stiffness_reaction_lhs_high

    U_high.set_allow_extrapolation(True)

    # Add actin cables if used 
    if exocyte_data != None and exocyte_data.use_cables == True:
        for i in range(exocyte_data.n_cables):
            Cdc42T_tmp, Cdc42D_tmp, BemGEF42_tmp, Cdc24_tmp, GapS_tmp, F_tmp, Cdc42TF_tmp, FS_tmp, S_tmp, P_tmp, Cdc42I_tmp = U_high.split(True)
            exocyte_data.calc_cable_attach_septin(Cdc42T_tmp, P_tmp, i, n_cables_use=i)

    # If location of cables are plotted 
    if exocyte_data != None and exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
        exocyte_data.plot_location_cables(U_high, _Cdc42T, file_loc, t_curr)

    Cdc42T_high, Cdc42D_high, BemGEF42_high, Cdc24_high, GapS_high, F_high, Cdc42TF_high, FS_high, S_high, P_high, Cdc42I_high = U_high.split(True)
    
    t_it = 1
    dt_use = 0.0
    delta_U = Function(H)
    for i in tqdm(range(sol_opt.term_max_it)):

        if exocyte_data != None and compute_time_exocyte == True:
            t_hit = exocyte_data.calc_time_hit()
            compute_time_exocyte = False
            t_hit = t_curr + t_hit
            print("Time-hit = {}".format(t_hit))

        while True: 
            dt_use = dt_low
            
            A_low.zero()
            A_low.axpy(1.0, A_low_mass, True)
            A_low.axpy(dt_low, A_low_stiff, True)
            
            b_low = assemble(mass_rhs_low + Constant(dt_low) * stiffness_reaction_rhs_low)
            solver.solve(A_low,  U_low.vector(), b_low)
            delta_U.vector()[:] = U_low.vector() - U_prev_low.vector()
            
            A_high.zero()
            A_high.axpy(1.0, A_high_mass, True)
            A_high.axpy(dt_high, A_high_stiff, True)

            # High accuracy update t -> t + dt_high
            b_high = assemble(mass_rhs_high + Constant(dt_high)*stiffness_reaction_rhs_high)
            solver.solve(A_high,  U_high.vector(), b_high)
            U_prev_high.assign(U_high)
            solve_ode_system(U_high, mod_param, exp_list_high, dt_high, eta_div_cell_area, dx)

            # # High accuracy update t + dt_high -> t + dt_high*2
            b_high = assemble(mass_rhs_high + Constant(dt_high)*stiffness_reaction_rhs_high)
            solver.solve(A_high,  U_high.vector(), b_high)
            U_prev_high.assign(U_high)
            
            U_norm.vector()[:] = U_low.vector() - U_high.vector()
            dt_low, dt_high, progress = Common.update_step_length(U_norm, dt_low, sol_opt.adaptive_rel_tol, sol_opt.adaptive_rho)
            if progress == False:
                update_exp_list(exp_list_high, exp_list_low)
                U_prev_high.vector()[:] = U_old.vector()
            else:
                # Only update ODE-solution if time-step is accepted
                solve_ode_system(U_high, mod_param, exp_list_high, dt_high, eta_div_cell_area, dx)
                break
        
        t_curr += dt_use

        # Check if exocytosis time should be computed 
        if t_curr > t_hit:
            print("Exocytosis occuring")
            # Add mass correction for exocytosis events 
            Cdc42_pract = exp_list_high[0].Cdc42c + assemble((Cdc42T_high + Cdc42D_high + Cdc42I_high + BemGEF42_high + Cdc42TF_high)*dx) * eta_div_cell_area
            Cdc24_pract = exp_list_high[1].Cdc24c + assemble((Cdc24_high + BemGEF42_high)*dx) * eta_div_cell_area
            GapS_pract = exp_list_high[2].GapSc + assemble(GapS_high * dx) * eta_div_cell_area
            F_pract = exp_list_high[3].Fc + assemble((F_high + Cdc42TF_high)*dx) * eta_div_cell_area
            S_pract = exp_list_high[4].Sc + assemble((FS_high + S_high + P_high)*dx) * eta_div_cell_area
            FS_pract = exp_list_high[5].FSc + assemble(FS_high*dx) * eta_div_cell_area
            Ic_pract = exp_list_high[6].Ic + assemble(Cdc42I_high*dx) * eta_div_cell_area

            Cdc42T_tmp, Cdc42D_tmp, BemGEF42_tmp, Cdc24_tmp, GapS_tmp, F_tmp, Cdc42TF_tmp, FS_tmp, S_tmp, P_tmp = U_high.split(True)
            if exocyte_data.use_cables == False:
                i_hit = exocyte_data.calc_site_hit_septin(Cdc42T_tmp, P_tmp)
            else:
                i_hit = exocyte_data.calc_site_hit_cable()

            coord_hit_sphere = exocyte_data.get_coord_sphere_hit(i_hit)
            
            for j in range(exocyte_data.n_coord):
                exocyte_data.calc_coord_interpolate(j, coord_hit_sphere)
            
            exocyte_data.set_values_post_exocyte(U_high)
            U_high.vector()[:] = exocyte_data.U_new.vector()
            U_prev_high.assign(U_high)
            compute_time_exocyte = True

            exp_list_high[0].Cdc42c = Cdc42_pract - assemble((Cdc42T_high + Cdc42I_high + Cdc42D_high + BemGEF42_high + Cdc42TF_high)*dx) * eta_div_cell_area
            exp_list_high[1].Cdc24c = Cdc24_pract - assemble((Cdc24_high + BemGEF42_high)*dx) * eta_div_cell_area
            exp_list_high[2].GapSc = GapS_pract - assemble(GapS_high * dx) * eta_div_cell_area
            exp_list_high[3].Fc = F_pract - assemble((F_high + Cdc42TF_high)*dx) * eta_div_cell_area
            exp_list_high[4].Sc = S_pract - assemble((FS_high + S_high + P_high)*dx) * eta_div_cell_area
            exp_list_high[5].FSc = FS_pract - assemble(FS_high*dx) * eta_div_cell_area
            exp_list_high[6].Ic = Ic_pract - assemble(Cdc42I_high*dx) * eta_div_cell_area

            vtkfile_Cdc42T << (_Cdc42T, t_curr)
            vtkfile_S << (_S, t_curr)
            vtkfile_P << (_P, t_curr)
            vtkfile_GapS << (_GapS, t_curr)

            if exocyte_data != None and exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
                exocyte_data.plot_location_cables(U_high, _Cdc42T, file_loc, t_curr)

        # Calculate new locations for cables if required 
        if exocyte_data != None and exocyte_data.use_cables == True:
            Cdc42T_tmp, Cdc42D_tmp, BemGEF42_tmp, Cdc24_tmp, GapS_tmp, F_tmp, Cdc42TF_tmp, FS_tmp, S_tmp, P_tmp, Cdc42I_tmp = U_high.split(True)
            exocyte_data.update_cable_pos_septin(t_curr, Cdc42T_tmp, P_tmp)

        if t_it % sol_opt.print_pwd_it == 0:
            vtkfile_Cdc42T << (_Cdc42T, t_curr)
            vtkfile_S << (_S, t_curr)
            vtkfile_P << (_P, t_curr)
            vtkfile_GapS << (_GapS, t_curr)
            if exocyte_data != None and exocyte_data.plot_cables == True and exocyte_data.use_cables == True:
                exocyte_data.plot_location_cables(U_high, _Cdc42T, file_loc, t_curr)
        t_it += 1

        U_prev_low.vector()[:] = U_high.vector()
        U_old.vector()[:] = U_high.vector()
        update_exp_list(exp_list_low, exp_list_high)

    # Write end coordinates
    vtkfile_Cdc42T << (_Cdc42T, t_curr)
    vtkfile_S << (_S, t_curr)
    vtkfile_P << (_P, t_curr)

    
def run_okada_model(file_loc, mod_param_list, sol_opt_list, times_run, n_threads, exocyte_data=None):
    
    input_list = []
    for mod_param in mod_param_list:
        for sol_opt in sol_opt_list:
            file_loc_tmp = Common.FileLocations(file_loc.mesh_name, file_loc.save_name, mod_param.calc_save_tag(), name_store=file_loc.name_store)
            input_list.append([copy.deepcopy(file_loc_tmp), copy.deepcopy(mod_param), copy.deepcopy(sol_opt), copy.deepcopy(exocyte_data)])
    
    len_input_list = len(input_list)
    for i in range(times_run-1):
        for j in range(len_input_list):
            input_list.append(copy.deepcopy(input_list[j]))
    
    print("Length of input argument list = {}".format(len(input_list)))
    
    print("Starting simulations for 2d system")
    l = mp.Lock()
    pool = mp.Pool(initializer=init, initargs=(l,), processes=n_threads)
    result = pool.starmap(solve_okada, input_list)
    pool.close()
    pool.join()
    print("Done with simulations for 2d system")