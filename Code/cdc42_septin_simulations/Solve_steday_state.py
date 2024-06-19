#!/usr/bin/env python

from dolfin import *
from tqdm import tqdm
import numpy as np
import pandas as pd
import scipy as sp
import itertools
import Common
import Config


def cost_function(x, mod_param):
    
    Cdc42T, Cdc42D, Cdc42c, BemGEF42, BemGEFm, BemGEFc = x
    
    k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, Dm_gdp = mod_param.return_param()
    
    # Cdc42T
    f1 = (k2a*BemGEFm + k3*BemGEF42) * Cdc42D - (k2b + k4a*BemGEFm + k7*BemGEFc)*Cdc42T + k4b*BemGEF42
    
    # Cdc42D
    f2 = k2b*Cdc42T - (k2a*BemGEFm + k3*BemGEF42)*Cdc42D - k5b*Cdc42D + k5a*Cdc42c
    
    # BemGEF42
    f3 = (k4a*BemGEFm + k7*BemGEFc)*Cdc42T - k4b*BemGEF42
    
    # BemGEF
    f4 = k1a*BemGEFc - k1b*BemGEFm + k4b*BemGEF42 - k4a*BemGEFm*Cdc42T
    
    ode1 = k1b*BemGEFm - k1a*BemGEFc - k7*BemGEFc*Cdc42T
    ode2 = k5b*Cdc42D - k5a*Cdc42c
    
    cost = f1**2 + f2**2 + f3**2 + f4**2 + ode1**2 + ode2**2
    
    return cost 


def solve_for_steady_state(mod_param):
    
    A_mat = np.matrix([[mod_param.eta*1.0, mod_param.eta*1.0, 1.0, mod_param.eta*1.0, 0.0, 0.0], 
                       [0.0, 0.0, 0.0, mod_param.eta*1.0, mod_param.eta*1.0, 1.0]])
    lb = np.array([mod_param.Cdc42_tot - mod_param.Cdc42_tot*0.05, mod_param.BemGEF_tot - mod_param.BemGEF_tot*0.02])
    ub = np.array([mod_param.Cdc42_tot + mod_param.Cdc42_tot*0.05, mod_param.BemGEF_tot + mod_param.BemGEF_tot*0.02])
    constraint = sp.optimize.LinearConstraint(A_mat, lb, ub)
    
    # Start guess based on Lew-lab paper 
    Cdc42T0 = 32.3633
    Cdc42D0 = 22.1990
    Cdc42c0 = 0.4008
    BemGEF420 = 1.4381
    BemGEFm0 = 0.0431
    BemGEFc0 = 0.0013
    x0 = np.array([Cdc42T0, Cdc42D0, Cdc42c0, BemGEF420, BemGEFm0, BemGEFc0])
    res = sp.optimize.minimize(cost_function, x0, args=(mod_param), constraints=(constraint), tol=1e-10)
    x_opt = res.x
    
    mod_param.ic_Cdc42T = x_opt[0]
    mod_param.ic_Cdc42D = x_opt[1]
    mod_param.ic_Cdc42I = x_opt[2]
    mod_param.ic_BemGEF42 = x_opt[3]
    mod_param.ic_BemGEFm = x_opt[4]
    mod_param.ic_BemGEFc = x_opt[5]

    eta = mod_param.eta
    mod_param.Cdc42_tot_pract = mod_param.ic_Cdc42T*eta + mod_param.ic_Cdc42D*eta + mod_param.ic_BemGEF42*eta + mod_param.ic_Cdc42I
    mod_param.BemGEF_tot_pract = mod_param.ic_BemGEFm*eta + mod_param.ic_BemGEFc + mod_param.ic_BemGEF42*eta

def solve_steady_state_axl2(mod_param, skip_opt=False):
    
    A_mat = np.matrix([[mod_param.eta*1.0, mod_param.eta*1.0, 1.0, mod_param.eta*1.0, 0.0, 0.0], 
                       [0.0, 0.0, 0.0, mod_param.eta*1.0, mod_param.eta*1.0, 1.0]])
                       
    lb = np.array([mod_param.Cdc42_tot - 0.05, 0.015])
    ub = np.array([mod_param.Cdc42_tot + 0.05, 0.019])
    constraint = sp.optimize.LinearConstraint(A_mat, lb, ub)
    
    # Start guess based on Lew-lab paper 
    Cdc42T0 = 32.3633
    Cdc42D0 = 22.1990
    Cdc42c0 = 0.4008
    BemGEF420 = 1.4381
    BemGEFm0 = 0.0431
    BemGEFc0 = 0.0013
    x0 = np.array([Cdc42T0, Cdc42D0, Cdc42c0, BemGEF420, BemGEFm0, BemGEFc0])
    if skip_opt == False:
        res = sp.optimize.minimize(cost_function_simple_septin, x0, args=(mod_param, "Axl2"), constraints=(constraint), tol=1e-10)
        x_opt = res.x
        mod_param.ic_Cdc42T = x_opt[0]
        mod_param.ic_Cdc42D = x_opt[1]
        mod_param.ic_Cdc42I = x_opt[2]
        mod_param.ic_BemGEF42 = x_opt[3]
        mod_param.ic_BemGEFm = x_opt[4]
        mod_param.ic_BemGEFc = x_opt[5]
    else:
        mod_param.ic_Cdc42T = Cdc42T0
        mod_param.ic_Cdc42D = Cdc42D0
        mod_param.ic_Cdc42I = Cdc42c0
        mod_param.ic_BemGEF42 = BemGEF420
        mod_param.ic_BemGEFm = BemGEFm0
        mod_param.ic_BemGEFc = BemGEFc0

    mod_param.ic_Sc =  mod_param.S_tot
    mod_param.ic_GapSc = mod_param.GapS_tot
    mod_param.ic_Axl2c = mod_param.Axl2_tot
    mod_param.ic_BemGEFm_S = 0.0
    mod_param.ic_BemGEF42_S = 0.0
    mod_param.ic_S = 0.0
    mod_param.ic_GapS = 0.0
    mod_param.ic_P = 0.0
    mod_param.ic_Axl2 = 0.0
    mod_param.ic_Axl2_S = 0.0
    
    eta = mod_param.eta
    mod_param.Cdc42_tot_pract = mod_param.ic_Cdc42T*eta + mod_param.ic_Cdc42D*eta + mod_param.ic_BemGEF42*eta + mod_param.ic_Cdc42I + mod_param.ic_BemGEF42_S*eta
    mod_param.BemGEF_tot_pract = mod_param.ic_BemGEFm*eta + mod_param.ic_BemGEFc + mod_param.ic_BemGEF42*eta + mod_param.ic_BemGEFm_S*eta + mod_param.ic_BemGEF42_S*eta
    mod_param.S_tot_pract = mod_param.ic_S*eta + mod_param.ic_Sc*eta + mod_param.ic_Sc + mod_param.ic_BemGEFm_S*eta + mod_param.ic_BemGEF42_S*eta
    mod_param.GapS_tot_pract = mod_param.ic_GapS*eta + mod_param.ic_GapSc
    mod_param.Axl2_tot_pract = mod_param.ic_Axl2c

def solve_steady_state_axl2_v4(mod_param, skip_opt=False):
    
    # Start guess based on Lew-lab paper 
    Cdc42T0 = 32.3633
    Cdc42D0 = 22.1990
    Cdc42c0 = 0.4008
    BemGEF420 = 1.4381
    BemGEFm0 = 0.0431
    BemGEFc0 = 0.0013
    mod_param.ic_Cdc42T = Cdc42T0
    mod_param.ic_Cdc42D = Cdc42D0
    mod_param.ic_Cdc42I = Cdc42c0
    mod_param.ic_BemGEF42 = BemGEF420
    mod_param.ic_BemGEFm = BemGEFm0
    mod_param.ic_BemGEFc = BemGEFc0

    mod_param.ic_Sc =  mod_param.S_tot
    mod_param.ic_GapSc = mod_param.GapS_tot
    mod_param.ic_Axl2c = mod_param.Axl2_tot
    mod_param.ic_BemGEFm_S = 0.0
    mod_param.ic_BemGEF42_S = 0.0
    mod_param.ic_S = 0.0
    mod_param.ic_GapS = 0.0
    mod_param.ic_P = 0.0
    mod_param.ic_Axl2 = 0.0
    mod_param.ic_Axl2_S = 0.0
    mod_param.ic_FOO = 0.0
    
    eta = mod_param.eta
    mod_param.Cdc42_tot_pract = mod_param.ic_Cdc42T*eta + mod_param.ic_Cdc42D*eta + mod_param.ic_BemGEF42*eta + mod_param.ic_Cdc42I + mod_param.ic_BemGEF42_S*eta
    mod_param.BemGEF_tot_pract = mod_param.ic_BemGEFm*eta + mod_param.ic_BemGEFc + mod_param.ic_BemGEF42*eta + mod_param.ic_BemGEFm_S*eta + mod_param.ic_BemGEF42_S*eta
    mod_param.S_tot_pract = mod_param.ic_S*eta + mod_param.ic_Sc*eta + mod_param.ic_Sc + mod_param.ic_BemGEFm_S*eta + mod_param.ic_BemGEF42_S*eta
    mod_param.GapS_tot_pract = mod_param.ic_GapS*eta + mod_param.ic_GapSc
    mod_param.Axl2_tot_pract = mod_param.ic_Axl2c    

def cost_function_nf(x, mod_param, double=False, print_i=False):

    Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star, Cdc42c, BemGEFc, BemGEFc_star = x
    if double == False:
        k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, k8max, k8h, k8n, k9max, k9h, k9n, Dm_gdp = mod_param.return_param()
    else:
        k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, k8max, k8h, k8n, k9max, k9h, k9n, k8max_2, k8h_2, k9max_2, k9h_2, Dm_gdp = mod_param.return_param()            
    k8 = k8max * (BemGEF42 + BemGEF42_star)**k8n / (k8h**k8n + (BemGEF42 + BemGEF42_star)**k8n)
    k9 = k9max * (BemGEFc_star**k9n) / (k9h**k9n + BemGEFc_star**k9n)
    k8 *= mod_param.strength_feedback
    k9 *= mod_param.strength_feedback

    k7_term = k7 * (BemGEFc_star  + BemGEFc)
    k4a_term = k4a * (BemGEFm + BemGEFm_star)
    k4b_term = k4b * (BemGEF42 + BemGEF42_star)
    k8_term = k8 * (BemGEF42 + BemGEF42_star) * BemGEF42
    k9_term = k9 * BemGEFc_star

    # Cdc42T
    f1 = (k2a*BemGEFm + k3*BemGEF42) * Cdc42D - (k2b + k4a_term + k7_term)*Cdc42T + k4b_term
        
    # Cdc42D
    f2 = k2b*Cdc42T - (k2a*BemGEFm + k3*BemGEF42)*Cdc42D - k5b*Cdc42D + k5a*Cdc42c
    
    # BemGEF42
    f3 = (k4a*BemGEFm + k7*BemGEFc)*Cdc42T - k4b*BemGEF42 - k8_term

    # BemGEF
    f4 = k1a*BemGEFc - k1b*BemGEFm + k4b*BemGEF42 - k4a*BemGEFm * Cdc42T

    # BemGEF42*
    f5 = (k4a * BemGEFm_star + k7*BemGEFc_star)*Cdc42T - k4b*BemGEF42_star + k8_term

    # BemGEF*
    f6 = k1a*BemGEFc_star - k1b*BemGEFm_star + k4b*BemGEF42_star - k4a*BemGEFm_star * Cdc42T

    # Cdc42c 
    ODE1 = k5b*Cdc42D - k5a*Cdc42c

    # BemGEFc
    ODE2 = eta * (k1b * BemGEFm - k1a * BemGEFc - k7*BemGEFc*Cdc42T) + k9_term

    # BemGEFc*
    ODE3 = eta * (k1b * BemGEFm_star - k1a * BemGEFc_star - k7*BemGEFc_star*Cdc42T) - k9_term

    cost = f1**2 + f2**2 + f3**2 + f4**2 + f5**2 + f6**2 + ODE1**2 + ODE2**2 + ODE3**2

    if print_i:
        print("f1 = {:.3e}, f2 = {:.3e}, f3 = {:.3e}, f4 = {:.3e}, f5 = {:.3e}, f6 = {:.3e}".format(f1, f2, f3, f4, f5, f6))
        print("ODE1 = {:.3e}, ODE2 = {:.3e}, ODE3 = {:.3e}".format(ODE1, ODE2, ODE3))

    return cost 
    

def solve_for_steady_state_nf(mod_param, double=False):
    
    # Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star, Cdc42c, BemGEFc, BemGEFc_star
    A_mat = np.matrix([[mod_param.eta*1.0, mod_param.eta*1.0, mod_param.eta*1.0, 0.0, mod_param.eta*1.0, 0.0, 1.0, 0.0, 0.0], 
                       [0.0, 0.0, mod_param.eta*1.0, mod_param.eta*1.0, mod_param.eta*1.0, mod_param.eta*1.0, 0.0, 1.0, 1.0]])
    lb = np.array([mod_param.Cdc42_tot - 0.05, mod_param.BemGEF_tot - 0.005])
    ub = np.array([mod_param.Cdc42_tot + 0.05, mod_param.BemGEF_tot + 0.005])
    constraint = sp.optimize.LinearConstraint(A_mat, lb, ub)
    
    # Start guess based on Lew-lab paper 
    Cdc42T0 = 32.3633
    Cdc42D0 = 22.1990
    BemGEF420 = 1.4381 / 2.0
    BemGEFm0 = 0.0431 / 2.0
    BemGEF42_star0 = 1.4381 / 2.0
    BemGEFm_star0 = 0.0431 / 2.0
    Cdc42c0 = 0.4008
    BemGEFc0 = 0.0013 / 2.0
    BemGEFc_star0 = 0.0013 / 2.0
    x0 = np.array([Cdc42T0, Cdc42D0, BemGEF420, BemGEFm0, BemGEF42_star0, BemGEFm_star0, Cdc42c0, BemGEFc0, BemGEFc_star0])
    bounds = tuple(itertools.repeat((0, None), 9))
    res = sp.optimize.minimize(cost_function_nf, x0, args=(mod_param, double), constraints=(constraint), tol=1e-7, 
        options={"maxiter" : 1000}, bounds=bounds)
    x_opt = res.x
    
    mod_param.ic_Cdc42T = x_opt[0]
    mod_param.ic_Cdc42D = x_opt[1]
    mod_param.ic_BemGEF42 = x_opt[2]
    mod_param.ic_BemGEFm = x_opt[3]
    mod_param.ic_BemGEF42_star = x_opt[4]
    mod_param.ic_BemGEFm_star = x_opt[5]
    mod_param.ic_Cdc42I = x_opt[6]
    mod_param.ic_BemGEFc = x_opt[7]
    mod_param.ic_BemGEFc_star = x_opt[8]

    eta = mod_param.eta
    mod_param.Cdc42_tot_pract = mod_param.ic_Cdc42T*eta + mod_param.ic_Cdc42D*eta + mod_param.ic_BemGEF42*eta + mod_param.ic_BemGEF42_star*eta + mod_param.ic_Cdc42I + mod_param.ic_Cdc42I
    mod_param.BemGEF_tot_pract = mod_param.ic_BemGEFm*eta + mod_param.ic_BemGEFc + mod_param.ic_BemGEF42*eta
    mod_param.BemGEF_star_tot_pract = mod_param.ic_BemGEFm_star*eta + mod_param.ic_BemGEFc_star + mod_param.ic_BemGEF42_star*eta
    
def solve_for_steady_state_nf_small(mod_param):
    

    # Cdc42T, Cdc42D, BemGEF42, BemGEFm, BemGEF42_star, BemGEFm_star, Cdc42c, BemGEFc, BemGEFc_star
    A_mat = np.matrix([[mod_param.eta*1.0, mod_param.eta*1.0, mod_param.eta*1.0, 0.0, mod_param.eta*1.0, 0.0, 1.0, 0.0, 0.0], 
                       [0.0, 0.0, mod_param.eta*1.0, mod_param.eta*1.0, mod_param.eta*1.0, mod_param.eta*1.0, 0.0, 1.0, 1.0]])
    lb = np.array([mod_param.Cdc42_tot - 0.05, mod_param.Cdc24_tot - 0.005])
    ub = np.array([mod_param.Cdc42_tot + 0.05, mod_param.Cdc24_tot + 0.005])
    constraint = sp.optimize.LinearConstraint(A_mat, lb, ub)
    
    # Start guess based on Lew-lab paper 
    Cdc42T0 = 32.3633
    Cdc42D0 = 22.1990
    BemGEF420 = 1.4381 / 2.0
    BemGEFm0 = 0.0431 / 2.0
    BemGEF42_star0 = 1.4381 / 2.0
    BemGEFm_star0 = 0.0431 / 2.0
    Cdc42c0 = 0.4008
    BemGEFc0 = 0.0013 / 2.0
    BemGEFc_star0 = 0.0013 / 2.0
    x0 = np.array([Cdc42T0, Cdc42D0, BemGEF420, BemGEFm0, BemGEF42_star0, BemGEFm_star0, Cdc42c0, BemGEFc0, BemGEFc_star0])
    bounds = tuple(itertools.repeat((0, None), 9))

    mod_param_new = Common.ModelParametersNf(L=2.5)
    
    res = sp.optimize.minimize(cost_function_septin_small, x0, args=(mod_param), constraints=(constraint), tol=1e-10, 
        options={"maxiter" : 1000}, bounds=bounds)
    x_opt = res.x

    mod_param.ic_Cdc42T = x_opt[0]
    mod_param.ic_Cdc42D = x_opt[1]
    mod_param.ic_BemGEF42 = x_opt[2]
    mod_param.ic_Cdc24 = x_opt[3]
    mod_param.ic_Cdc24m_star = x_opt[5]
    mod_param.ic_Cdc24c_star = x_opt[8]
    mod_param.ic_BemGEF42_star = x_opt[4]
    mod_param.ic_Cdc42I = x_opt[6]
    mod_param.ic_Cdc24c = x_opt[7]
    
def cost_function_okada(x, mod_param):
    
    (k1a, k1b, k2, k3, k4, k5, k6a, k6b, k7, k8a, k8b, k9a, k9b, k10, k11, k12a, k12b, k13a, k13b,
    k14, k15, k16, k17, k18, k19, k20, k23, k24a, k24b, k25, k26a, k26b, k27, eta, Dm_gdp) = mod_param.return_param()
    
    Cdc42T, Cdc42D, BemGEF42, Cdc24, GapS, F, Cdc42TF, FS, S, P, Cdc42I_exp, Cdc24c_exp, GapSc_exp, Fc_exp, Sc_exp, FSc_exp, Cdc42I, Ic_exp = x
    
    v1a = k1a * Cdc24c_exp
    v1b = k1b * Cdc24
    v2 = k2 * Cdc24 * Cdc42D
    v3 = k3 * Cdc42T
    v4 = k4 * GapS * Cdc42T
    v5 = k5 * GapSc_exp * Cdc42T
    v6a = k6a * Cdc42T * Cdc24
    v6b = k6b * BemGEF42
    v7 = k7 * Cdc42D * BemGEF42
    v8a = k8a * Cdc42I_exp
    v8b = k8b * Cdc42I
    v9a = k9a * Ic_exp * Cdc42D
    v9b = k9b * Cdc42I
    v10 = k10 * Cdc24c_exp * Cdc42T
    v11 = k11 * F
    v12a = k12a * Cdc42T * F
    v12b = k12b * Cdc42TF
    v13a = k13a * Cdc42T * Fc_exp
    v13b = k13b * Cdc42TF
    v14 = k14 * S
    v15 = k15 * S * S
    v16 = k16 * S
    v17 = k17 * P
    v18 = k18 * GapSc_exp * P
    v19 = k19 * GapS
    v20 = k20 * 0.9
    v23 = k23 * FS
    v24a = k24a * Fc_exp * Sc_exp
    v24b = k24b * FSc_exp
    v25 = k25 * Cdc42T * FSc_exp
    v26a = k26a * F * S
    v26b = k26b * FS
    v27 = k27 * Cdc42T * FS
    
    # Reaction terms
    # Cdc42T
    f1 = v2 - v3 - v4 - v5 - v6a + v6b + v7- v10 - v12a+ v12b-v13a + v13b - v20 - v25 - v27
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
    
    # Cdc24c
    ODE1 = eta * (v1b - v1a - v10) 
    # Cdc42c
    ODE2 = v8b - v8a
    # GapSc
    ODE3 = v19 - v18
    # Fc, surface integral integrates to the area times the concentrations 
    ODE4 = eta * (v11 - v13a + v13b) - v24a + v24b
    # Sc
    ODE5 = eta * (v14) - v24a + v24b
    # FSc
    ODE6 = eta * (v23 - v25) + v24a - v24b
    # Ic
    ODE7 = v9b - v9a
    
    cost = f1**2 + f2**2 + f3**2 + f4**2 + f5**2 + f6**2 + f7**2 + f8**2 + f9**2 + f10**2 + f11**2
    cost += ODE1**2 + ODE2**2 + ODE3**2 + ODE4**2 + ODE5**2 + ODE6**2 + ODE7**2
    
    return cost 

# Total concentrations for constraining the initial values 
def solve_steady_state_okada(mod_param):
    
    eta = mod_param.eta
    Cdc42_tot = mod_param.Cdc42_tot
    Cdc24_tot = mod_param.Cdc24_tot
    F_tot = mod_param.F_tot
    S_tot = mod_param.S_tot
    GapS_tot = mod_param.GapS_tot
    I_tot = mod_param.I_tot
    
    # Order of states: Cdc42T, Cdc42D, BemGEF42, Cdc24, GapS, F, Cdc42TF, FS, S, P, Cdc42I_exp, Cdc24c_exp, GapSc_exp, Fc_exp, Sc_exp, FSc_exp, Cdc42I, Ic
    A_mat = np.matrix([[0, 0, 1*eta, 1*eta, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], # For BemGEF-complex
                       [1*eta, 1*eta, 1*eta, 0, 0, 0, 1*eta, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1*eta, 0], # For Cdc42
                       [0, 0, 0, 0, 0, 1*eta, 1*eta, 1*eta, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], # For F 
                       [0, 0, 0, 0, 0, 0, 0, 1*eta, 1*eta, 1*eta, 0, 0, 0, 0, 1, 1, 0, 0], # For S 
                       [0, 0, 0, 0, 1*eta, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], # For G
                       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1*eta, 1]])
                       
    # Constrain 10 % within the reported value 
    lb = np.array([Cdc24_tot - Cdc24_tot*0.01, Cdc42_tot - Cdc42_tot*0.01, F_tot - F_tot*0.01, S_tot - S_tot*0.01, GapS_tot - GapS_tot*0.01, I_tot - I_tot*0.01])
    ub = np.array([Cdc24_tot + Cdc24_tot*0.01, Cdc42_tot + Cdc42_tot*0.04, F_tot + F_tot*0.01, S_tot + S_tot*0.01, GapS_tot + GapS_tot*0.01, I_tot + I_tot*0.01])
    constraint = sp.optimize.LinearConstraint(A_mat, lb, ub)
    
    # Setting start-guesses based on Lew-Lab model. For the Septin module I lump the pool into the cytosolic
    # component and derive the steady-state thereafter
    Cdc42T0 = 32.3633
    Cdc42D0 = 22.1990
    Cdc42c0 = 0.4008
    BemGEF420 = 1.4381
    Cdc240 = 0.0431
    Cdc24c0 = 0.0013
    start_guess = np.zeros(18)
    start_guess[0] = Cdc42T0
    start_guess[1] = Cdc42D0 / 6
    start_guess[2] = BemGEF420
    start_guess[3] = Cdc240
    start_guess[10] = Cdc42c0
    start_guess[11] = Cdc24c0
    start_guess[12] = GapS_tot
    start_guess[13] = F_tot
    start_guess[14] = S_tot
    start_guess[16] = Cdc42D0 + Cdc42T0
    start_guess[17] = I_tot
    
    # Only accept positive parameters 
    bounds = tuple(itertools.repeat((0, None), 18))
    res = sp.optimize.minimize(cost_function_okada, start_guess, args=(mod_param), constraints=(constraint), tol=1e-10,
                               options={"maxiter" : 1000}, bounds=bounds)
    x_opt = res.x
    
    print(res)
    mod_param.ic_Cdc42T = x_opt[0]
    mod_param.ic_Cdc42D = x_opt[1]
    mod_param.ic_BemGEF42 = x_opt[2]
    mod_param.ic_Cdc24 = x_opt[3]
    mod_param.ic_GapS = x_opt[4]
    mod_param.ic_F = x_opt[5]
    mod_param.ic_Cdc42TF = x_opt[6]
    mod_param.ic_FS = x_opt[7]
    mod_param.ic_S = x_opt[8]
    mod_param.ic_P = x_opt[9]
    mod_param.ic_Cdc42Ic = x_opt[10]
    mod_param.ic_Cdc24c = x_opt[11]
    mod_param.ic_GapSc = x_opt[12]
    mod_param.ic_Fc = x_opt[13]
    mod_param.ic_Sc = x_opt[14]
    mod_param.ic_FSc = x_opt[15]
    mod_param.ic_Cdc42I = x_opt[16]
    mod_param.ic_Ic = x_opt[17]

    state_list = ["Cdc42T", "Cdc42D", "BemGEF42", "Cdc24", "GapS", "F", "Cdc42TF", "FS", "S", 
        "P", "Cdc42Ic", "Cdc24C", "GapSc", "Fc", "Sc", "FSc", "Cdc42I", "Ic"]
    for i in range(len(state_list)):
        print("{} = {:.3e}".format(state_list[i], x_opt[i]))
    
    return mod_param
