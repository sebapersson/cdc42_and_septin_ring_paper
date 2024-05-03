'''
This file contains functions required for finding the symmetry 
breaking steady states for u and v given a parameter set of 
non dimensional parameters. To actually obtain parameter 
values, run the function find_init_u_and_v
'''


import numpy as np
import warnings
from scipy.optimize import fsolve


# Function which solves for the steady state in (u, v) by finding
# where the null-clines for f(u, v) and q(u, v) intersect. Note,
# as q is a quadratic function the function will return
# two vectors, the first being for the positive null-cline
# of u, and the second being for the negative null-cline of
# q as a function of u. Note that the steady state is found
# numerically by solving for the part where the null-clines,
# as functions of u, are equal. 
# Args:
#     param, a model parameter object
#     u_start, a start-guess for the u-parameter 
# Returns:
#     ss_u, ss_v for both positive and negative
#         null-clines of q 
def solve_for_steady_state(param, u_start):
    
    # Nullcline for function f as function of u 
    null_cline_f = lambda u: (u / (param.c2 + u*u))
    
    # Nullcline for function q as function of u (note quadratic)
    A = lambda u: param.c_max - u
    B = lambda u: (param.V0_init - (param.a * u) )
    K = lambda u: ( 1/param.a * ((param.a*A(u) + B(u)) + (param.c_1 / param.c1 )))
    L = lambda u: ( (A(u) * B(u))/ param.a)
    null_cline_q_neg = lambda u: (0.5 * (K(u) - np.sqrt((K(u)*K(u)) - (4*L(u)))))
    null_cline_q_pos = lambda u: (0.5 * (K(u) + np.sqrt((K(u)*K(u)) - (4*L(u)))))
    
    # Solve for the steady states for both q-null-clines
    ss_u_pos = fsolve(lambda x : null_cline_f(x) - null_cline_q_pos(x), u_start, xtol=2.2204e-16)[0]
    ss_u_neg = fsolve(lambda x : null_cline_f(x) - null_cline_q_neg(x), u_start, xtol=2.2204e-16)[0]
    ss_v_pos, ss_v_neg = null_cline_f(ss_u_pos), null_cline_f(ss_u_neg)
    
    return [[ss_u_pos, ss_v_pos], [ss_u_neg, ss_v_neg]]


# Function that will check if a steady actually fulfills the steady state
# criteria that q = f = 0
# Args:
#     ss, a steady state (array on form [u, v])
#     param, the model parameters 
#     tol, tolerance when checking if u and v are sufficiently small 
# Returns:
#     true if is steady state, else false
def is_steady_state(ss, param, tol=1e-3):
    # f and q functions 
    f = lambda u,v: (param.c2 * v - u + u*u*v)
    q = lambda u,v: (( param.c1 * (param.c_max - (u+v)) * (param.V0_init - (param.a*(u+v))) ) -
                 (param.c1*v))
    # Check if steady state 
    u_ss, v_ss = ss
    if np.absolute(f(u_ss, v_ss)) < tol and np.absolute(q(u_ss, v_ss)):
        return True
    else:
        return False


# Function that will check if symmetry breaking conditions are fulfilled
# for a given steady state.
# Args:
#     ss, a steady state (array on form [u, v])
# Returns:
#     true or false depending on whether or not symmetry breaking
#     (str) classic, non_classic or no_symmetry depending on which conditions 
def check_if_symmetry_breaking(ss, param):
    
    # The partial derivatives
    f_u = lambda u,v: (2*u*v - 1)
    f_v = lambda u,v: (param.c2 + u*u )
    q_u = lambda u,v: ((-param.c1 * param.a) * (m - (u + v)))
    q_v = lambda u,v: (q_u(u, v) - param.c_1)
    q_V = lambda u,v:  (param.c1 * (param.c_max - (u + v )))
    # Parameters for checking conditions
    v_prime = -param.a
    m = param.V0_init / param.a
    
    # Check the conditions 
    u_ss, v_ss = ss
    
    # Negative trace of homogeneous system
    cond1 = (f_u(u_ss, v_ss) -  f_v(u_ss, v_ss) +  q_v(u_ss, v_ss) + (q_V(u_ss, v_ss) * v_prime)) < 0
    
    # Positive determinant of homogeneous system
    term1 = f_u(u_ss, v_ss) * ( q_v(u_ss, v_ss) + (q_V(u_ss, v_ss) * v_prime) );
    term2 = f_v(u_ss, v_ss) * ( q_u(u_ss, v_ss) + (q_V(u_ss, v_ss) * v_prime) );
    cond2 = (term1 - term2 ) > 0
    
    # Classic Turing 
    cond3 = (f_u(u_ss, v_ss)*q_v(u_ss, v_ss) - f_v(u_ss, v_ss)*q_u(u_ss, v_ss)) >= 0
    cond4 = (param.d*f_u(u_ss, v_ss) -f_v(u_ss, v_ss) + q_v(u_ss, v_ss)) > 0
    term3 = param.d*f_u(u_ss, v_ss) -f_v(u_ss, v_ss) + q_v(u_ss, v_ss)
    term4 = f_u(u_ss, v_ss)*q_v(u_ss, v_ss) - f_v(u_ss, v_ss)*q_u(u_ss, v_ss) 
    Q = term3*term3 - 4*param.d*term4
    cond5 = (term3*term3 - 4*param.d*term4) >= 0
    
    # Non classic Turing
    cond7 = term4 < 0
    cond8 = ((1/(2*param.d)) * (term3 + np.sqrt(Q))) > 0
    
    # Check if any symmetry breaking conditions are fulfilled
    if cond1 and cond2 and cond3 and cond4 and cond5:
        # Classic Turing
        return True, "classic"
    elif cond1 and cond2 and cond7 and cond8:
        # Non classic Turing
        return True, "non_classic"
    else:
        return False, "non_symmetry"


# Function that given a steady-state and parameter value checks that
# everything is positive
# Args:
#     param, the non-dimensional parameters for a model
#     ss, the initial values on the form [u, v]
# Returns:
#     true if everything is positive, else false 
def is_positive_ss(param, ss):
    u0, v0 = ss
    cond1 = u0 > 0
    cond2 = v0 > 0
    cond3 = (param.V0_init - param.a * (u0 + v0)) > 0
    
    return cond1 and cond2 and cond3 


# Function that given a parameter set for the model will find the
# initial values for u and v that fulfill symmetry breaking
# criteria. If no such values are found, the function returns
# false, which can be picked up by other functions
# Args:
#     param, the non-dimensional parameters for a model
# Returns:
#     ss, the initial values on the form [u, v]
#     info_ss, info of symmetry breaking condition
#     false for both argument if steady state not found
def find_init_u_and_v(param):
    
    # Limits where to search for parameters 
    u_min_limit = np.sqrt(param.c2)
    u_max_limit = np.min([param.c_max, param.V0_init / param.a])
    u_start_vec = np.linspace(u_min_limit, u_max_limit, num=1000)
    
    warnings.filterwarnings("ignore")
    for u_start in u_start_vec:
        ss_pos, ss_neg = solve_for_steady_state(param, u_start)
        sym_break_pos, info_ss_pos = check_if_symmetry_breaking(ss_pos, param)
        sym_break_neg, info_ss_neg = check_if_symmetry_breaking(ss_neg, param)
        
        # Return if Turing and symmetry breaking
        if is_steady_state(ss_pos, param) and sym_break_pos and is_positive_ss(param, ss_pos):
            warnings.filterwarnings("default")
            return ss_pos, info_ss_pos
        elif is_steady_state(ss_neg, param) and sym_break_neg and is_positive_ss(param, ss_neg):
            warnings.filterwarnings("default")
            return ss_neg, info_ss_neg
    
    # In case no initial value is found 
    return False, False 

