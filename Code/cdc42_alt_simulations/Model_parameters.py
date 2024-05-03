"""
This file contains functions and classes related to the model-rate parameters 
and initial conditions. 

Overall this file contains:
  Classes that stores parameters (with or without dim) and class that stores bud-scar parameters 
  Functions for non-dimensionalization, calculating parameter tag for file-system, converting parameters 
    to Fenics-constants, and creating parameters-lists when running simulations. 
"""


from dolfin import *
import numpy as np
import itertools
import copy

import Find_steady_state


# Class holding the dimensionless parameters for a model. Note, standard values
# are assumed, which can be changed by the user. Furthermore, note
# that these are non-dimensional parameters. 
# Args:
#     a, proportion of the surface area over the volume
#     D, cytosolic diffusion component
#     d, relative diffusion coefficient
#     c1, relative influx of cdc42 from the cytsol
#     c_1, relative dissociation of inactive cdc42 from membrane
#     c2, relative activation rate of cdc42
#     gamma, reaction strength parameters 
#     c_max, maximum concentration in the membrane
#     V0_init, initial value for calculating V0
#     u0, initial value for u
#     v0, initial value for v
#     sym_break_cond, the symmetry breaking condition for a model
#     find_mesh, if the mesh is to be found for the parameter combination
#     scale_sigma, whether or not the sigma value should be scaled (if
#       false standard sigma is used)
class ModelParametersNonDim:
    def __init__(self, a=3, D=10000, d=10, c1=0.05, c_1=0.04, c2=0.45,
                 gamma=20, c_max=3, V0_init=6, u0=1.2263, v0=0.6276,
                 sym_break_cond='turing', scale_sigma=False, time_conversion_factor=1.0):
        self.a = a
        self.D = D
        self.d = d
        self.c1 = c1
        self.c_1 = c_1
        self.c2 = c2
        self.gamma = gamma
        self.c_max = c_max
        self.u0 = u0
        self.v0 = v0
        # Note, V0 is calculated
        self.V0_init = V0_init
        self.V0 = V0_init - (a*(u0 + v0))
        self.sym_break_cond = sym_break_cond
        self.scale_sigma = scale_sigma
        self.time_conversion_factor = time_conversion_factor
    
    # Method that prints the parameters values
    def print_param(self):
        print("a={:.3e}, D={:.3e}, d={:.3e}, u0={:.3e}, v0={:.3e}, V0_init={:.3e}, V0={:.3e}".
              format(self.a, self.D, self.d, self.u0, self.v0, self.V0_init, self.V0))
        print("c1={:.3e}, c_1={:.3e}, c2={:.3e}, c_max={:.3e}, gamma={:.3e}".
              format(self.c1, self.c_1, self.c2, self.c_max, self.gamma))
        print("Symmetry breaking condition = {}".format(self.sym_break_cond))


# Class holding the parameters, with dimensions, for the model.
# Note that standard values exist. The scale sigma is an indicator
# whether or not the sigma value should be scaled according to the
# concentration of a certain specie. 
class ModelParametersWithDim:
    def __init__(self, k1=0.28, k_1=0.133, k2=0.001368, k_2=0.01733*1.6, k3=1,
                 k_max=2.5400, DA=0.011, DI=0.023, DG=10, R=3.95, G0=1,
                 d=None, scale_sigma=False):
        self.k1 = k1
        self.k_1 = k_1
        self.k2 = k2
        self.k_2 = k_2
        self.k3 = k3
        self.k_max = k_max
        self.DA = DA
        self.DG = DG
        self.R = R
        self.G0 = G0
        self.scale_sigma = scale_sigma
        
        # If value for d (DI / DA) is provided set DI accordingly
        if d != None:
            self.d = d
            self.DI = self.DA * d
        else:
            self.DI = DI
            self.d = self.DI/self.DA
        
    # Method that prints the parameter values
    def print_param(self):
        print("k1={:.3e}, k_1={:.3e}, k2={:.3e}, k_2={:.3e}, k3={:.3e}, k_max={:.3e}".
              format(self.k1, self.k_1, self.k2, self.k_2, self.k3, self.k_max))
        print("DA={:.3e}, DI={:.3e}, DG={:.3e}, R={:.3e}, G0={:.3e}, d={:.3e}".
              format(self.DA, self.DI, self.DG, self.R, self.G0, self.DI/self.DA))


# Class holding the initial values for a model
# Args:
#     ic_u, initial value of u (membrane component)
#     ic_v, initial value of v (membrane component)
#     ic_V, initial value of V (cytosol component)
#     sigma, noise-level to perturb steady state with 
#     rc, radius for cytosol nodes 
class ic(UserExpression):
    def __init__(self, *args, **kwargs):
        self.ic_u = kwargs.pop('ic_u')
        self.ic_v = kwargs.pop('ic_v')
        self.ic_V = kwargs.pop('ic_V')
        self.sigma = kwargs.pop('sigma')
        self.find_mesh = kwargs.pop('find_mesh')
        self.scale_sigma = kwargs.pop('scale_sigma')
        super(ic,self).__init__(*args, **kwargs)
    
    def eval(self, values, x):
        if self.scale_sigma == False:
            values[0] = self.ic_u + np.random.normal(scale=self.sigma)
            values[1] = self.ic_v + np.random.normal(scale=self.sigma)
            values[2] = self.ic_V + np.random.normal(scale=self.sigma)
        else:
            sigma_u = self.ic_u * self.scale_sigma
            sigma_v = self.ic_v * self.scale_sigma
            sigma_V = self.ic_V * self.scale_sigma
            values[0] = self.ic_u + np.random.normal(scale=sigma_u)
            values[1] = self.ic_v + np.random.normal(scale=sigma_v)
            values[2] = self.ic_V + np.random.normal(scale=sigma_V)
    
    def value_shape(self):
        return(3,)        


# Class holding information regarding which surfaces are the bud-scars,
# and what values the parameters should take in those bud-scars.
# Note that there is a difference between the older bud-scars
# and the newest, as the newest has a GAP sitting in it. This class
# holds the indices for the budding scars, and two model parameter
# objects, one for the newer bud-scars and one for the older bud-scars.
# For the newest bud-scar, this class also holds information regarding
# the bud-genes proteins which circle the latest bud scar. 
# Args:
#     decrease_k2, the decrease factor for k2
#     increase_k_2, the increase factor for k_2
#     increase_k2, the increase in k2 around the latest bud-scar
#     increase_k3, the increase of k3 in the ring 
#     i_new, the index for the newest budding scar
#     i_old, a list of the indices for the old budding scars
#     i_ring
#     d, the value of the relative diffusion in the bud-scar 
class ParamBudScars:
    def __init__(self, decrease_k2, increase_k_2, increase_k2, increase_k3, i_new, i_old, i_ring, param_stand, i_barrier=[],
                 decrease_DA_barrier=0.0, d=110, alt_ring=False):
        # The indices for physical surfaces holding bud-scars and the ring 
        self.i_old = i_old
        self.i_new = i_new
        self.i_barrier = i_barrier
        
        # In case we have double ring from creating mesh
        if type(i_ring) is list:
            self.i_ring = i_ring
        else:
            self.i_ring = [i_ring]
        
        # Standard parameter values
        k2_standard = param_stand.k2
        k_2_standard = param_stand.k_2
        k3_standard = param_stand.k3
        DA_standard = param_stand.DA
        
        # Create parameter objects 
        k2_new_scar = k2_standard * decrease_k2
        k_2_new_scar = k_2_standard * increase_k_2
        k2_new_ring = k2_standard * increase_k2
        k3_new_ring = k3_standard * increase_k3
        
        # Ensure that correct parameter values are used for example for gamma
        if alt_ring == False:
            param_dim_old = copy.deepcopy(param_stand)
            param_dim_old.k_2 = k_2_new_scar; param_dim_old.d = d
            param_dim_new = copy.deepcopy(param_stand)
            param_dim_new.k_2 = k_2_new_scar; param_dim_new.d = d
            param_dim_ring = copy.deepcopy(param_stand)
            param_dim_ring.k2 = k2_new_ring; param_dim_ring.k3 = k3_new_ring; param_dim_ring.d = d        
            self.param_old = make_param_non_dimension(param_dim_old, find_ss=False)
            self.param_new = make_param_non_dimension(param_dim_new, find_ss=False)
            self.param_ring = make_param_non_dimension(param_dim_ring, find_ss=False)
        else:
            param_dim_old = copy.deepcopy(param_stand)
            param_dim_old.k_2 = k_2_new_scar; param_dim_old.d = d
            param_dim_new = copy.deepcopy(param_stand)
            param_dim_new.k_2 = k_2_new_scar; param_dim_new.d = d
            param_dim_ring = copy.deepcopy(param_stand)
            param_dim_ring.k2 = k2_new_ring; param_dim_ring.k3 = k3_new_ring; param_dim_ring.d = d; param_dim_ring.k_2 = k_2_new_scar
            self.param_old = make_param_non_dimension(param_dim_old, find_ss=False)
            self.param_new = make_param_non_dimension(param_dim_new, find_ss=False)
            self.param_ring = make_param_non_dimension(param_dim_ring, find_ss=False)
        
        # For barrier around bud-scars
        DA_barrier = decrease_DA_barrier * DA_standard
        param_dim_barrier = ModelParametersWithDim(DA = DA_barrier, d=d)
        self.param_barrier = make_param_non_dimension(param_dim_barrier, find_ss=False)
        
        # For tagging files with the change in k2
        self.decrease_k2 = decrease_k2
        self.increase_k_2 = increase_k_2
        self.increase_k2 = increase_k2
        self.increase_k3 = increase_k3
        self.decrease_DA_barrier = decrease_DA_barrier


def number_to_str(number):
    return str(round(number, 2)).replace(".", "P") 


# Function that will convert the parameters with dimension to
# parameters without dimensions.
# Args:
#     param_dim, parameter object with dimensions
#     find_ss, whether or not a steady state should be found
#       not the case for parameters in subdomain 
# Returns:
#     param_non_dim, parameter object without dimensions
def make_param_non_dimension(param_dim, find_ss=True):
    
    # Express radius in dm to get units correct
    R = param_dim.R
    k1, k_1 = param_dim.k1, param_dim.k_1
    k2, k_2 = param_dim.k2, param_dim.k_2
    k3, k_max = param_dim.k3, param_dim.k_max
    DA, DI, DG = param_dim.DA, param_dim.DI, param_dim.DG
    G0 = param_dim.G0
    
    # Non dimensional parameters 
    c1 = (k1/k_2 * np.sqrt(k_2/k3) / R) * 2.5
    c_1 = k_1/k_2
    c2 = k2/k_2 
    c_max = k_max * np.sqrt(k3/k_2)
    d = DI / DA
    gamma = ((R**2) * k_2) / DA
    V0_init = G0 * np.sqrt(k3/k_2) * (R / 2.5) 
    D = DG / DA
    scale_sigma = param_dim.scale_sigma
    time_conversion_factor = R**2 / (DA*60) 
    
    # Create temporary param object for finding steady state
    if find_ss:
        param_temp = ModelParametersNonDim(c1=c1, c_1=c_1, c2=c2, c_max=c_max,
                                           d=d, gamma=gamma, V0_init=V0_init, D=D)
        ss, info_ss = Find_steady_state.find_init_u_and_v(param_temp)
        if ss == False:
            print("Could not find any steady state -> Will discard parameters")
            return False
        u0, v0 = ss
    else:
        u0 = 0.0
        v0 = 0.0
        info_ss = "turing"
    
    param_mod = ModelParametersNonDim(c1=c1, c_1=c_1, c2=c2, c_max=c_max,
                                      d=d, gamma=gamma, V0_init=V0_init, D=D,
                                      u0=u0, v0=v0, sym_break_cond=info_ss,
                                      scale_sigma=scale_sigma,
                                      time_conversion_factor=time_conversion_factor)
    return param_mod    


# Function that will produce a tag naming the parameters
# that are not standard values for a specific experiment.
# Note that the tag is produced in terms of the
# parameters with dimension.  
# Args:
#     mod_param, the model parameter object
#     param_bud_scar, the parameters for the bud-scar 
# Returns:
#     tag_param, the parameter value tag 
def calc_tag_param(param_dim, param_bud_scar=False):
    
    # The standard non-dimensional parameters
    param_dim_stand = ModelParametersWithDim()
    
    # User provided standard values for the parameters 
    k1_standard = param_dim_stand.k1
    k_1_standard = param_dim_stand.k_1
    k2_standard = param_dim_stand.k2
    k_2_standard = param_dim_stand.k_2
    k3_standard = param_dim_stand.k3
    k_max_standard = param_dim_stand.k_max
    DA_standard = param_dim_stand.DA
    DI_max_standard = param_dim_stand.DI
    DG_init_standard = param_dim_stand.DG
    R_standard = param_dim_stand.R
    G0_standard = param_dim_stand.G0
    d_standard = param_dim_stand.d
    sigma_scale_standard = param_dim_stand.scale_sigma
    
    tag = ""
    if param_dim.d !=  d_standard:
        tag += "d" + number_to_str(param_dim.d)
    if param_dim.k1 !=  k1_standard:
        tag += "k1" + Solve_rd_system.number_to_str(param_dim.k1 / k1_standard)
    if param_dim.k_1 !=  k_1_standard:
        tag += "k_1" + number_to_str(param_dim.k_1 / k_1_standard)
    if param_dim.k2 !=  k2_standard:
        tag += "k2" + number_to_str(param_dim.k2 / k2_standard)
    if param_dim.k_2 !=  k_2_standard:
        tag += "k_2" + number_to_str(param_dim.k_2 / k_2_standard)
    if param_dim.k3 !=  k3_standard:
        tag += "k3" + number_to_str(param_dim.k3 / k3_standard)
    if param_dim.k_max !=  k_max_standard:
        tag += "k_max" + number_to_str(param_dim.k_max / k_max_standard)
    if param_dim.R !=  R_standard:
        tag += "R" + number_to_str(param_dim.R / R_standard)
    if param_dim.k1 !=  k1_standard:
        tag += "a" + number_to_str(param_dim.k1 / k1_standard)
    if param_dim.scale_sigma !=  sigma_scale_standard:
        tag += "ssig" + number_to_str(param_dim.scale_sigma)
    if param_dim.DG != DG_init_standard:
        tag += "DG" + number_to_str(param_dim.DG / DG_init_standard)
    if param_dim.G0 != G0_standard:
        tag += "G0" + number_to_str(param_dim.G0 / G0_standard)
    
    # Check if there are bud-scars
    if param_bud_scar != False:
        tag += "_obud_k2" + number_to_str(param_bud_scar.decrease_k2)
        tag += "_k_2" + number_to_str(param_bud_scar.increase_k_2)
        tag += "_ring_k2" + number_to_str(param_bud_scar.increase_k2)
        tag += "_k3" + number_to_str(param_bud_scar.increase_k3)
        if param_bud_scar.decrease_DA_barrier != 1.0:
            tag += "_DA" + number_to_str(param_bud_scar.decrease_DA_barrier)
    
    return tag        


# Function that converts the model parameters into FeniCS-constants
# Args:
#     mod_param, a model parameters objects
# Returns:
#     c1, c_1, c2, gamma, d, D, c_max as FeniCS constants 
def convert_to_fenics_constants(mod_param):
    
    c1 = Constant(mod_param.c1)
    c_1 = Constant(mod_param.c_1)
    c2 = Constant(mod_param.c2)
    gamma = Constant(mod_param.gamma)
    d = Constant(mod_param.d)
    D = Constant(mod_param.D)
    c_max= Constant(mod_param.c_max)
    
    return c1, c_1, c2, gamma, d, D, c_max


# Function that will create the param_list for the none-budding scar experimental design.
# The function returns all combinations of the parameters to run 
# Args:
#     k2_change_list, the factors to change k2 with
#     k_2_change_list, the factors to change k_2 with
#     r_list, the factors to change the radius with
# Returns:
#     param_list, the list with the parameter objects
def calc_param_list_no_bud_scar(k2_change_list, k_2_change_list, r_list, d=100):
    
    # Get the standard values
    param_standard = ModelParametersWithDim()
    k2_standard = param_standard.k2
    k_2_standard = param_standard.k_2
    
    # Calculate values for k2 and k_2
    k2_val, k_2_val = [], []
    for k2 in k2_change_list:
        k2_val.append(k2 * k2_standard)
    for k_2 in k_2_change_list:
        k_2_val.append(k_2 * k_2_standard)
    
    # Create the parameters list, order k2, k_2 and r 
    param_pre_list = [k2_val, k_2_val, r_list]
    param_list = list(itertools.product(*param_pre_list))
    
    # Create the list with parameter objects
    param_obj_list = []
    for param in param_list:
        param_no_dim = ModelParametersWithDim(k2=param[0], k_2=param[1], R=param[2], d=d)
        param_obj_list.append(param_no_dim)
    
    return param_obj_list


# Function that calculates the bud-scar parameter change list
# when attempting to find the parameters for the bud scar.
# Args:
#    decrease_k2_list, the decrease factor in k2 for the bud-scars
#    increase_k_2_list, the increase factor in k2 for the bud-scars
#    increase_k2_list, the increase in k2 for the ring
#    increase_k3_list, the increase in k3 for the ring
#    d, relative diffusion in the bud-scar 
# Returns:
#    param_bud_scar_list, the list with bud-scar parameter objects
def create_bud_scar_list(decrease_k2_list, increase_k_2_list, increase_k2_list, increase_k3_list, param_standard_list, d=110,
                         i_new=3, i_old=[2], i_ring=4, i_barrier=[], decrease_DA_list=[1.0], alt_ring=False):
    
    print("i_new = {}".format(i_new))
    print("i_ring = {}".format(i_ring))
    
    param_pre_list = [decrease_k2_list,
                      increase_k_2_list,
                      increase_k2_list,
                      increase_k3_list,
                      decrease_DA_list]
    bud_scar_list = list(itertools.product(*param_pre_list))
    param_bud_scar_list = []
    for param_standard in param_standard_list:
        for param in bud_scar_list:
            param_bud = ParamBudScars(decrease_k2=param[0],
                                      increase_k_2=param[1],
                                      increase_k2=param[2],
                                      increase_k3=param[3],
                                      i_new=i_new, i_old=i_old, i_ring=i_ring,
                                      param_stand=param_standard,
                                      i_barrier=i_barrier, decrease_DA_barrier=param[4],
                                      d=d,
                                      alt_ring=alt_ring)
            param_bud_scar_list.append(param_bud)
    
    return param_bud_scar_list

