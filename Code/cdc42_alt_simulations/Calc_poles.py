'''
This file contains functions for calculating the number 
of poles using the heruistic DB-scan algorithm (see supplementary part paper). 
Also, it contains files for calculating pole-information, and writing pole-results 
to file. 

Note, the standard settings for the DB-scan are adapted to the specific 
geometry of a unit-sphere. 
'''


import numpy as np
import pandas as pd
from dolfin import *
from sklearn.cluster import DBSCAN
from sklearn.metrics import accuracy_score
import os
import Config 


# Class that holds pole-data, that is data on the progression on forming a pole.
# The pole-data object also has a method allowing it to write the pole data
# to a file, where the result is appended if the file already exists. Note,
# this data is required 
# Args:
#     t, the current time-value 
#     u_max, the maximum u-concentration
#     pole_surface_ratio, the pole surface ratio
#     pole_surface_ratio_hard, the inner pole surface ratio
#     pole_percent_binary, the percent that is binary on the surface 
class PoleData:
    def __init__(self, t, u_max, pole_surface_ratio, pole_surface_ratio_hard, pole_percent_binary):
        self.t = t
        self.u_max = u_max
        self.pole_surface_ratio = pole_surface_ratio
        self.pole_surface_ratio_hard = pole_surface_ratio_hard
        self.pole_percent_binary = pole_percent_binary


def write_pole_to_file(pole_data, file_loc):
    # Write to file, if file exists append
    should_terminate = False
    
    Config.lock.acquire()
    try:
        data_save = pd.DataFrame({"t": [pole_data.t],
                                  "u_max": [pole_data.u_max],
                                  "ratio_hard": [pole_data.pole_surface_ratio_hard],
                                  "ratio": [pole_data.pole_surface_ratio],
                                  "pre_binary": [pole_data.pole_percent_binary],
                                  "index" : [file_loc.value_write_pole]})
        if not os.path.isfile(file_loc.pole_data):
            data_save.to_csv(file_loc.pole_data)
        else:
            data_save.to_csv(file_loc.pole_data, header=False, mode='a')
    except:
        should_terminate = True
    Config.lock.release()
    
    return should_terminate


# Function that calculates pole data at a given time-step t.
# More specifically, ration pole is returned (both small
# and big ratio) och the percentage of nodes that
# are segregated in the simulations. The end result
# is written to a file. Note that print pole data option
# must be provided for the fem-options for this
# function to actually be triggered. 
# Args:
#    U, the values of the current solution
#    bV, the boundary mesh for the cell surface 
#    fem_opt, the options for the FEM-solution
#    t, the current time-step 
# Returns:
#    ratios and percentage of nodes that are covered with pole 
def calc_pole_data(U, bV, fem_opt, t):
    # Calculate number of min and max-dofs
    _u, _v, _V = U.split()
    u_mem = interpolate(_u, bV)
    
    # See if data is segmented
    u_vec = u_mem.vector()
    u_min = u_vec.min()
    u_max = u_vec.max()
    u_len = len(u_vec)
    num_max_dofs = 0
    num_max_dofs_hard = 0
    num_min_dofs = 0
    # Tolerances if concentration is sufficiently close to minimum 
    min_TOL = fem_opt.binary_percent_pole_min*(u_max - u_min)
    max_TOL = fem_opt.binary_percent_pole_max*(u_max - u_min)
    max_TOL_hard = 0.2095*(u_max - u_min)
    # Loop through and see number of nodes that are minimum or maximum
    for index_u in range(0, u_len):
        if abs(u_vec[index_u] - u_max) < max_TOL:
            num_max_dofs += 1
        elif abs(u_vec[index_u] - u_min) < min_TOL:
            num_min_dofs += 1
        
        # Hard max limit (size of most intense part)
        if abs(u_vec[index_u] - u_max) < max_TOL_hard:
            num_max_dofs_hard += 1
    
    pole_data = PoleData(t=t,
                         u_max=u_max,
                         pole_surface_ratio=100.0*num_max_dofs/u_len,
                         pole_surface_ratio_hard=100.0*num_max_dofs_hard/u_len,
                         pole_percent_binary=(num_max_dofs + num_min_dofs) / u_len)
    return pole_data        


# Function that given a data-frame on format x, y, z, u will 
# calculate the number of poles using the dbscan cluster algorithm.
# For the clustering, the euclidian norm is used, with a minimum 
# distance eps (can be set by the user) and number of neighbours 
# being equal to 5. Note, the data should only be for one 
# individual. 
# Args:
#     data, the data to cluster 
#     u_max_filt, the filtering level for u_max, should be in [0, 1]
#     eps, the neighbour tolerance for dbscan 
#     min_samples, the number of minimal neighbours for the clustering
# Returns:
#     n_poles, the number of poles with current parameters 
def calc_n_poles(data, u_max_filt, eps=0.1, min_samples=5):
    # Filter the data to get clusters 
    u_max = np.max(data["u"])
    filter_limit = u_max * u_max_filt
    data_filt = data[data["u"] > filter_limit]
    
    # np-array is required for dbscan 
    data_array = data_filt[["x", "y", "z"]].to_numpy()
    # Cluster and get the number of poles 
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(data_array)
    db_labels = db.labels_
    n_poles = len(set(db_labels)) - (1 if -1 in db_labels else 0) 
    
    return n_poles 


# Function that checks the number of poles by a heuristic-scheme. More precisely, 
# if number of poles is calculated to 1 using the standard value of 
# u_max = 0.55, then u_max is increaed to 0.90 and decreased to 0.40. 
# This aims to avoid the special cases, where one pole might be very
# weak, or two poles are almost joined togehter. 
# Args:
#     data, the data to cluster 
#     eps, the neighbour tolerance for dbscan 
#     min_samples, the number of minimal neighbours for the clustering
# Returns:
#     n_poles, the number of poles, where 1 is returned if there is not a
#       weak pole or two poles are that are almost joined 
def calc_n_poles_her(data, eps=0.1, min_samples=5):
    u_max_filt = 0.55
    n_poles = calc_n_poles(data, u_max_filt, eps=eps, min_samples=min_samples)
    
    # See if should trigger heuristic 
    if n_poles == 1:
        u_max_inc = 0.90
        n_poles_inc = calc_n_poles(data, u_max_inc, eps=eps, 
                                   min_samples=min_samples)
        u_max_dec = 0.4
        n_poles_dec = calc_n_poles(data, u_max_dec, eps=eps, 
                                   min_samples=min_samples)
        if n_poles_inc != n_poles:
            return n_poles_inc
        elif n_poles_dec != n_poles:
            return n_poles_dec
        
    return n_poles
