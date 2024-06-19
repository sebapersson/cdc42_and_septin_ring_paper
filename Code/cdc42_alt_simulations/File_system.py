""""
This file contains functions and classes for organizing the file-system. 

More specifically it contains the FileLocations class and helper functions for this. 
A FileLocations-object contains all the file-locations for a model-run. More specifically:

1. Path to mesh-files (msh-files, xdmf-files). 
2. Path to directory where to store results. The directory in which results 
   are stored is tagged by nucleus, number of bud-scars, parameters for bud-scars 
   and which parameter values are changed from the standard values. For pvd-files the 
   file-location to the result folder is created so that each run with the same parameters 
   get a unique directory under a parent directory. Note, each run gets a result directory in 
   the intermediate and results folders. 
"""


import Model_parameters
import os
import numpy as np 
import pandas as pd
import FEM 


# Class holding the file-locations for a model.
# Args:
#    nucleus (bool), if the mesh has a nucleus
#    mod_param, the model parameters 
#    n_bud_scar="zero", the number of bud-scars
#    find_mesh, whether or not the run is investigating different mesh sizes
#    param_bud_scar, parameters for the bud-scars
#    tag_exp, a tag for the experiment for better saving the result
#    tag_solver, tag for the linear solver used for the equation system
#    fem_opt, fem options for the experiment 
class FileLocations:
    def __init__(self, nucleus, mod_param_dim, n_bud_scar="zero",
                 find_mesh=False, param_bud_scar=False, tag_exp="", fem_opt=FEM.FemOpt()):
        self.nucleus = nucleus
        self.n_bud_scar = n_bud_scar
        
        # Add if there is a nucleus or not
        if nucleus:
            tag_nuc = "tn_"
        else:
            tag_nuc = "fn_"
        
        # If finding mesh mode is one
        self.find_mesh = False
        if find_mesh != False:
            self.find_mesh = True
            dens_tag = "_" + str(find_mesh).replace(".", "P")
            tag_mesh = "Find_mesh/Dens" + dens_tag + "/"
        else:
            tag_mesh = ""
        
        # Tag for number of bud-scars
        tag_bud_scars = 'b' + n_bud_scar
        
        # Mesh file-location paths
        if find_mesh != False:
            self.msh_file = "./Generate_mesh/Find_mesh/Mesh_" + tag_nuc + tag_bud_scars + dens_tag + ".msh"
            self.mesh_folder = ("../../Intermediate/Mesh/Find_mesh/Mesh_" + tag_nuc
                                + tag_bud_scars + dens_tag + "/")
        else:
            self.msh_file = "./Generate_mesh/Mesh_" + tag_nuc + tag_bud_scars + ".msh"
            self.mesh_folder = ("../../Intermediate/Mesh/" + tag_nuc
                                + tag_bud_scars + "/")
        
        # The tag for different parameters
        tag_param = Model_parameters.calc_tag_param(mod_param_dim, param_bud_scar)
        if tag_param != "":
            temp = tag_param
            tag_param = "_" + temp
        
        # For the fem-opt solver
        if fem_opt.solver != 'gmres':
            tag_solver = "_" + fem_opt.solver
        else:
            tag_solver = ""
        if fem_opt.pre_cond != 'ilu':
            tag_pre_cond = "_" + fem_opt.pre_cond
        else:
            tag_pre_cond = ""
        if fem_opt.tol_step_diff != 10:
            tag_step_diff = "sdiff" + Model_parameters.number_to_str(fem_opt.tol_step_diff)
        else:
            tag_step_diff = ""
        if fem_opt.seed != False:
            tag_seed = "_seed" + Model_parameters.number_to_str(fem_opt.seed)
        else:
            tag_seed = ""
        tag_fem_opt = tag_solver + tag_pre_cond + tag_seed + tag_step_diff
        
        # Result folders 
        self.result_folder = ("../../Result/" + tag_exp + tag_mesh + tag_nuc +
                              tag_bud_scars + tag_param + tag_fem_opt + "/")
        # Ensure exists
        if not os.path.exists(self.result_folder):
            os.makedirs(self.result_folder)
        self.pvd_folder = self.result_folder + "pvd_folder/"
        if not os.path.exists(self.pvd_folder):
            os.mkdir(self.pvd_folder)
        
        # Intermediate result folders
        self.intermediate_folder = ("../../Intermediate/Experiments/" + tag_exp + tag_mesh 
                                    + tag_nuc + tag_bud_scars + tag_param + tag_fem_opt + "/")
        if not os.path.exists(self.intermediate_folder):
            os.makedirs(self.intermediate_folder)
        self.intermediate_file = self.intermediate_folder + "Simulation_data.csv"
        
        # The case pole-data is written to file
        if fem_opt.write_pole_data != False:
            self.pole_data = self.intermediate_folder + "Pole_data.csv"
            # Ensure that the correct index is used when writing to pole_data
            self.pole_data_index = self.intermediate_folder + "Pole_data_index.csv"
            self.value_write_pole = 1
        else:
            self.pole_data = False
            self.pole_data_index = False
            self.value_write_pole = 1
        
        # If the u-surface values should be saved in the end
        if fem_opt.save_u_surface != False:
            self.u_surface = self.intermediate_folder + "u_surface.csv"
        else:
            self.u_surface = ""
    
    # Method that correctly calculates the pole-index value
    # when writing to save pole-index file
    def calc_pole_index(self):
        if not os.path.isfile(self.pole_data_index):
            value_write = 1
            data_to_save = pd.DataFrame({"Value": [value_write]})
            data_to_save.to_csv(self.pole_data_index)
        else:
            data_file = pd.read_csv(self.pole_data_index)["Value"]
            max_val = np.max(data_file)
            value_write = max_val + 1
            data_to_save = pd.DataFrame({"Value": [value_write]})
            data_to_save.to_csv(self.pole_data_index, header=False, mode='a')
        self.value_write_pole = value_write
    
    # Method that will print the file folders
    def print_dir_loc(self):
        print("Intermediate folder = {}".format(self.intermediate_folder))
        print("Result folder = {}".format(self.result_folder))
        print("Mesh folder = {}".format(self.mesh_folder))


# Function that will, for a run, create the correct pvd-folder.
# Note that the folders are named as 1, 2, etc..., where
# the number represents the run (multiple runs are used
# as the rd-system is stochastic).
# Args:
#     file_loc, file-locations object for a model
# Returns:
#     pvd_folder, path to the pvd-folder
def create_pwd_folder(file_loc):
    path_pvd_files = file_loc.pvd_folder
    dirs_pvd_dir = os.listdir(path_pvd_files)
    i_str = str(file_loc.value_write_pole)
    
    pvd_folder = path_pvd_files + i_str + '/'
    return pvd_folder
