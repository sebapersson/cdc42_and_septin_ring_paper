from dolfin import *
import meshio
import os
import sys
import itertools 
import copy 
import numpy as np
import pandas as pd
import pickle
from sklearn.cluster import DBSCAN

class FileLocations:
    def __init__(self, mesh_name, save_name, tag_save="_", name_store="", name_save=""):
        
        self.msh_file = "./Mesh_files/" + mesh_name + ".msh"
        self.mesh_folder = "../../Intermediate/Mesh/" + mesh_name + "/"
        self.save_name = save_name
        self.mesh_name = mesh_name
        
        # Directory for storing csv-files 
        self.dir_intermediate = "../../Intermediate/Experiments/" + save_name + tag_save + "/"
        if not os.path.isdir(self.dir_intermediate):
            os.makedirs(self.dir_intermediate)
        
        # Keeping track of which index a run has 
        self.pole_data_index = self.dir_intermediate + "Pole_data_index.csv"
        self.run_index = 1
        self.tag_save = tag_save
        
        # Folder storing pvd-files
        self.dir_pvd = ""

        # Folder where the end result is stored, or read from 
        self.name_store = name_store
        self.name_save = name_save
        if name_save == "":
            self.dir_save_run = "../../Intermediate/Experiments/Saved_runs/" + name_store + "/"
        else:
            self.dir_save_run = "../../Intermediate/Experiments/Saved_runs/" + name_save + "/"
        self.dir_read_run = "../../Intermediate/Experiments/Saved_runs/" + name_store + "/"
        if not os.path.isdir(self.dir_save_run):
            os.makedirs(self.dir_save_run)
    
    # Method that correctly calculates the pole-index value, allowing
    # tracking of simulations 
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
        self.run_index = value_write
        self.dir_pvd = "../../Results/Pvd_folders/" + self.save_name + self.tag_save + "/" + self.name_store + "/" + str(self.run_index) + "/"
    
    
    def process_msh_file_to_xdmf_2d(self):
        
        # For compacter code 
        dir_save_mesh = self.mesh_folder
        path_msh_file = self.msh_file
        
        if not os.path.isdir(dir_save_mesh):
            os.makedirs(dir_save_mesh)
        
        if not os.path.isfile(path_msh_file):
            print("Error: Msh-file does not exist")
            sys.exit(1)
        
        # Conversion using meshio version 3.2
        msh = meshio.read(path_msh_file)
        meshio.write(dir_save_mesh + "mesh.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}))
        # Write the triangle subdomains 
        meshio.write(dir_save_mesh + "subdomains.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                                                              cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
    
    
    def read_mesh_2d(self):
        
        mesh_folder = self.mesh_folder
        
        # Reading the mesh 
        mesh = Mesh()
        with XDMFFile(mesh_folder + "mesh.xdmf") as infile:
            infile.read(mesh)
        # Subdomains
        sub_domains_pre = MeshValueCollection("size_t", mesh, 3)
        with XDMFFile(mesh_folder + "subdomains.xdmf") as infile:
            infile.read(sub_domains_pre, "name_to_read")
        subdomains = cpp.mesh.MeshFunctionSizet(mesh, sub_domains_pre)
        
        return mesh, subdomains
    
    
    def process_msh_file_to_xdmf_3d(self):
        
        dir_save_mesh = self.mesh_folder
        path_msh_file = self.msh_file
        
        if not os.path.isdir(dir_save_mesh):
            os.makedirs(dir_save_mesh)
        
        if not os.path.isfile(path_msh_file):
            print("Error: Msh-file does not exist")
            sys.exit(1)
        
        # Conversion using meshio version 3.2
        msh = meshio.read(path_msh_file)
        meshio.write(dir_save_mesh + "mesh.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]}))
        meshio.write(dir_save_mesh + "subdomains.xdmf",
                     meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]},
                                 cell_data={"tetra": {"name_to_read": msh.cell_data["tetra"]["gmsh:physical"]}}))
        meshio.write(dir_save_mesh + "surfaces.xdmf",
                     meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                                 cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
    
    
    def read_mesh_3d(self):
        
        mesh_folder = self.mesh_folder
        
        # Reading the mesh 
        mesh = Mesh()
        with XDMFFile(mesh_folder + "mesh.xdmf") as infile:
            infile.read(mesh)
        # Subdomains
        sub_domains_pre = MeshValueCollection("size_t", mesh, 3)
        with XDMFFile(mesh_folder + "subdomains.xdmf") as infile:
            infile.read(sub_domains_pre, "name_to_read")
        subdomains = cpp.mesh.MeshFunctionSizet(mesh, sub_domains_pre)
        surface_domain_pre = MeshValueCollection("size_t", mesh, 2)
        with XDMFFile(mesh_folder + "surfaces.xdmf") as infile:
            infile.read(surface_domain_pre, "name_to_read")
        surfaces = cpp.mesh.MeshFunctionSizet(mesh, surface_domain_pre)
        
        return mesh, subdomains, surfaces

class ModelParametersOkada:
    def __init__(self, k1a=10, k1b=10, k2=0.009, k3=0.0001, k4=0.024, k5=0.024, k6a=10, k6b=10, k7=0.021,
                 k8a=0.9, k8b=0.13, k9a=1.5, k9b=0.5, k10=10, k11=0.5, k12a=10, k12b=6.2, k13a=1, k13b=0.62,
                 k14=0.02, k15=0.05, k16=0.05, k17=0.1, k18=0.05, k19=0.01, k20=0.024, k23=2.0, k24a=10, k24b=0.5,
                 k25=1.0, k26a=10, k26b=0.5, k27=1.0, r = 2.5, L = 2.5,
                 Cdc42_tot=5.0, Cdc24_tot = 0.017, F_tot=0.06, GapS_tot=0.005, S_tot=0.2, I_tot=5.0,
                 Dm=0.01, Dm_gdp=0.01, Dms = 0.0025, Dmss = 0.00025, Kt=1300, crowding=False, barrier=False, Kd=150):
        
        self.k1a = k1a
        self.k1b = k1b
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5
        self.k6a = k6a
        self.k6b = k6b
        self.k7 = k7
        self.k8a = k8a
        self.k8b = k8b
        self.k9a = k9a
        self.k9b = k9b
        self.k10 = k10
        self.k11 = k11
        self.k12a = k12a
        self.k12b = k12b
        self.k13a = k13a
        self.k13b = k13b
        self.k14 = k14
        self.k15 = k15
        self.k16 = k16
        self.k17 = k17
        self.k18 = k18
        self.k19 = k19
        self.k20 = k20
        self.k23 = k23
        self.k24a = k24a
        self.k24b = k24b
        self.k25 = k25
        self.k26a = k26a
        self.k26b = k26b
        self.k27 = k27
        
        # If septin crowds the membrane 
        self.Kt = Kt
        self.crowding = crowding

        self.Kd = Kd
        self.barrier = barrier
        
        self.Dm = Dm / L**2
        self.Dms = Dms / L**2
        self.Dmss = Dmss / L**2
        self.Dm_gdp = Dm_gdp / L**2
        
        dx = 10e-9
        rad_m = r * 1e-6
        self.eta = ((rad_m + dx)**3 - rad_m**3) / rad_m**3
        self.L = L
        self.r = r
        
        self.Cdc42_tot = Cdc42_tot
        self.Cdc24_tot = Cdc24_tot
        self.F_tot = F_tot
        self.GapS_tot = GapS_tot
        self.S_tot = S_tot
        self.I_tot = I_tot
        
        self.ic_Cdc42T = 0.0
        self.ic_Cdc42D = 0.0
        self.ic_Cdc42I = 0.0
        self.ic_BemGEF42 = 0.0
        self.ic_Cdc24 = 0.0
        self.ic_GapS = 0.0
        self.ic_F = 0.0
        self.ic_Cdc42TF = 0.0
        self.ic_FS = 0.0
        self.ic_S = 0.0
        self.ic_P = 0.0
        
        self.ic_Cdc42Ic = 0.0
        self.ic_Cdc24c = 0.0
        self.ic_GapSc = 0.0
        self.ic_Fc = 0.0
        self.ic_Sc = 0.0
        self.ic_FSc = 0.0
        
        # Compensate and ensure Cdc42c is reduced
        self.tag = ""
        
    def make_param_constants(self):
        
        k1a = Constant(self.k1a)
        k1b = Constant(self.k1b)
        k2 = Constant(self.k2)
        k3 = Constant(self.k3)
        k4 = Constant(self.k4)
        k5 = Constant(self.k5)
        k6a = Constant(self.k6a)
        k6b = Constant(self.k6b)
        k7 = Constant(self.k7)
        k8a = Constant(self.k8a)
        k8b = Constant(self.k8b)
        k9a = Constant(self.k9a)
        k9b = Constant(self.k9b)
        k10 = Constant(self.k10)
        k11 = Constant(self.k11)
        k12a = Constant(self.k12a)
        k12b = Constant(self.k12b)
        k13a = Constant(self.k13a)
        k13b = Constant(self.k13b)
        k14 = Constant(self.k14)
        k15 = Constant(self.k15)
        k16 = Constant(self.k16)
        k17 = Constant(self.k17)
        k18 = Constant(self.k18)
        k19 = Constant(self.k19)
        k20 = Constant(self.k20)
        k23 = Constant(self.k23)
        k24a = Constant(self.k24a)
        k24b = Constant(self.k24b)
        k25 = Constant(self.k25)
        k26a = Constant(self.k26a)
        k26b = Constant(self.k26b)
        k27 = Constant(self.k27)
        Dm = Constant(self.Dm)
        Dm_gdp = Constant(self.Dm_gdp)
        Dms = Constant(self.Dms)
        Dmss = Constant(self.Dmss)
        
        return (k1a, k1b, k2, k3, k4, k5, k6a, k6b, k7, k8a, k8b, k9a, k9b, k10, k11, k12a, k12b, k13a, k13b,
                k14, k15, k16, k17, k18, k19, k20, k23, k24a, k24b, k25, k26a, k26b, k27, Dm, Dms, Dmss, Dm_gdp)
        
    def return_param(self):
        return (self.k1a, self.k1b, self.k2, self.k3, self.k4, self.k5, self.k6a, self.k6b, self.k7, self.k8a, self.k8b, self.k9a, self.k9b, self.k10,
                self.k11, self.k12a, self.k12b, self.k13a, self.k13b, self.k14, self.k15, self.k16, self.k17, self.k18, self.k19,
                self.k20, self.k23, self.k24a, self.k24b, self.k25, self.k26a, self.k26b, self.k27, self.eta, self.Dm_gdp)


    def calc_save_tag(self, bud_scar_param=None):
        
        def number_to_str(number):
            return str(round(number, 3)).replace(".", "P") 
        
        mod_param_stand = ModelParametersSeptin()
        
        tag = "_"
        if self.k3 != mod_param_stand.k3:
            tag += "k3" + number_to_str(self.k3 / mod_param_stand.k3)
        if self.k4 != mod_param_stand.k4:
            tag += "k4" + number_to_str(self.k4 / mod_param_stand.k4)
        if self.k5 != mod_param_stand.k5:
            tag += "k5" + number_to_str(self.k5 / mod_param_stand.k5)
        if self.k8b != mod_param_stand.k8b:
            tag += "k8b" + number_to_str(self.k8b / mod_param_stand.k8b)
        if self.k14 != mod_param_stand.k14:
            tag += "k14" + number_to_str(self.k14 / mod_param_stand.k14)
        if self.k17 != mod_param_stand.k17:
            tag += "k17" + number_to_str(self.k17 / mod_param_stand.k17)
        if self.k18 != mod_param_stand.k18:
            tag += "k18" + number_to_str(self.k18 / mod_param_stand.k18)
        if self.k19 != mod_param_stand.k19:
            tag += "k19" + number_to_str(self.k19 / mod_param_stand.k19)
        if self.r != mod_param_stand.r:
            tag += "r" + number_to_str(self.r / mod_param_stand.r)
        if self.k25 != mod_param_stand.k25:
            tag += "k25" + number_to_str(self.k25 / mod_param_stand.k25)
        if self.k6a != mod_param_stand.k6a:
            tag += "k6a" + number_to_str(self.k6a / mod_param_stand.k6a)
        if self.k10 != mod_param_stand.k10:
            tag += "k10" + number_to_str(self.k10 / mod_param_stand.k10)
        if self.Cdc24_tot != mod_param_stand.Cdc24_tot:
            tag += "Cdc24_tot" + number_to_str(self.Cdc24_tot / mod_param_stand.Cdc24_tot)

        return tag 

# 2d positive feedback model
class ModelParametersDim:
    def __init__(self, k1a=10, k1b=10, k2a=0.16, k2b=0.35, k3=0.35, k4a=10, k4b=10, k5a=36,
                 k5b=0.65, k7=10, Dm=0.011, Dm_gdp=0.023, Dc=10, eta=0.01, r = 2.5, L=1.0, 
                 ic_Cdc42T=31.35, ic_Cdc42D=21.76, ic_BemGEF42=1.42, ic_BemGEFm=0.04395, ic_BemGEFc=0.001359, ic_Cdc42I=0.39,
                 BemGEF_tot=0.017, Cdc42_tot=1.0, old_way=True, save_run=None, read_run=None, cap_CDC42T=None):
        
        self.old_way = True
        self.save_run = save_run
        self.read_run = read_run

        self.k1a = k1a
        self.k1b = k1b
        self.k2a = k2a
        self.k2b = k2b
        self.k3 = k3
        self.k4a = k4a
        self.k4b = k4b
        self.k5a = k5a
        self.k5b = k5b
        self.k7 = k7
        self.Dm = Dm / L**2
        self.Dm_gdp = Dm_gdp / L**2
        self.Dc = Dc / L**2
        self.L = L
        dx = 10e-9
        rad_m = r * 1e-6
        self.eta = ((rad_m + dx)**3 - rad_m**3) / rad_m**3
        self.r = r
        
        self.ic_Cdc42T = ic_Cdc42T
        self.ic_Cdc42D = ic_Cdc42D
        self.ic_BemGEF42 = ic_BemGEF42
        self.ic_BemGEFm = ic_BemGEFm
        self.ic_BemGEFc = ic_BemGEFc
        self.ic_Cdc42I = ic_Cdc42I
        
        self.Cdc42_tot = Cdc42_tot
        self.BemGEF_tot = BemGEF_tot

        self.cap_CDC42T = cap_CDC42T
        
        # Compensate and ensure Cdc42c is reduced
        self.tag = ""
        
    
    def make_param_constants(self):
        k1a = Constant(self.k1a)
        k1b = Constant(self.k1b)
        k2a = Constant(self.k2a)
        k2b = Constant(self.k2b)
        k3 = Constant(self.k3)
        k4a = Constant(self.k4a)
        k4b = Constant(self.k4b)
        k5a = Constant(self.k5a)
        k5b = Constant(self.k5b)
        k7 = Constant(self.k7) 
        Dm = Constant(self.Dm) 
        Dc = Constant(self.Dc) 
        eta = Constant(self.eta)
        Dm_gdp = Constant(self.Dm_gdp)
        
        return k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, Dm_gdp
    
    def return_param(self):
        return self.k1a, self.k1b, self.k2a, self.k2b, self.k3, self.k4a, self.k4b, self.k5a, self.k5b, self.k7, self.Dm, self.Dc, self.eta, self.Dm_gdp
    
    def calc_save_tag(self, bud_scar_param=None, exocyte_data=None):
        
        def number_to_str(number):
            return str(round(number, 3)).replace(".", "P") 
        
        mod_param_stand = ModelParametersDim()
        
        tag = "_"
        if self.k1a != mod_param_stand.k1a:
            tag += "k1a" + number_to_str(self.k1a / mod_param_stand.k1a)
        if self.k1b != mod_param_stand.k1b:
            tag += "k1b" + number_to_str(self.k1b / mod_param_stand.k1b)
        if self.k2a != mod_param_stand.k2a:
            tag += "k2a" + number_to_str(self.k2a / mod_param_stand.k2a)
        if self.k2b != mod_param_stand.k2b:
            tag += "k2b" + number_to_str(self.k2b / mod_param_stand.k2b)
        if self.k3 != mod_param_stand.k3:
            tag += "k3" + number_to_str(self.k3 / mod_param_stand.k3)
        if self.k4a != mod_param_stand.k4a:
            tag += "k4a" + number_to_str(self.k4a / mod_param_stand.k4a)
        if self.k4b != mod_param_stand.k4b:
            tag += "k4b" + number_to_str(self.k4b / mod_param_stand.k4b)
        if self.k5a != mod_param_stand.k5a:
            tag += "k5a" + number_to_str(self.k5a / mod_param_stand.k5a)
        if self.k5b != mod_param_stand.k5b:
            tag += "k5b" + number_to_str(self.k5b / mod_param_stand.k5b)
        if self.k7 != mod_param_stand.k7:
            tag += "k7" + number_to_str(self.k7 / mod_param_stand.k7)
        if self.Dm != mod_param_stand.Dm:
            tag += "Dm" + number_to_str(self.Dm / mod_param_stand.Dm)
        if self.Dc != mod_param_stand.Dc:
            tag += "Dc" + number_to_str(self.Dc / mod_param_stand.Dc)
        if self.Dm_gdp != mod_param_stand.Dm_gdp:
            tag += "DmGdp" + number_to_str(self.Dm_gdp / mod_param_stand.Dm_gdp)
        if self.r != mod_param_stand.r:
            tag += "r" + number_to_str(self.r / mod_param_stand.r)
        if self.Cdc42_tot != mod_param_stand.Cdc42_tot:
            tag += "Cdc42Tot" + number_to_str(self.Cdc42_tot / mod_param_stand.Cdc42_tot)
        if self.BemGEF_tot != mod_param_stand.BemGEF_tot:
            tag += "BemGefTot" + number_to_str(self.BemGEF_tot / mod_param_stand.BemGEF_tot)

        if bud_scar_param != None:
            tag += "ringk1a" + number_to_str(bud_scar_param.param_ring.k1a / mod_param_stand.k1a)
            tag += "scark2b" + number_to_str(bud_scar_param.param_bud_scar.k2b / mod_param_stand.k2b)

        if self.cap_CDC42T != None:
            tag += "CapCdc42T" + number_to_str(self.cap_CDC42T)

        if exocyte_data != None:
            tag += "Lhit" + number_to_str(exocyte_data.lambda_hit)
            tag += "Lcab" + number_to_str(exocyte_data.lambda_cable)
            tag += "alpha" + number_to_str(exocyte_data.alpha)            
        
        self.tag = tag
        return tag 

class ModelParametersAxl2V1:
    def __init__(self, k1a=10, k1b=10, k2a=0.16, k2b=0.35, k3=0.35, k4a=10, k4b=10, k5a=36,
                 k5b=0.65, k7=10, k12a=0.05, k12b=0.01, k13=0.024, k14=0.024, 
                 k15=1.0, k16=0.05, k17=0.05, k18=0.1, k19=3.0, k20=1.0, k21=0.02, k22=0.5, k23=5.0,
                 Dm=0.01, Dms=0.0025, Dmss=0.00025, Dm_gdp=0.01, Dc=10, eta=0.01, r = 2.5, L=1.0, 
                 BemGEF_tot=0.017, Cdc42_tot=1.0, S_tot = 0.2, GapS_tot=0.005, Kd=150, Axl2_tot=0.2,
                 crowding=False, Kt=1500, exo_axl2=100, recruited=False, use_k23_exp=False, tol_k23=1e-1, 
                 P_pos=False, k17_alt=False, use_gamma=False, gamma_d=0.1, gamma_k=0.3, gamma_push=0.1, 
                 read_run=None, crowding_p=False):

        self.read_run = read_run

        self.crowding_p = crowding_p

        self.k17_alt = k17_alt
        self.P_pos = P_pos
        self.Kd = Kd
        self.crowding=crowding
        self.Kt = Kt
        self.exo_axl2 = exo_axl2
        self.recruited = recruited
        self.use_k23_exp = use_k23_exp
        self.tol_k23 = tol_k23

        self.gamma_d = gamma_d
        self.gamma_k = gamma_k
        self.gamma_push = gamma_push
        self.use_gamma = use_gamma

        self.k1a = k1a
        self.k1b = k1b
        self.k2a = k2a
        self.k2b = k2b
        self.k3 = k3
        self.k4a = k4a
        self.k4b = k4b
        self.k5a = k5a
        self.k5b = k5b
        self.k7 = k7
        self.k12a = k12a
        self.k12b = k12b
        self.k13 = k13
        self.k14 = k14
        self.k15 = k15
        self.k16 = k16
        self.k17 = k17
        self.k18 = k18
        self.k19 = k19
        self.k20 = k20
        self.k21 = k21
        self.k22 = k22
        self.k23 = k23

        self.Dm = Dm / L**2
        self.Dms = Dms / L**2
        self.Dmss = Dmss / L**2
        self.Dm_gdp = Dm_gdp / L**2
        self.Dc = Dc / L**2
        self.L = L
        dx = 10e-9
        rad_m = r * 1e-6
        self.eta = ((rad_m + dx)**3 - rad_m**3) / rad_m**3
        self.r = r
        
        self.ic_Cdc42T = 0.0
        self.ic_Cdc42D = 0.0
        self.ic_BemGEF42 = 0.0
        self.ic_BemGEFm = 0.0
        self.ic_BemGEFc = 0.0
        self.ic_Cdc42I = 0.0
        self.ic_Sc = 0.0
        self.ic_S = 0.0
        self.ic_P = 0.0
        self.ic_GapS = 0.0
        self.ic_GapSc = 0.0
        self.ic_Axl2c = 0.0
        self.ic_Axl2 = 0.0
        self.ic_Axl2_S = 0.0
        
        self.Cdc42_tot = Cdc42_tot
        self.BemGEF_tot = BemGEF_tot
        self.S_tot = S_tot
        self.GapS_tot = GapS_tot
        self.Axl2_tot = Axl2_tot

        # Compensate and ensure Cdc42c is reduced
        self.tag = ""
    
    def make_param_constants(self):
        k1a = Constant(self.k1a)
        k1b = Constant(self.k1b)
        k2a = Constant(self.k2a)
        k2b = Constant(self.k2b)
        k3 = Constant(self.k3)
        k4a = Constant(self.k4a)
        k4b = Constant(self.k4b)
        k5a = Constant(self.k5a)
        k5b = Constant(self.k5b)
        k7 = Constant(self.k7) 
        k12a = Constant(self.k12a)
        k12b = Constant(self.k12b)
        k13 = Constant(self.k13)
        k14 = Constant(self.k14)
        k15 = Constant(self.k15)
        k16 = Constant(self.k16) 
        k17 = Constant(self.k17) 
        k18 = Constant(self.k18) 
        k19 = Constant(self.k19) 
        k20 = Constant(self.k20) 
        k21 = Constant(self.k21) 
        k22 = Constant(self.k22) 
        Dm = Constant(self.Dm) 
        Dms = Constant(self.Dms) 
        Dmss = Constant(self.Dmss) 
        Dc = Constant(self.Dc) 
        eta = Constant(self.eta)
        Dm_gdp = Constant(self.Dm_gdp)
        
        return k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, k12a, k12b, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, Dm, Dms, Dmss, Dc, eta, Dm_gdp
    
    def return_param(self):
        return self.k1a, self.k1b, self.k2a, self.k2b, self.k3, self.k4a, self.k4b, self.k5a, self.k5b, self.k7, self.k12a, self.k12b, self.k13, self.k14, self.k15, self.k16, self.k17, self.k18, self.k19, self.k20, self.k21, self.k22, self.Dm, self.Dms, self.Dmss, self.Dc, self.eta, self.Dm_gdp
    
    def calc_save_tag(self, bud_scar_param=None, exocyte_data=None):
        
        def number_to_str(number):
            return str(round(number, 3)).replace(".", "P") 
        
        mod_param_stand = ModelParametersAxl2V1()
        
        tag = "_"
        if self.k1a != mod_param_stand.k1a:
            tag += "k1a" + number_to_str(self.k1a / mod_param_stand.k1a)
        if self.k1b != mod_param_stand.k1b:
            tag += "k1b" + number_to_str(self.k1b / mod_param_stand.k1b)
        if self.k2a != mod_param_stand.k2a:
            tag += "k2a" + number_to_str(self.k2a / mod_param_stand.k2a)
        if self.k2b != mod_param_stand.k2b:
            tag += "k2b" + number_to_str(self.k2b / mod_param_stand.k2b)
        if self.k3 != mod_param_stand.k3:
            tag += "k3" + number_to_str(self.k3 / mod_param_stand.k3)
        if self.k4a != mod_param_stand.k4a:
            tag += "k4a" + number_to_str(self.k4a / mod_param_stand.k4a)
        if self.k4b != mod_param_stand.k4b:
            tag += "k4b" + number_to_str(self.k4b / mod_param_stand.k4b)
        if self.k5a != mod_param_stand.k5a:
            tag += "k5a" + number_to_str(self.k5a / mod_param_stand.k5a)
        if self.k5b != mod_param_stand.k5b:
            tag += "k5b" + number_to_str(self.k5b / mod_param_stand.k5b)
        if self.k7 != mod_param_stand.k7:
            tag += "k7" + number_to_str(self.k7 / mod_param_stand.k7)
        if self.k12a != mod_param_stand.k12a:
            tag += "k12a" + number_to_str(self.k12a / mod_param_stand.k12a)
        if self.k12b != mod_param_stand.k12b:
            tag += "k12b" + number_to_str(self.k12b / mod_param_stand.k12b)
        if self.k13 != mod_param_stand.k13:
            tag += "k13" + number_to_str(self.k13 / mod_param_stand.k13)
        if self.k14 != mod_param_stand.k14:
            tag += "k14" + number_to_str(self.k14 / mod_param_stand.k14)
        if self.k15 != mod_param_stand.k15:
            tag += "k15" + number_to_str(self.k15 / mod_param_stand.k15)
        if self.k16 != mod_param_stand.k16:
            tag += "k16" + number_to_str(self.k16 / mod_param_stand.k16)
        if self.k17 != mod_param_stand.k17:
            tag += "k17" + number_to_str(self.k17 / mod_param_stand.k17)
        if self.k18 != mod_param_stand.k18:
            tag += "k18" + number_to_str(self.k18 / mod_param_stand.k18)
        if self.k19 != mod_param_stand.k19:
            tag += "k19" + number_to_str(self.k19 / mod_param_stand.k19)
        if self.k20 != mod_param_stand.k20:
            tag += "k20" + number_to_str(self.k20 / mod_param_stand.k20)
        if self.k21 != mod_param_stand.k21:
            tag += "k21" + number_to_str(self.k21 / mod_param_stand.k21)
        if self.k22 != mod_param_stand.k22:
            tag += "k22" + number_to_str(self.k22 / mod_param_stand.k22)
        if self.k23 != mod_param_stand.k23:
            tag += "k23" + number_to_str(self.k23 / mod_param_stand.k23)
        if self.Kt != mod_param_stand.Kt:
            tag += "Kt" + number_to_str(self.Kt / mod_param_stand.Kt)
        if self.Dm != mod_param_stand.Dm:
            tag += "Dm" + number_to_str(self.Dm / mod_param_stand.Dm)
        if self.Dms != mod_param_stand.Dms:
            tag += "Dms" + number_to_str(self.Dms / mod_param_stand.Dms)
        if self.Dmss != mod_param_stand.Dmss:
            tag += "Dmss" + number_to_str(self.Dmss / mod_param_stand.Dmss)
        if self.Dc != mod_param_stand.Dc:
            tag += "Dc" + number_to_str(self.Dc / mod_param_stand.Dc)
        if self.Dm_gdp != mod_param_stand.Dm_gdp:
            tag += "DmGdp" + number_to_str(self.Dm_gdp / mod_param_stand.Dm_gdp)
        if self.r != mod_param_stand.r:
            tag += "r" + number_to_str(self.r / mod_param_stand.r)
        if self.Cdc42_tot != mod_param_stand.Cdc42_tot:
            tag += "Cdc42Tot" + number_to_str(self.Cdc42_tot / mod_param_stand.Cdc42_tot)
        if self.BemGEF_tot != mod_param_stand.BemGEF_tot:
            tag += "BemGefTot" + number_to_str(self.BemGEF_tot / mod_param_stand.BemGEF_tot)
        if self.Axl2_tot != mod_param_stand.Axl2_tot:
            tag += "Axl2_tot" + number_to_str(self.Axl2_tot / mod_param_stand.Axl2_tot)
        if self.exo_axl2 != mod_param_stand.exo_axl2:
            tag += "Exo" + number_to_str(self.exo_axl2 / mod_param_stand.exo_axl2)
        if self.tol_k23 != mod_param_stand.tol_k23:
            tag += "Tol23" + number_to_str(self.tol_k23 / mod_param_stand.tol_k23)
        if self.gamma_d != mod_param_stand.gamma_d:
            tag += "Gd" + number_to_str(self.gamma_d / mod_param_stand.gamma_d)
        if self.gamma_k != mod_param_stand.gamma_k:
            tag += "Gk" + number_to_str(self.gamma_k / mod_param_stand.gamma_k)
        if self.gamma_push != mod_param_stand.gamma_push:
            tag += "Gpush" + number_to_str(self.gamma_push / mod_param_stand.gamma_push)

        if self.k17_alt == True:
            tag += "k17_alt"

        if self.crowding_p == False:
            tag += "_noCrowding"            

        if exocyte_data != None:
            tag += "sigma2" + number_to_str(exocyte_data.sigma2) 
            tag += "Kg" + number_to_str(exocyte_data.kG)
            tag += "Lhit" + number_to_str(exocyte_data.lambda_hit)
            tag += "Lcab" + number_to_str(exocyte_data.lambda_cable)
            tag += "dnew" + number_to_str(exocyte_data.factor_dnew)

        self.tag = tag
        return tag 

class ModelParametersAxl2V4:
    def __init__(self, k1a=10, k1b=10, k2a=0.16, k2b=0.35, k3=0.35, k4a=10, k4b=10, k5a=36,
                 k5b=0.65, k7=10, k12a=0.05, k12b=0.01, k13=0.024, k14=0.024, 
                 k15=1.0, k16=0.05, k17=0.05, k18=0.1, k19=3.0, k20=1.0, k21=0.02, k22=0.5, k23=5.0,
                 k24=5.0, k25=1.0, DmFOO=0.0045, GapSlim=False,
                 Dm=0.01, Dms=0.0025, Dmss=0.00045, Dm_gdp=0.01, Dc=10, eta=0.01, r = 2.5, L=1.0, 
                 BemGEF_tot=0.017, Cdc42_tot=1.0, S_tot = 0.2, GapS_tot=0.005, Kd=150, Axl2_tot=0.2,
                 crowding=False, Kt=1500, exo_axl2=100, recruited=False, use_k23_exp=False, tol_k23=1e-1, 
                 P_pos=False, k17_alt=False, use_gamma=False, gamma_d=0.1, gamma_k=0.3, gamma_push=0.1, 
                 read_run=None, crowding_p=False):

        self.read_run = read_run

        self.crowding_p = crowding_p

        self.k17_alt = k17_alt
        self.P_pos = P_pos
        self.Kd = Kd
        self.crowding=crowding
        self.Kt = Kt
        self.exo_axl2 = exo_axl2
        self.recruited = recruited
        self.use_k23_exp = use_k23_exp
        self.tol_k23 = tol_k23
        self.GapSlim = GapSlim

        self.gamma_d = gamma_d
        self.gamma_k = gamma_k
        self.gamma_push = gamma_push
        self.use_gamma = use_gamma

        self.k1a = k1a
        self.k1b = k1b
        self.k2a = k2a
        self.k2b = k2b
        self.k3 = k3
        self.k4a = k4a
        self.k4b = k4b
        self.k5a = k5a
        self.k5b = k5b
        self.k7 = k7
        self.k12a = k12a
        self.k12b = k12b
        self.k13 = k13
        self.k14 = k14
        self.k15 = k15
        self.k16 = k16
        self.k17 = k17
        self.k18 = k18
        self.k19 = k19
        self.k20 = k20
        self.k21 = k21
        self.k22 = k22
        self.k23 = k23
        self.k24 = k24
        self.k25 = k25

        self.Dm = Dm / L**2
        self.DmFOO = DmFOO / L**2
        self.Dms = Dms / L**2
        self.Dmss = Dmss / L**2
        self.Dm_gdp = Dm_gdp / L**2
        self.Dc = Dc / L**2
        self.L = L
        dx = 10e-9
        rad_m = r * 1e-6
        self.eta = ((rad_m + dx)**3 - rad_m**3) / rad_m**3
        self.r = r
        
        self.ic_Cdc42T = 0.0
        self.ic_Cdc42D = 0.0
        self.ic_BemGEF42 = 0.0
        self.ic_BemGEFm = 0.0
        self.ic_BemGEFc = 0.0
        self.ic_Cdc42I = 0.0
        self.ic_Sc = 0.0
        self.ic_S = 0.0
        self.ic_P = 0.0
        self.ic_GapS = 0.0
        self.ic_GapSc = 0.0
        self.ic_Axl2c = 0.0
        self.ic_Axl2 = 0.0
        self.ic_Axl2_S = 0.0
        self.ic_FOO = 0.0
        
        self.Cdc42_tot = Cdc42_tot
        self.BemGEF_tot = BemGEF_tot
        self.S_tot = S_tot
        self.GapS_tot = GapS_tot
        self.Axl2_tot = Axl2_tot

        # Compensate and ensure Cdc42c is reduced
        self.tag = ""
    
    def make_param_constants(self):
        k1a = Constant(self.k1a)
        k1b = Constant(self.k1b)
        k2a = Constant(self.k2a)
        k2b = Constant(self.k2b)
        k3 = Constant(self.k3)
        k4a = Constant(self.k4a)
        k4b = Constant(self.k4b)
        k5a = Constant(self.k5a)
        k5b = Constant(self.k5b)
        k7 = Constant(self.k7) 
        k12a = Constant(self.k12a)
        k12b = Constant(self.k12b)
        k13 = Constant(self.k13)
        k14 = Constant(self.k14)
        k15 = Constant(self.k15)
        k16 = Constant(self.k16) 
        k17 = Constant(self.k17) 
        k18 = Constant(self.k18) 
        k19 = Constant(self.k19) 
        k20 = Constant(self.k20) 
        k21 = Constant(self.k21) 
        k22 = Constant(self.k22)
        k24 = Constant(self.k24) 
        k25 = Constant(self.k25)         
        Dm = Constant(self.Dm) 
        Dms = Constant(self.Dms) 
        Dmss = Constant(self.Dmss) 
        Dc = Constant(self.Dc) 
        eta = Constant(self.eta)
        Dm_gdp = Constant(self.Dm_gdp)
        
        return k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, k12a, k12b, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k24, k25, Dm, Dms, Dmss, Dc, eta, Dm_gdp
    
    def return_param(self):
        return self.k1a, self.k1b, self.k2a, self.k2b, self.k3, self.k4a, self.k4b, self.k5a, self.k5b, self.k7, self.k12a, self.k12b, self.k13, self.k14, self.k15, self.k16, self.k17, self.k18, self.k19, self.k20, self.k21, self.k22, self.k24, self.k25, self.Dm, self.Dms, self.Dmss, self.Dc, self.eta, self.Dm_gdp
    
    def calc_save_tag(self, bud_scar_param=None, exocyte_data=None):
        
        def number_to_str(number):
            return str(round(number, 3)).replace(".", "P") 
        
        mod_param_stand = ModelParametersAxl2V4()
        
        tag = "_"
        if self.k1a != mod_param_stand.k1a:
            tag += "k1a" + number_to_str(self.k1a / mod_param_stand.k1a)
        if self.k1b != mod_param_stand.k1b:
            tag += "k1b" + number_to_str(self.k1b / mod_param_stand.k1b)
        if self.k2a != mod_param_stand.k2a:
            tag += "k2a" + number_to_str(self.k2a / mod_param_stand.k2a)
        if self.k2b != mod_param_stand.k2b:
            tag += "k2b" + number_to_str(self.k2b / mod_param_stand.k2b)
        if self.k3 != mod_param_stand.k3:
            tag += "k3" + number_to_str(self.k3 / mod_param_stand.k3)
        if self.k4a != mod_param_stand.k4a:
            tag += "k4a" + number_to_str(self.k4a / mod_param_stand.k4a)
        if self.k4b != mod_param_stand.k4b:
            tag += "k4b" + number_to_str(self.k4b / mod_param_stand.k4b)
        if self.k5a != mod_param_stand.k5a:
            tag += "k5a" + number_to_str(self.k5a / mod_param_stand.k5a)
        if self.k5b != mod_param_stand.k5b:
            tag += "k5b" + number_to_str(self.k5b / mod_param_stand.k5b)
        if self.k7 != mod_param_stand.k7:
            tag += "k7" + number_to_str(self.k7 / mod_param_stand.k7)
        if self.k12a != mod_param_stand.k12a:
            tag += "k12a" + number_to_str(self.k12a / mod_param_stand.k12a)
        if self.k12b != mod_param_stand.k12b:
            tag += "k12b" + number_to_str(self.k12b / mod_param_stand.k12b)
        if self.k13 != mod_param_stand.k13:
            tag += "k13" + number_to_str(self.k13 / mod_param_stand.k13)
        if self.k14 != mod_param_stand.k14:
            tag += "k14" + number_to_str(self.k14 / mod_param_stand.k14)
        if self.k15 != mod_param_stand.k15:
            tag += "k15" + number_to_str(self.k15 / mod_param_stand.k15)
        if self.k16 != mod_param_stand.k16:
            tag += "k16" + number_to_str(self.k16 / mod_param_stand.k16)
        if self.k17 != mod_param_stand.k17:
            tag += "k17" + number_to_str(self.k17 / mod_param_stand.k17)
        if self.k18 != mod_param_stand.k18:
            tag += "k18" + number_to_str(self.k18 / mod_param_stand.k18)
        if self.k19 != mod_param_stand.k19:
            tag += "k19" + number_to_str(self.k19 / mod_param_stand.k19)
        if self.k20 != mod_param_stand.k20:
            tag += "k20" + number_to_str(self.k20 / mod_param_stand.k20)
        if self.k21 != mod_param_stand.k21:
            tag += "k21" + number_to_str(self.k21 / mod_param_stand.k21)
        if self.k22 != mod_param_stand.k22:
            tag += "k22" + number_to_str(self.k22 / mod_param_stand.k22)
        if self.k23 != mod_param_stand.k23:
            tag += "k23" + number_to_str(self.k23 / mod_param_stand.k23)
        if self.k24 != mod_param_stand.k24:
            tag += "k24" + number_to_str(self.k24 / mod_param_stand.k24)
        if self.k25 != mod_param_stand.k25:
            tag += "k25" + number_to_str(self.k25 / mod_param_stand.k25)            
        if self.Kt != mod_param_stand.Kt:
            tag += "Kt" + number_to_str(self.Kt / mod_param_stand.Kt)
        if self.Dm != mod_param_stand.Dm:
            tag += "Dm" + number_to_str(self.Dm / mod_param_stand.Dm)
        if self.DmFOO != mod_param_stand.DmFOO:
            tag += "DmFOO" + number_to_str(self.DmFOO / mod_param_stand.DmFOO)        
        if self.Dms != mod_param_stand.Dms:
            tag += "Dms" + number_to_str(self.Dms / mod_param_stand.Dms)
        if self.Dmss != mod_param_stand.Dmss:
            tag += "Dmss" + number_to_str(self.Dmss / mod_param_stand.Dmss)
        if self.Dc != mod_param_stand.Dc:
            tag += "Dc" + number_to_str(self.Dc / mod_param_stand.Dc)
        if self.Dm_gdp != mod_param_stand.Dm_gdp:
            tag += "DmGdp" + number_to_str(self.Dm_gdp / mod_param_stand.Dm_gdp)
        if self.r != mod_param_stand.r:
            tag += "r" + number_to_str(self.r / mod_param_stand.r)
        if self.Cdc42_tot != mod_param_stand.Cdc42_tot:
            tag += "Cdc42Tot" + number_to_str(self.Cdc42_tot / mod_param_stand.Cdc42_tot)
        if self.BemGEF_tot != mod_param_stand.BemGEF_tot:
            tag += "BemGefTot" + number_to_str(self.BemGEF_tot / mod_param_stand.BemGEF_tot)
        if self.Axl2_tot != mod_param_stand.Axl2_tot:
            tag += "Axl2_tot" + number_to_str(self.Axl2_tot / mod_param_stand.Axl2_tot)
        if self.S_tot != mod_param_stand.S_tot:
            tag += "Stot" + number_to_str(self.S_tot / mod_param_stand.S_tot)            
        if self.exo_axl2 != mod_param_stand.exo_axl2:
            tag += "Exo" + number_to_str(self.exo_axl2 / mod_param_stand.exo_axl2)
        if self.tol_k23 != mod_param_stand.tol_k23:
            tag += "Tol23" + number_to_str(self.tol_k23 / mod_param_stand.tol_k23)
        if self.gamma_d != mod_param_stand.gamma_d:
            tag += "Gd" + number_to_str(self.gamma_d / mod_param_stand.gamma_d)
        if self.gamma_k != mod_param_stand.gamma_k:
            tag += "Gk" + number_to_str(self.gamma_k / mod_param_stand.gamma_k)
        if self.gamma_push != mod_param_stand.gamma_push:
            tag += "Gpush" + number_to_str(self.gamma_push / mod_param_stand.gamma_push)
        if self.GapSlim != False:
            tag += "GapSlim" + number_to_str(self.GapSlim)            
            
        if self.k17_alt == True:
            tag += "k17_alt"

        if self.crowding_p == False:
            tag += "_noCrowding"            

        if exocyte_data != None:
            tag += "sigma2" + number_to_str(exocyte_data.sigma2) 
            tag += "Kg" + number_to_str(exocyte_data.kG)
            tag += "Lhit" + number_to_str(exocyte_data.lambda_hit)
            tag += "Lcab" + number_to_str(exocyte_data.lambda_cable)
            tag += "dnew" + number_to_str(exocyte_data.factor_dnew)

        self.tag = tag
        return tag 
    
# Negative feedback model
class ModelParametersNf:
    def __init__(self, k1a=10, k1b=10, k2a=0.16, k2b=0.35, k3=0.35, k4a=10, k4b=10, k5a=36,
                 k5b=0.65, k7=10, k8max=0.0063, k8max2=0.0063, k8n=6, k8h=10, k8h2=250.0, k9max=0.0044, k9n=6, k9h=0.003,
                 Dm=0.011, Dm_gdp=0.023, Dc=10, eta=0.01, r = 2.5, L=1.0, BemGEF_tot=0.017, Cdc42_tot=1.0, 
                 strength_feedback=1.0, second_feedback=False):
        
        self.k1a = k1a
        self.k1b = k1b
        self.k2a = k2a
        self.k2b = k2b
        self.k3 = k3
        self.k4a = k4a
        self.k4b = k4b
        self.k5a = k5a
        self.k5b = k5b
        self.k7 = k7
        self.Dm = Dm / L**2
        self.Dc = Dc / L**2
        self.Dm_gdp = Dm_gdp / L**2
        self.L = L
        dx = 10e-9
        rad_m = r * 1e-6
        self.eta = ((rad_m + dx)**3 - rad_m**3) / rad_m**3
        self.r = r

        self.k8max = k8max
        self.k8h = k8h
        self.k8n = k8n
        self.k8h2 = k8h2
        self.k8max2 = k8max2
        self.second_feedback = second_feedback

        self.k9max = k9max
        self.k9h = k9h
        self.k9n = k9n
        
        self.ic_Cdc42T = 0.0
        self.ic_Cdc42D = 0.0
        self.ic_BemGEF42 = 0.0
        self.ic_BemGEFm = 0.0
        self.ic_BemGEFc = 0.0
        self.ic_Cdc42I = 0.0
        self.ic_BemGEF42_star = 0.0
        self.ic_BemGEFm_star = 0.0
        self.ic_BemGEFc_star = 0.0
        
        self.Cdc42_tot = Cdc42_tot
        self.BemGEF_tot = BemGEF_tot

        # Take the strength of the negative feedback into account 
        self.strength_feedback = strength_feedback
        
        # Compensate and ensure Cdc42c is reduced
        self.tag = ""
    
    
    def make_param_constants(self):
        k1a = Constant(self.k1a)
        k1b = Constant(self.k1b)
        k2a = Constant(self.k2a)
        k2b = Constant(self.k2b)
        k3 = Constant(self.k3)
        k4a = Constant(self.k4a)
        k4b = Constant(self.k4b)
        k5a = Constant(self.k5a)
        k5b = Constant(self.k5b)
        k7 = Constant(self.k7) 
        Dm = Constant(self.Dm) 
        Dc = Constant(self.Dc) 
        Dm_gdp = Constant(self.Dm_gdp)
        eta = Constant(self.eta)

        k8max = Constant(self.k8max )
        k8h = Constant(self.k8h)
        k8n = Constant(self.k8n)

        k9max = Constant(self.k9max)
        k9h = Constant(self.k9h)
        k9n = Constant(self.k9n)
        
        return k1a, k1b, k2a, k2b, k3, k4a, k4b, k5a, k5b, k7, Dm, Dc, eta, k8max, k8h, k8n, k9max, k9h, k9n, Dm_gdp
    
    def return_param(self):
        return (self.k1a, self.k1b, self.k2a, self.k2b, self.k3, self.k4a, self.k4b, self.k5a, self.k5b, self.k7, 
            self.Dm, self.Dc, self.eta, self.k8max, self.k8h, self.k8n, self.k9max, self.k9h, self.k9n, self.Dm_gdp)


    def print_ic(self):
        
        print("Cdc42T = {:.3e}".format(self.ic_Cdc42T))
        print("Cdc42D = {:.3e}".format(self.ic_Cdc42D))
        print("BemGEF42 = {:.3e}".format(self.ic_BemGEF42))
        print("BemGEFm = {:.3e}".format(self.ic_BemGEFm))
        print("BemGEF42_star = {:.3e}".format(self.ic_BemGEF42_star))
        print("BemGEFm_star = {:.3e}".format(self.ic_BemGEFm_star))
        print("Cdc42c = {:.3e}".format(self.ic_Cdc42I))
        print("BemGEFc = {:.3e}".format(self.ic_BemGEFc))
        print("BemGEFc_star = {:.3e}".format(self.ic_BemGEFc_star))


    def calc_save_tag(self, bud_scar_param=None):
        
        def number_to_str(number):
            return str(round(number, 3)).replace(".", "P") 
        
        mod_param_stand = ModelParametersNf()
        
        tag = "_"
        if self.k1a != mod_param_stand.k1a:
            tag += "k1a" + number_to_str(self.k1a / mod_param_stand.k1a)
        if self.k1b != mod_param_stand.k1b:
            tag += "k1b" + number_to_str(self.k1b / mod_param_stand.k1b)
        if self.k2a != mod_param_stand.k2a:
            tag += "k2a" + number_to_str(self.k2a / mod_param_stand.k2a)
        if self.k2b != mod_param_stand.k2b:
            tag += "k2b" + number_to_str(self.k2b / mod_param_stand.k2b)
        if self.k3 != mod_param_stand.k3:
            tag += "k3" + number_to_str(self.k3 / mod_param_stand.k3)
        if self.k4a != mod_param_stand.k4a:
            tag += "k4a" + number_to_str(self.k4a / mod_param_stand.k4a)
        if self.k4b != mod_param_stand.k4b:
            tag += "k4b" + number_to_str(self.k4b / mod_param_stand.k4b)
        if self.k5a != mod_param_stand.k5a:
            tag += "k5a" + number_to_str(self.k5a / mod_param_stand.k5a)
        if self.k5b != mod_param_stand.k5b:
            tag += "k5b" + number_to_str(self.k5b / mod_param_stand.k5b)
        if self.k7 != mod_param_stand.k7:
            tag += "k7" + number_to_str(self.k7 / mod_param_stand.k7)
        if self.k8h != mod_param_stand.k8h:
            tag += "k8h" + number_to_str(self.k8h / mod_param_stand.k8h)
        if self.k8h2 != mod_param_stand.k8h2:
            tag += "k8h2" + number_to_str(self.k8h2 / mod_param_stand.k8h2)            
        if self.k9h != mod_param_stand.k9h:
            tag += "k9h" + number_to_str(self.k9h / mod_param_stand.k9h)            
        if self.k8max != mod_param_stand.k8max:
            tag += "k8max" + number_to_str(self.k8max / mod_param_stand.k8max)
        if self.k8max2 != mod_param_stand.k8max2:
            tag += "k8max2" + number_to_str(self.k8max2 / mod_param_stand.k8max2)            
        if self.k9max != mod_param_stand.k9max:
            tag += "k9max" + number_to_str(self.k9max / mod_param_stand.k9max)                                    
        if self.Dm != mod_param_stand.Dm:
            tag += "Dm" + number_to_str(self.Dm / mod_param_stand.Dm)
        if self.Dc != mod_param_stand.Dc:
            tag += "Dc" + number_to_str(self.Dc / mod_param_stand.Dc)
        if self.r != mod_param_stand.r:
            tag += "r" + number_to_str(self.r / mod_param_stand.r)
        if self.Cdc42_tot != mod_param_stand.Cdc42_tot:
            tag += "Cdc42Tot" + number_to_str(self.Cdc42_tot / mod_param_stand.Cdc42_tot)
        if self.BemGEF_tot != mod_param_stand.BemGEF_tot:
            tag += "BemGefTot" + number_to_str(self.BemGEF_tot / mod_param_stand.BemGEF_tot)
        if self.Dm_gdp != mod_param_stand.Dm_gdp:
            tag += "DmGdp" + number_to_str(self.Dm_gdp / mod_param_stand.Dm_gdp)
        if self.strength_feedback != mod_param_stand.strength_feedback:
            tag += "SF" + number_to_str(self.strength_feedback / mod_param_stand.strength_feedback)

        if bud_scar_param != None:
            tag += "ringk1a" + number_to_str(bud_scar_param.param_ring.k1a / mod_param_stand.k1a)
            tag += "scark2b" + number_to_str(bud_scar_param.param_bud_scar.k2b / mod_param_stand.k2b)
        
        self.tag = tag
        return tag 
    
def calc_pole_ratios(cdc42):
    cdc42_mem = cdc42
         
    # To calculate size of the pole extract information of Cdc42-values
    cdc42_vec = cdc42_mem.vector()
    cdc42_min = cdc42_vec.min()
    cdc42_max = cdc42_vec.max()
    cdc42_len = len(cdc42_vec)

    # Tolerances for determining poles making out a pole
    min_TOL = 0.45 * (cdc42_max - cdc42_min)
    max_TOL = 0.45 * (cdc42_max - cdc42_min)
    max_TOL_strict = 0.2095 * (cdc42_max - cdc42_min)

    # Calculate number of nodes making out a pole
    num_max_dofs, num_min_dofs, num_max_dofs_strict = 0, 0, 0
    for index_cdc42 in range(0, cdc42_len):
        if abs(cdc42_vec[index_cdc42] - cdc42_max) < max_TOL:
            num_max_dofs += 1
        elif abs(cdc42_vec[index_cdc42] - cdc42_min) < min_TOL:
            num_min_dofs += 1            
        
        # Hard max limit (size of most intense part)
        if abs(cdc42_vec[index_cdc42] - cdc42_max) < max_TOL_strict:
            num_max_dofs_strict += 1

    pole_ratio = 100 * num_max_dofs / cdc42_len
    pole_ratio_strict = 100 * num_max_dofs_strict / cdc42_len
    
    return pole_ratio, pole_ratio_strict, cdc42_max

# Class for calculating pole data and writing said data to
# disk. 
class PoleData:
    # cdc42 : Current Cdc42 variable 
    # boundary_mesh : Mesh on the cell-surface 
    def __init__(self, cdc42, boundary_mesh, t, n_poles, term_crit="Max_it"):
        
        #cdc42_mem = interpolate(cdc42, boundary_mesh)
        pole_ratio, pole_ratio_strict, cdc42_max = calc_pole_ratios(cdc42)
        
        self.t = t
        self.cdc42_max = cdc42_max
        self.pole_ratio = pole_ratio
        self.pole_ratio_strict = pole_ratio_strict
        self.term_crit = term_crit
        self.n_poles = n_poles
    
    
    def save_data(self, file_name, file_loc):
        pole_data = pd.DataFrame({"t" : [self.t],
                                  "cdc42_max" : [self.cdc42_max],
                                  "pole_ratio" : [self.pole_ratio],
                                  "pole_ratio_strict" : [self.pole_ratio_strict],
                                  "run_index" : [file_loc.run_index],
                                  "term_crit" : [self.term_crit], 
                                  "n_poles" : [self.n_poles]})
        
        # Ensure that data is not written over
        if not os.path.isfile(file_name):
            pole_data.to_csv(file_name)
        else:
            pole_data.to_csv(file_name, header=False, mode='a')
    
    
    def print_f(self):
        print("Pole ratio = {:.3f}".format(self.pole_ratio))
        print("Pole ratio strict = {:.3f}".format(self.pole_ratio_strict))
        print("t = {:.3f}".format(self.t))
        print("Max Cdc42 = {:.3f}".format(self.cdc42_max))


class SolverOption:
    # cdc42 : Current Cdc42 variable 
    # boundary_mesh : Mesh on the cell-surface 
    def __init__(self, term_times_flat=10, term_times_decrease=10, term_max_it=200000, term_inc_cdc42=2, term_check=True,
                 term_eps_clust=0.2, term_min_samp_clust=5, term_tol_flat=1e-3, term_tol_flat_max=1e-3, term_time=False,
                 adaptive_rho=0.9, adaptive_rel_tol=1e-2, adaptive_check=True, print_pwd_it=1000, seed=False, save_train=False, 
                 save_data_t=False, save_data_t_it=100, save_dist=False, save_run=False, read_old_run=False,
                 solver="two_step"):
        # Termination criteria options 
        self.term_times_flat = term_times_flat
        self.term_times_decrease = term_times_decrease
        self.term_max_it = term_max_it
        self.term_check = term_check
        self.term_time = term_time
        self.term_inc_cdc42T = term_inc_cdc42
        # Adaptive step-length 
        self.adaptive_rho = adaptive_rho
        self.adaptive_rel_tol = adaptive_rel_tol
        # Printing options
        self.print_pwd_it = print_pwd_it
        # For reproducing results
        self.seed = seed
        self.save_train = save_train
        # For clustering and calculating number of poles
        self.term_eps_clust = term_eps_clust
        self.term_min_samp_clust = term_min_samp_clust
        self.save_data_t = save_data_t
        self.save_data_t_it = save_data_t_it
        self.solver = solver
        self.adaptive_check = adaptive_check
        self.term_tol_flat= term_tol_flat
        self.term_tol_flat_max = term_tol_flat_max

        self.save_run = save_run
        self.read_old_run = read_old_run

        self.save_dist = save_dist
        

def save_pole_dist_data(data_pole, t, file_loc, file_name, n_poles):

    data_pole["index"] = [file_loc.run_index for j in range(len(data_pole["x"]))]
    data_pole["t"] = [t for j in range(len(data_pole["x"]))]
    data_pole["n_poles"] = [n_poles for j in range(len(data_pole["x"]))]

    if not os.path.isfile(file_name):
        data_pole.to_csv(file_name)
    else:
        data_pole.to_csv(file_name, header=False, mode='a')


def calc_mod_param_list(k2a_change_list, k2b_change_list, r_list, unit_mesh=True, k2b_stand=0.35):
    
    param_standard = ModelParametersDim()
    k2a_standard = param_standard.k2a
    k2b_standard = k2b_stand
    r_standard = param_standard.r
    
    # Calculate values for k2 and k_2
    k2a_val, k2b_val, r_val = [], [], []
    for k2a in k2a_change_list:
        k2a_val.append(k2a * k2a_standard)
    for k2b in k2b_change_list:
        k2b_val.append(k2b * k2b_standard)
    for r in r_list:
        r_val.append(r * r_standard)
    
    # Create the parameters list, order k2, k_2 and r 
    param_pre_list = [k2a_val, k2b_val, r_list]
    param_list = list(itertools.product(*param_pre_list))
    
    # Create the list with parameter objects
    param_list_ret = []
    for param in param_list:
        if unit_mesh == True:
            param_list_ret.append(ModelParametersDim(k2a=param[0], k2b=param[1], r=param[2], L=param[2]))
        else:
            param_list_ret.append(ModelParametersDim(k2a=param[0], k2b=param[1], r=param[2]))
    
    return param_list_ret

def calc_mod_param_list_competition(k2b_list, Dm_list, BemGEFt_list, r,  unit_mesh=True):
    
    # Calculate values for k2 and k_2
    k2b_val, Dm_val ,BemGEFt_val = [], [], []
    for k2b in k2b_list:
        k2b_val.append(k2b)
    for Dm in Dm_list:
        Dm_val.append(Dm)
    for BemGEFt in BemGEFt_list:
        BemGEFt_val.append(BemGEFt)
    
    # Create the parameters list, order k2, k_2 and r 
    param_pre_list = [k2b_val, Dm_val, BemGEFt_val]
    param_list = list(itertools.product(*param_pre_list))
    
    # Create the list with parameter objects
    param_list_ret = []
    for param in param_list:
        if unit_mesh == True:
            param_list_ret.append(ModelParametersDim(k2b=param[0], Dm=param[1], BemGEF_tot=param[2], r=r, L=r))
        else:
            param_list_ret.append(ModelParametersDim(k2b=param[0], Dm=param[1], BemGEF_tot=param[2], r=r))
    
    return param_list_ret


def update_step_length_new(norm_error, dt_low, rel_tol, rho, p=1.0):

    rho = 0.99
    if norm_error == 0.0:
        dt_new = dt_low
    else:
        dt_new = (rel_tol * rho / norm_error)**(1/p) * dt_low
    # By taking the geometrical average we avoid ending up in zigzag 
    # jumping pattern when simulating. 
    #dt_new = 2.0 * dt_low * dt_new / (dt_low + dt_new)

    if norm_error <= rel_tol:
        progress = True 
    else:
        progress = False

    print("norm_error = {:.3e}, dt_new = {:.3e}, rel_tol = {:.3e}, progress = {}".format(norm_error, dt_new, rel_tol, progress)) 
    
    return dt_new, progress


def update_step_length(U_norm, dt_low, rel_tol, rho, p=1.0):
    
    norm_error = norm(U_norm)
    if norm_error == 0.0:
        dt_new = dt_low
    else:
        dt_new = (rel_tol * rho / norm_error)**(1/p) * dt_low
    # By taking the geometrical average we avoid ending up in zigzag 
    # jumping pattern when simulating. 
    dt_new = 2.0 * dt_low * dt_new / (dt_low + dt_new)

    if norm_error < rel_tol:
        progress = True 
    else:
        progress = False

    if dt_new > 0.1:
        dt_new = 0.1
    
    return dt_new, dt_new / 2.0, progress


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
def calc_n_poles(data, u_max_filt, eps=0.2, min_samples=5):
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
def calc_n_poles_her(data, eps=0.15, min_samples=5):
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


# Functions for reading and saving simultion results. Note, this function 
# assumes that the mesh does not change. 
def save_object(obj, filename):
    # Overwrites any existing file.
    with open(filename, 'wb') as outp:  
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def save_sim_result(U, exp_list_vector, dir_save):

    if not os.path.isdir(dir_save):
        os.makedirs(dir_save)

    U_vector = np.zeros(len(U.vector()))
    j = 0
    for val in U.vector():
        U_vector[j] = val
        j += 1
    save_object(U_vector, dir_save + "/U_val")
    save_object(exp_list_vector, dir_save + "/Exp_list")

def read_sim_result(dir_save):
    with open(dir_save + "U_val", "rb") as f:
        U_val = pickle.load(f)
    with open(dir_save + "Exp_list", "rb") as f:
        exp_val = pickle.load(f)

    return U_val, exp_val


def create_param_bud_scar(increase_import, increase_gap, mod_param, i_bud_scar, i_ring):

    k2b_stand = mod_param.k2b
    k1a_stand = mod_param.k1a

    param_bud_scar = copy.deepcopy(mod_param)
    param_ring = copy.deepcopy(mod_param)
    param_bud_scar.k2b = k2b_stand * increase_gap
    param_ring.k1a = k1a_stand * increase_import

    bud_scar_val = BudScarParam(param_ring, param_bud_scar, i_bud_scar, i_ring)

    return bud_scar_val
