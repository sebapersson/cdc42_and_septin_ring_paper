'''
This file contains functions for converting gmsh msh-files into 
xdmf-formats, and then reading these xdmf-formats info Fenics meshes. 
The conversion from msh to xdmf (required by Fenics) is performed using 
meshio. 

More specifically, the conversion to xdmf-files performed by
the process_msh_file_to_xdmf-function creates three files 
    msh.xdmf: The mesh-file
    subdomains.xdmf: File containing the subdomains 
    surfaces.xdmf: File containing the surfaces (facets)
These files are processed into fenics meshes by read_mesh-function. 
'''


import os
import meshio
from dolfin import *


# Function that will check if a certain mesh file exists, given the number
# of bud-scars.
# Args:
#     n_bud_scars, the number of bud-scars
#     nucleus, bool if the mesh has nucleus
# Returns:
#     true or false depending if file exists 
def msh_file_exists(n_bud_scars, dens_val, nucleus=True):
    # Add if there is a nucleus or not
    if nucleus:
        tag_nuc = "tn_"
    else:
        tag_nuc = "fn_"
    tag_bud_scars = 'b' + n_bud_scars
    
    if dens_val == False:
        msh_file = "./Generate_mesh/Mesh_" + tag_nuc + tag_bud_scars + ".msh"    
        if not os.path.isfile(msh_file):
            if nucleus:
                print("msh-file with nucleus and {} scars does not exist".format(n_bud_scars))
            else:
                print("msh-file without nucleus and {} scars does not exist".format(n_bud_scars))
                return False
        else:
            return True
    elif dens_val != False:
        tag_dens = "_" + str(dens_val).replace(".", "P")
        msh_file = "./Generate_mesh/Find_mesh/Mesh_" + tag_nuc + tag_bud_scars + tag_dens + ".msh"
        if not os.path.isfile(msh_file):
            print("Error: file = {} does not exist when testing mesh size".format(msh_file))
            return False
        return True 


# Function that will read a gmsh msh-file and process it into
# xdmf-format (readable for fenics). Note, output is written
# to file. Conversion is performed using meshio (correct meshio
# version is important). 
# Args:
#     file_locations, file-locations object for an experiment 
# Outputs:
#     msh.xdmf: The mesh-file
#     subdomains.xdmf: File containing the subdomains 
#     surfaces.xdmf: File containing the surfaces (facets)
def process_msh_file_to_xdmf(file_locations):
    
    # For compacter code 
    dir_save_mesh = file_locations.mesh_folder
    path_msh_file = file_locations.msh_file
    
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


# Function that reads the xdmf mesh files and returns
# the mesh, subdomains and surfaces.
# Args:
#     file_locations, file locations for experiment
# Returns:
#     mesh, the mesh object
#     subdomains, the subdomains object
#     surfaces, the surfaces object 
def read_mesh(file_loc):
    
    mesh_folder = file_loc.mesh_folder
    
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
