from dolfin import *
import copy 
import sys 
import numpy as np
from numpy.lib.arraysetops import unique 

class ExocytosisData:
    # cdc42 : Current Cdc42 variable 
    # boundary_mesh : Mesh on the cell-surface 
    def __init__(self, r_exocyte, n_states, values_exocyte, r_old, L = 1.0, lambda_hit = 50 / 60, kd_par_mult=0.5, 
                 n_cables=0, lambda_cable = 1.0 / 60.0, use_cables=False, centroid_cables=False, plot_cables=False, 
                 cutoff=np.inf, factor_dnew=1.0, endocytosis=False, r_endocyte=0.025, 
                 Kg = 5.0, sigma2=0.06**2, cables_window=False, only_print_eps_hit=False,
                 interpolate_new=True, alpha=0.5, gamma=1.0, gamma_P=10.0, i_P=np.inf, cutoff_p_move=20, 
                 i_S_states=None, eps_hit=0.99):
        
        # If P cannot be dragged given a certain concentration 
        self.cutoff_p_move = cutoff_p_move

        # Vector of values on the exocyte 
        if n_states != len(values_exocyte):
            print("Warning : Exocyte value vector is to long compared with number of states")
            print("Will exit.")
            sys.exit(1)

        self.cables_windows = cables_window
        self.only_print_eps_hit = only_print_eps_hit
        self.interpolate_new = interpolate_new
        self.alpha = alpha
        self.gamma = gamma
        self.gamma_P = gamma_P
        self.i_P = i_P
        self.i_S_states = i_S_states
        self.eps=eps_hit

        # Parameters for cables
        self.kG = Kg
        self.sigma2 = sigma2
        self.factor_dnew = factor_dnew

        self.values_exocyte = values_exocyte
        self.n_states = n_states
        self.n_coords = 0
        self.U_new = 0
        self.kd_par_mult = kd_par_mult

        # Set up actin cables 
        self.use_cables = use_cables
        self.n_cables = n_cables
        self.lambda_cable = lambda_cable
        self.cable_index = np.zeros(self.n_cables)
        self.x_coord_cables = np.zeros(self.n_cables)
        self.y_coord_cables = np.zeros(self.n_cables)
        self.z_coord_cables = np.zeros(self.n_cables)
        self.cable_tau = np.zeros(self.n_cables)
        for i in range(self.n_cables):
           self.cable_tau[i] = self.calc_cable_tau(i)
        # For centroid of cables 
        self.centroid_cables = centroid_cables

        # Calculate new radius after inserting exocyte 
        self.endocytosis = endocytosis
        r_new = np.sqrt(r_old**2 + r_exocyte**2)
        self.r_new = r_new / L
        self.r_old = r_old / L
        r_exocyte /= L

        A_exocyte = r_exocyte**2 * 4 * np.pi
        phi = np.arccos(1 - A_exocyte / (2*np.pi*r_new**2))
        d_new = r_new * phi
        self.d_new = d_new *factor_dnew
        self.d_new_old = d_new
        self.A_exocyte = A_exocyte

        print("d_new = {:.3e}".format(self.d_new))

        self.lambda_hit = lambda_hit
        self.cable_file = None
        self.plot_cables = plot_cables

        # Cutoff limit for excytosis rate 
        self.cutoff = cutoff / L

        # If endocytosis occurs its rate matches exocytosis rate 
        self.lambda_endocyte = lambda_hit * 4.0
        self.r_endocyte = r_endocyte / L
        A_endocyte = self.r_endocyte**2 * np.pi * 4
        phi = np.arccos(1 - A_endocyte / (2*np.pi*self.r_old**2))
        self.d_loss = self.r_old * phi
        self.r_new_endocyte =  np.sqrt(self.r_old**2 - self.r_endocyte**2)
        self.A_endocyte = A_endocyte
        
    
    def calc_coords(self, function_space):
        # Arrays with euclidean (x, y, z) and spherical (phi, theta, r) coordinates 
        coord_eu = function_space.tabulate_dof_coordinates()
        coord_eu = coord_eu[range(0, len(coord_eu), self.n_states)] # Want to avoid duplicating 
        coord_sphere = np.zeros((len(coord_eu), 3))
        for i in range(len(coord_sphere)):
            coord_sphere[i, :] = eu_coord_arr_to_sphere(coord_eu[i, :])

        self.coord_eu = coord_eu
        self.coord_sphere = coord_sphere
        self.n_coord = len(coord_eu)
        self.coord_interpolate = np.zeros((self.n_coord, 3))
        self.coord_interpolate_P = np.zeros((self.n_coord, 3))

        U_new = Function(function_space)
        U_new.set_allow_extrapolation(True)
        self.U_new = U_new


    def get_coord_sphere_hit(self, i_hit):
        return copy.deepcopy(self.coord_sphere[i_hit, :])


    def calc_coord_interpolate(self, i, coord_hit_sphere, exocyte=True, P_vec=None, mod_param=None):

        dist_hit = geodist(sphere_coord_arr_to_eu(coord_hit_sphere), self.coord_eu[i], self.r_old)    
        if self.interpolate_new == False:    
            # Calculating for exocytosis 
            if exocyte == True:
                if dist_hit < self.cutoff:
                    self.coord_interpolate[i, :] = calc_backward_point(coord_hit_sphere, self.coord_sphere[i, :], self.r_new, self.d_new, self.d_new_old)
                else:
                    self.coord_interpolate[i, :] = self.coord_eu[i]
            # For endocytosis 
            else:
                if dist_hit < self.cutoff:
                    coord_hit_sphere[1] += 0.01
                    self.coord_interpolate[i, :] = calc_backward_point_endo(coord_hit_sphere, self.coord_sphere[i, :], self.r_new_endocyte, self.d_loss)
                    coord_hit_sphere[1] -= 0.01
                else:
                    self.coord_interpolate[i, :] = self.coord_eu[i]

            move_eu_coord_in_r(self.coord_interpolate[i, :], self.r_old)
        else:
            coord_hit_eu = sphere_coord_arr_to_eu(coord_hit_sphere)
            if exocyte == True:
                if dist_hit < self.cutoff:
                    self.coord_interpolate[i, :] = find_point_b_exo(coord_hit_eu, self.coord_eu[i, :], self.A_exocyte, 
                        self.r_old, self.alpha, self.gamma)
                    # Make P resistent to pushing (and subsequently S)
                    if P_vec is None:
                        self.coord_interpolate_P[i, :] = find_point_b_exo(coord_hit_eu, self.coord_eu[i, :], self.A_exocyte, 
                            self.r_old, 1.0, self.gamma_P)
                    else:
                        self.coord_interpolate_P[i, :] = find_point_b_exo(coord_hit_eu, self.coord_eu[i, :], self.A_exocyte, 
                            self.r_old, 1.0, self.gamma_P * mod_param.gamma_push * (P_vec[i] + 1.0))
                else:
                    self.coord_interpolate[i, :] = self.coord_eu[i]
            # Endocytosis 
            else:
                if dist_hit < self.cutoff:
                    self.coord_interpolate[i, :] = find_point_c_endo(coord_hit_eu, self.coord_eu[i, :], self.A_endocyte, 
                        self.r_old, self.alpha, self.gamma)
                    self.coord_interpolate_P[i, :] = find_point_c_endo(coord_hit_eu, self.coord_eu[i, :], self.A_endocyte, 
                        self.r_old, 1.0, self.gamma_P)
                else:
                    self.coord_interpolate[i, :] = self.coord_eu[i]

    
    def calc_site_hit(self, Cdc42T):
        
        Cdc42T_vec = Cdc42T.vector().get_local()
        kd_par = np.max(Cdc42T_vec) * self.kd_par_mult 
        prob_vec = Cdc42T_vec**4 / (kd_par**4 + Cdc42T_vec**4)
        prob_vec /= np.sum(prob_vec)

        i_hit = np.random.choice(range(self.n_coord), p=prob_vec)

        return i_hit

    def calc_site_hit_endo(self, Cdc42T):
        Cdc42T_vec = Cdc42T.vector().get_local()
        kd_par = np.max(Cdc42T_vec) * 0.5
        prob_vec = Cdc42T_vec**4 / (kd_par**4 + Cdc42T_vec**4)
        prob_vec /= np.sum(prob_vec)

        i_hit = np.random.choice(range(self.n_coord), p=prob_vec)

        return i_hit

    def calc_site_hit_septin(self, Cdc42T, P):
        
        Cdc42T_vec = Cdc42T.vector().get_local()
        P_vec = P.vector().get_local()
        prob_vec = np.abs(Cdc42T_vec**4 / (1 + P_vec))
        prob_vec /= np.sum(prob_vec)

        i_hit = np.random.choice(range(self.n_coord), p=prob_vec)

        return i_hit

    def calc_cable_tau(self, index):
        return -np.log(np.random.rand()) / self.lambda_cable

    def calc_cable_attach(self, Cdc42T, i_cable, n_cables_use=0, first=False):

        Cdc42T_vec = Cdc42T.vector().get_local()
        if first == True:
            i_max = np.argmax(Cdc42T_vec)
            print(self.coord_eu[i_max])
            print("i_max = {}".format(i_max))
            self.cable_index[i_cable] = i_max
            return 
        
        kd_par = np.max(Cdc42T_vec) * self.kd_par_mult 
        prob_vec = Cdc42T_vec**4 / (kd_par**4 + Cdc42T_vec**4)
        prob_vec /= np.sum(prob_vec)
        
        if self.centroid_cables == True and n_cables_use > 0:
            centroid = self.calc_cable_centroid(n_cables_use)
            geodist_vec = np.array([geodist(centroid, self.coord_eu[i], self.r_old) for i in range(len(self.coord_eu))])
            prob_vec_p2 =  np.exp(-1 * geodist_vec**2 / (2*self.sigma2))
            prob_vec_p2 /= np.sum(prob_vec_p2)
            prob_vec += self.kG * prob_vec_p2
        
        prob_vec /= np.sum(prob_vec)
        i_hit = np.random.choice(range(self.n_coord), p=prob_vec)
        self.cable_index[i_cable] = i_hit

    def calc_cable_attach_window(self, Cdc42T, i_cable):

        Cdc42T_vec = Cdc42T.vector().get_local()
        Cdc42T_max = Cdc42T_vec.max()
        eps = self.eps
        
        prob_vec = Cdc42T_vec > Cdc42T_max*eps
        prob_vec = prob_vec.astype(float)
        prob_vec /= np.sum(prob_vec)
        i_hit = np.random.choice(range(self.n_coord), p=prob_vec)
        self.cable_index[i_cable] = i_hit

        print("Nodes accepted = {}".format(np.sum(Cdc42T_vec > Cdc42T_max*eps)))

    def calc_hit_window(self, Cdc42T):

        Cdc42T_vec = Cdc42T.vector().get_local()
        Cdc42T_max = Cdc42T_vec.max()
        eps = self.eps
        
        prob_vec = Cdc42T_vec > Cdc42T_max*eps
        prob_vec = prob_vec.astype(float)
        prob_vec /= np.sum(prob_vec)
        i_hit = np.random.choice(range(self.n_coord), p=prob_vec)
        print("Nodes accepted = {} / {}".format(np.sum(Cdc42T_vec > Cdc42T_max*eps), len(Cdc42T_vec)))        
        return i_hit
 
    def calc_cable_attach_septin(self, Cdc42T, P, i_cable, n_cables_use=0, first=False):
        
        Cdc42T_vec = Cdc42T.vector().get_local()
        P_vec = P.vector().get_local()
        prob_vec = np.abs(Cdc42T_vec**4 / (1 + P_vec))
        prob_vec /= np.sum(prob_vec)

        if first == True:
            i_max = np.argmax(Cdc42T_vec)
            print(self.coord_eu[i_max])
            print("i_max = {}".format(i_max))
            self.cable_index[i_cable] = i_max
            return 
         
        if self.centroid_cables == True and n_cables_use > 0:
            centroid = self.calc_cable_centroid(n_cables_use)
            geodist_vec = np.array([geodist(centroid, self.coord_eu[i], self.r_old) for i in range(len(self.coord_eu))])
            prob_vec_p2 =  np.exp(-1 * geodist_vec**2 / (2*self.sigma2))
            prob_vec_p2 /= np.sum(prob_vec_p2)
            prob_vec += self.kG * prob_vec_p2
        
        prob_vec /= np.sum(prob_vec)
        i_hit = np.random.choice(range(self.n_coord), p=prob_vec)

        self.cable_index[i_cable] = i_hit

    def calc_cable_centroid(self, n_cables_use):
        # Get cable coordinates 
        cable_coord = np.zeros((3, n_cables_use))
        for i in range(n_cables_use):
            cable_coord[:, i] = self.coord_eu[int(self.cable_index[i])]
            move_eu_coord_in_r(cable_coord[:, i], 1.0)

        if n_cables_use == 1 or len(np.unique(cable_coord, axis=1)) == 3:
            centroid = cable_coord[:, 0]
        else:
            centroid = spherical_centroid(cable_coord)
        move_eu_coord_in_r(centroid, self.r_old)
        return centroid
        

    def update_cable_pos(self, t_curr, Cdc42T):
        for i in range(self.n_cables):
            if self.cable_tau[i] < t_curr:
                self.cable_tau[i] += self.calc_cable_tau(self.lambda_cable)
                if self.cables_windows == False:
                    self.calc_cable_attach(Cdc42T, i, n_cables_use=self.n_cables)
                else:
                    self.calc_cable_attach_window(Cdc42T, i)
                print("Updating position for cable {}, time new cable {}".format(i, self.cable_tau[i]))

    def update_cable_pos_septin(self, t_curr, Cdc42T, P):
        for i in range(self.n_cables):
            if self.cable_tau[i] < t_curr:
                self.cable_tau[i] += self.calc_cable_tau(self.lambda_cable)
                if self.cables_windows == False:
                    self.calc_cable_attach_septin(Cdc42T, P, i, n_cables_use=self.n_cables)
                else:
                    self.calc_cable_attach_window(Cdc42T, i)
                print("Updating position for cable {}, time new cable {}".format(i, self.cable_tau[i]))

    def calc_site_hit_cable(self):
        cable_hit = np.random.randint(low=0, high=self.n_cables)
        i_hit = self.cable_index[cable_hit]
        return int(i_hit)

    def calc_time_hit(self):
        return -np.log(np.random.rand()) / self.lambda_hit

    def calc_time_endo(self):
        return -np.log(np.random.rand()) / self.lambda_endocyte
    
    def set_values_post_exocyte(self, U_old):

        val_coord = np.zeros(self.n_states)

        for i in range(self.n_coord):

            coord = self.coord_interpolate[i, :]
            coord_P = self.coord_interpolate_P[i, :]
            x, y, z = coord
            
            # In case of hit
            hit_cord = x == 0 and y == 0 and z == 0
            if hit_cord:
                val_coord = self.values_exocyte
                val_coord_P = self.values_exocyte
            else:
                val_coord = U_old(coord)
                val_coord_P = U_old(coord_P)
            
            i1 = i * self.n_states
            i2 = i1 + (self.n_states)
            k = 0
            for j in range(i1, i2):
                if self.i_P == (k + 1):
                    # P cannot be dragged if at to high concentration 
                    if val_coord_P[k] > self.cutoff_p_move and not hit_cord:
                        print("Coord = ", coord)
                        print("Coord_P = ", coord_P)
                        print("val_coord_P[k] = ", val_coord_P[k])
                        self.U_new.vector()[j] = U_old.vector()[j]
                    else:
                        self.U_new.vector()[j] = val_coord_P[k] # Polymerised septin can dependent on i_P be allowed to move less 
                # Do not drag S or Axl2_S if P is to high 
                elif self.i_S_states != None and (k + 1) in self.i_S_states: 
                    for i_S in self.i_S_states:
                        if i_S == (k + 1) and not hit_cord:
                            # P cannot be dragged if at to high concentration 
                            if val_coord_P[k] > self.cutoff_p_move:
                                self.U_new.vector()[j] = U_old.vector()[j]
                            else:
                                self.U_new.vector()[j] = val_coord_P[k] # Polymerised septin can dependent on i_P be allowed to move less 
                else:
                    self.U_new.vector()[j] = val_coord[k]
                k += 1

    def set_values_post_endocyte(self, U_old):
        val_coord = np.zeros(self.n_states)

        for i in range(self.n_coord):

            coord = self.coord_interpolate[i, :]
            coord_P = self.coord_interpolate_P[i, :]
            val_coord = U_old(coord)
            val_coord_P = U_old(coord_P)
            i1 = i * self.n_states
            i2 = i1 + (self.n_states)
            k = 0
            for j in range(i1, i2):
                if self.i_P == (k + 1):
                    self.U_new.vector()[j] = val_coord_P[k] # Polymerised septin can dependent on i_P be allowed to move less 
                else:
                    self.U_new.vector()[j] = val_coord[k]
                k += 1

    def plot_location_cables(self, U_old, state_print, file_loc, t_curr):

        if self.cable_file == None:
            self.cable_file = File(file_loc.dir_pvd + "Loc_cables" + ".pvd")

        k = 0
        # Avoid duplications 
        cable_arr = unique(self.cable_index)
        old_val = np.zeros(self.n_cables * self.n_states)
        for i in range(len(cable_arr)):
            i1 = int(cable_arr[i]) * self.n_states
            i2 = i1 + (self.n_states)
            for j in range(i1, i2):
                old_val[k] = U_old.vector()[j]
                U_old.vector()[j] = 1000
                k += 1

        self.cable_file << (state_print, t_curr)

        # Reset old values 
        k = 0        
        for i in range(len(cable_arr)):
            i1 = int(cable_arr[i]) * self.n_states
            i2 = i1 + (self.n_states)
            for j in range(i1, i2):
                U_old.vector()[j] = old_val[k]
                k += 1

# From an array with euclidean coordinates calcualte the corresponding spherical 
# coordinates. 
def eu_coord_arr_to_sphere(coord):

    x, y, z = coord
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    phi = np.arccos(z / r)
    theta = np.arctan2(y, x)

    return np.array([phi, theta, r])


# From spherical coordinates calculate the corresponding euclidean coordinates. 
def sphere_coord_arr_to_eu(coord):
    phi, theta, r = coord
    return np.array([r * np.cos(theta) * np.sin(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(phi)])


# Rescale the vector coord_arr_eu to have length r_new. 
def move_eu_coord_in_r(coord_arr_eu, r_new):
    len_vec = np.sqrt(np.sum(coord_arr_eu**2))

    if len_vec == 0:
        return 

    coord_arr_eu *= r_new / len_vec


# Calculate which coordinate c-old should correspond to if an exocyte hits 
# at the point a-old. r_new is the radius of the sphere after an exocyte 
# has merged with the membrane. Note, a_old and c_old are given in spherical 
# coordinates. Returns, a, b, c, where b is location where c should retreive 
# its concentrations. 
def calc_backward_point(a_old, c_old, r_new, d_new, d_new_old):

    x, y, z = sphere_coord_arr_to_eu(c_old)

    a_old = copy.deepcopy(a_old)
    c_old = copy.deepcopy(c_old)
    # a = [phi, theta, r] (spherical coordinates)
    r_old = a_old[2]

    a_old[2] = r_new
    c_old[2] = r_new
    # Expand a and b radially 
    a = sphere_coord_arr_to_eu(a_old)
    c = sphere_coord_arr_to_eu(c_old)

    if np.abs(np.dot(a, c) / r_new**2) > 1.0:
        arg_arccos = 1.0
    else:
        arg_arccos = np.dot(a, c) / r_new**2

    d_ac = r_new * np.arccos(arg_arccos)
    d_old = d_ac - d_new

    if d_ac < d_new_old:
        return np.array([0.0, 0.0, 0.0])

    epsilon = d_old / r_old
    alpha = np.arccos(np.dot(a, c) / r_new**2)
    phi = alpha - epsilon

    # Right hand side equation system 
    A1 = r_new**2 * np.cos(epsilon) - np.dot(a, c)
    A2 = r_new**2 * np.cos(phi) - np.dot(c, c)
    A3 = r_new**2

    # Extract individual components of a and b for additional calculations 
    a1, a2, a3 = a
    c1, c2, c3 = c

    Dt = (a1*c2-a2*c1)
    Q1 = (A1*c2 - A2*a2) + Dt*c1
    P1 = (a2*c3 - a3*c2) 

    Q2 = (-A1*c1 + A2*a1) + Dt*c2 
    P2 = (-a1*c3 + a3*c1) 

    # Double root. The expression under the square root does not have to 
    # be calculated. Check if denominator is zero for stabillity 
    root = (-Dt**2*c3 - P1*Q1 - P2*Q2) / (Dt**2 + P1**2 + P2**2)
    d3 = root

    nominator = (a1*c2 - a2*c1)
    
    # Case where demnominator is zero (avoid numerical errors)
    if nominator == 0.0:
        d1, d2 = 0.0, 0.0
    else:
        d1 = (A1*c2 - A2*a2 + a2*c3*d3 - a3*c2*d3) / nominator
        d2 = (-A1*c1 + A2*a1 - a1*c3*d3 + a3*c1*d3) / nominator

    # For few cases, numerical errors can propegate 
    if np.abs(d1) > 1e1:
        d1 = 0.0
    if np.abs(d2) > 1e1:
        d2 = 0.0
    
    d = np.array([d1, d2, d3])
    b = c + d

    return b


# Calculate the centroid of points in sphere, when the points in stored in an 
# array where the rows are coordinates, and the columns are points. 
def geodist(x, y, r):
    dotprod = np.dot(x, y) / r**2
    if dotprod > 1.0:
        dotprod = 1.0
    elif dotprod < -1.0:
        dotprod = -1.0
    
    return r * np.arccos(dotprod)
def distance(x, y, axis=0):
    return sqrt((np.power(x - y, 2)).sum(axis=axis))
def improve_centroid(c, ps):
    ans = np.sum([ps[:, i] / np.sqrt(1 - np.power(ps[:, i] @ c, 2)) for i in range(ps.shape[1])], axis=0)
    norm = sqrt(ans @ ans)
    return ans / norm 
def fixpoint(f, x0, eps=1e-5, maxiter=1000, **kwargs):
    for _ in range(maxiter):
        x = f(x0, **kwargs)
        if distance(x, x0) < eps:
            return x
        x0 = x
    raise Exception("Did not converge")
def spherical_centroid(ps, eps=1e-5, maxiter=10000):
    return fixpoint(improve_centroid, np.zeros((3,)), ps=ps, eps=eps, maxiter=maxiter)


# Calculate which coordinate c-old should correspond to if an endocyte hits 
# at the point a-old. r_new is the radius of the sphere after an exocyte 
# has merged with the membrane. Note, a_old and c_old are given in spherical 
# coordinates. Returns, c-new, where c_new is coordinates c_old should retrive 
# its concentration from 
def calc_backward_point_endo(a_old, b_old, r_new, d_loss):

    # Need for further calculations 
    a_old_cp = copy.deepcopy(a_old)
    b_old_cp = copy.deepcopy(b_old)
    r_old = a_old[2]

    # Shrink radially to new radius, a = [phi, theta, r] (spherical coordinates)
    a_old_cp[2] = r_new
    b_old_cp[2] = r_new
    a_small = sphere_coord_arr_to_eu(a_old_cp)
    b_small = sphere_coord_arr_to_eu(b_old_cp)
    a_old = sphere_coord_arr_to_eu(a_old)
    b_old = sphere_coord_arr_to_eu(b_old)

    d_small_cell = geodist(a_small, b_small, r_new)
    d_ac_old = d_small_cell + d_loss
    d_ab = geodist(a_old, b_old, r_old)
    d_bc = d_ac_old - d_ab

    epsilon = d_bc / r_old # Angle between b and c 
    theta = d_ac_old / r_old # Angle between a and c 
    phi = d_ab / r_old
    
    # Right hand side equation system 
    A1 = r_old**2 * np.cos(theta) 
    A2 = r_old**2 * np.cos(epsilon)
    A3 = r_old**2

    a1, a2, a3 = a_old
    b1, b2, b3 = b_old

    Dt = a1*b2 - a2*b1
    Q1 = (A1*b2 - A2*a2) 
    P1 = (a2*b3 - a3*b2) 
    Q2 = -A1*b1 + A2*a1
    P2 = -a1*b3 + a3*b1

    c3 = (-P1*Q1 - P2*Q2) / (Dt**2 + P1**2 + P2**2)    
    if Dt == 0:
        c1 = 0.0
        c2 = 0.0
    else:
        c1 = (A1*b2 - A2*a2 + a2*b3*c3 - a3*b2*c3) / (a1*b2 - a2*b1)
        c2 = (-A1*b1 + A2*a1 - a1*b3*c3 + a3*b1*c3) / (a1*b2 - a2*b1)

    # For few cases, numerical errors can propegate 
    if np.abs(c1) > 1e1:
        c1 = 0.0
    if np.abs(c2) > 1e1:
        c2 = 0.0
    
    c = np.array([c1, c2, c3])

    return c


def slerp(t, p0, p1):
    Omega = np.arccos(np.dot(p0, p1)) / np.dot(p0, p0)    
    return np.sin((1-t)*Omega) / np.sin(Omega) * p0 + np.sin(t*Omega)/ np.sin(Omega)*p1


def calc_dr_exo(p, p0, A_exo, R, alpha=1.0, gamma=1.0):
    R_prim = np.sqrt(R**2 + A_exo/(4*np.pi))
    s = geodist(p, p0, R)

    arg_arccos = 1.0 - (R/R_prim)**2*(1.0 - np.cos(s/R)) - A_exo/(2*np.pi*R_prim**2*alpha)
    if arg_arccos <= -1.0:
        return 0.0
    else:
        dr = 1.0 / gamma * (R_prim*np.arccos(arg_arccos) - s)

    return dr


def calc_dr_endo(point_a, point_b, A_endo, R, alpha=1.0, gamma=1.0):
    R_prim = np.sqrt(R**2 - A_endo/(4*np.pi))
    s = geodist(point_b, point_a, R)

    arg_arccos = 1.0 - (R/R_prim)**2*(1.0 - np.cos(s/R)) + A_endo/(2*np.pi*R_prim**2*alpha)

    if arg_arccos <= -1.0 or arg_arccos == np.inf:
        return 0.0
    elif arg_arccos >= 1.0:
        dr = s / gamma
    else:
        dr = 1.0 / gamma * (-R_prim*np.arccos(arg_arccos) + s)

    return dr


def f_endo(dr, point_a, point_b, point_c, R):
    dist_ab = geodist(point_a, point_b, R)
    dist_ac = geodist(point_a, point_c, R)

    return -((dist_ac - dr) - dist_ab)**2


def f_exo(dr, point_a, point_b, point_c, R):
    dist_ab = geodist(point_a, point_b, R)
    dist_ac = geodist(point_a, point_c, R)

    return -(dist_ac - (dist_ab + dr))**2


def find_point_b_exo(point_a, point_c, A_exo, R, alpha=1.0, gamma=1.0, eps=1e-6):
    a = 1e-18
    b = 1.0 - a

    # Check special cases with the opposite end of the cell 
    dr_c = calc_dr_exo(point_a, point_c, A_exo, R, alpha, gamma)
    if dr_c == 0.0:
        return copy.deepcopy(point_c)

    # Check special case when point_c is inside the insertion zone 
    d_ac = geodist(point_a, point_c, R)
    if d_ac <= dr_c:
        return np.array([0.0, 0.0, 0.0])
    
    k = 1
    tol_list = [1e-10, 1e-3, 5e-2]
    for tol in tol_list:
        k = 1
        while k < 10000:

            t1 = (a + b)*0.5 - eps
            t2 = (a + b)*0.5 + eps

            # Sherical interpolation to get potential points 
            point_b_t1 = slerp(t1, point_a, point_c)  
            point_b_t2 = slerp(t2, point_a, point_c)  

            # Calc cost f 
            dr_t1 = calc_dr_exo(point_a, point_b_t1, A_exo, R, alpha, gamma)
            dr_t2 = calc_dr_exo(point_a, point_b_t2, A_exo, R, alpha, gamma)
            f_t1 = f_exo(dr_t1, point_a, point_b_t1, point_c, R)
            f_t2 = f_exo(dr_t2, point_a, point_b_t2, point_c, R)

            point_b_mean = (point_b_t1 + point_b_t2)*0.5
            dr_mean = calc_dr_exo(point_a, point_b_mean, A_exo, R, alpha, gamma)
            f_mean = f_exo(dr_mean, point_a, point_b_mean, point_c, R)

            if np.abs(f_mean) < tol:
                return (point_b_t1 + point_b_t2) * 0.5
            if f_t1 > f_t2:
                b = t2
            else:
                a = t1
            k += 1
            
    print("Warning : Could not find interpolation point exocytosis")
    print("f_mean = {:.3e}".format(f_mean))


def find_point_c_endo(point_a, point_b, A_endo, R, alpha=1.0, gamma=1.0, eps=1e-6):

    if np.sum(point_a - point_b) == 0.0:
        return copy.deepcopy(point_b)
    
    # Calculate a point beyond point c, and use to find point c
    d_ab = geodist(point_a, point_b, R)
    d_bc = calc_dr_endo(point_a, point_b, A_endo*4, R, alpha, gamma) # A_endo * 4 = A_exo

    angle_ac = (d_ab+d_bc) / R
    angle_bc = d_bc / R
    A1 = R**2 * np.cos(angle_ac)
    A2 = R**2 * np.cos(angle_bc)
    A3 = R**2

    a1, a2, a3 = point_a
    b1, b2, b3 = point_b
    Dt = a1*b2 - a2*b1
    Q1 = (A1*b2 - A2*a2) 
    P1 = (a2*b3 - a3*b2) 
    Q2 = -A1*b1 + A2*a1
    P2 = -a1*b3 + a3*b1
    # If opposite end of the cell 
    if Dt**2 + P1**2 + P2**2 == 0.0:
        return copy.deepcopy(point_b)
    c3 = (-P1*Q1 - P2*Q2) / (Dt**2 + P1**2 + P2**2)    
    if Dt == 0:
        c1 = 0.0
        c2 = 0.0
    else:
        c1 = (A1*b2 - A2*a2 + a2*b3*c3 - a3*b2*c3) / (a1*b2 - a2*b1)
        c2 = (-A1*b1 + A2*a1 - a1*b3*c3 + a3*b1*c3) / (a1*b2 - a2*b1)

    # For few cases, numerical errors can propegate 
    if np.abs(c1) > 1e1:
        c1 = 0.0
    if np.abs(c2) > 1e1:
        c2 = 0.0
    point_c_tmp = np.array([c1, c2, c3])

    # Try to find point c by solving optimisation problem 
    a = 1e-18
    b = 1.0 - a
    d_ab = geodist(point_a, point_b, R)
    d_ac = geodist(point_a, point_c_tmp, R)

    # Check if opposite end of cell 
    dr_c = calc_dr_endo(point_a, point_c_tmp, A_endo, R, alpha=1.0, gamma=1.0)
    if dr_c == 0.0:
        return copy.deepcopy(point_b)

    # Check if inside endocytosis zone 
    d_ac = geodist(point_a, point_c_tmp, R)
    if d_ac <= dr_c:
        return np.array([0.0, 0.0, 0.0])
    
    k = 1
    tol_list = [1e-10, 1e-3, 5e-2]
    for tol in tol_list:
        k = 1
        while k < 100:

            t1 = (a + b)*0.5 - eps
            t2 = (a + b)*0.5 + eps

            # Sherical interpolation to get potential points 
            point_c_t1 = slerp(t1, point_b, point_c_tmp)  
            point_c_t2 = slerp(t2, point_b, point_c_tmp)  

            # Calc cost f 
            dr_t1 = calc_dr_endo(point_a, point_c_t1, A_endo, R, alpha, gamma)
            dr_t2 = calc_dr_endo(point_a, point_c_t2, A_endo, R, alpha, gamma)
            f_t1 = f_endo(dr_t1, point_a, point_b, point_c_t1, R)
            f_t2 = f_endo(dr_t2, point_a, point_b, point_c_t2, R)

            point_c_mean = (point_c_t1 + point_c_t2)*0.5
            dr_mean = calc_dr_endo(point_a, point_c_mean, A_endo, R, alpha, gamma)
            f_mean = f_endo(dr_mean, point_a, point_b, point_c_mean, R)

            if np.abs(f_mean) < tol:
                return (point_c_t1 + point_c_t2) * 0.5
            if f_t1 > f_t2:
                b = t2
            else:
                a = t1

            k += 1

    print("Warning : Could not find interpolation point endocytosis")
    print("f_mean = {:.3e}, tol {:.3e}".format(np.abs(f_mean), tol))
