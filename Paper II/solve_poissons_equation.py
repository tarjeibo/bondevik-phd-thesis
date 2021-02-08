<<<<<<< HEAD
import numpy as np
from math import sqrt
from matplotlib import pyplot as plt
from scipy.integrate import solve_bvp



class PoissonEquationBVPsolver(object):
    def __init__(self, epsilon_r, c_v_bulk, c_oh_bulk, c_acc, T, unit_cell_volume, phi_core, core_charge, phi_core_iterator = 0.001, scl_length = 5e-9):
        self.epsilon_r = epsilon_r
        self.c_v_bulk = c_v_bulk
        self.c_oh_bulk = c_oh_bulk
        self.c_acc = c_acc
        self.T = T
        self.unit_cell_volume = unit_cell_volume
        self.phi_core = phi_core
        self.core_charge = core_charge
        self.phi_core_iterator = phi_core_iterator
        self.scl_length = scl_length

        self.set_physical_constants()
        self.set_scl_array()
        self.iterate_poisson_until_convergence()


    def set_physical_constants(self):
        self.e = 1.6e-19
        self.kBev = 8.6173325e-5
        self.epsilon_0 = 8.854187817e-12

    def set_scl_array(self):
        self.n_scl_elements = 1000
        self.scl_array = np.linspace(0, self.scl_length, self.n_scl_elements)

    def iterate_poisson_until_convergence(self):
        self.scl_charge = 0 #dummy value
        if self.core_charge > 0:
            while self.core_charge + self.scl_charge > 0:
                self.solve_poisson_equation()
                self.scl_charge = self.get_scl_charge_per_area()
                self.phi_core = self.phi_core + self.phi_core_iterator
        else:
            while self.core_charge + self.scl_charge < 0:
                self.solve_poisson_equation()
                self.scl_charge = self.get_scl_charge_per_area()
                self.phi_core = self.phi_core - self.phi_core_iterator

    def solve_poisson_equation(self):
        y = np.zeros((2, self.scl_array.size))
        y[0][0] = self.phi_core
        y[1][0] = -1
        res = solve_bvp(self.bvp_solver, self.boundary_conditions, self.scl_array, y)
        self.phi_scl_list = res.sol(self.scl_array)[0]

    def boundary_conditions(self, y_start, y_end):
        return np.array([y_start[0] - self.phi_core, y_end[1]])

    def bvp_solver(self, x, y):
        kBT = self.kBev*self.T
        self.c_v_scl = 3*self.c_v_bulk*np.exp(-2*y[0]/kBT)/(3 + self.c_v_bulk*(np.exp(-2*y[0]/kBT) - 1) + self.c_oh_bulk*(np.exp(-y[0]/kBT) - 1))
        self.c_oh_scl = 3*self.c_oh_bulk*np.exp(-y[0]/kBT)/(3 + self.c_v_bulk*(np.exp(-2*y[0]/kBT) - 1) + self.c_oh_bulk*(np.exp(-y[0]/kBT) - 1))
        self.chg_density = self.e*(2*self.c_v_scl + self.c_oh_scl - self.c_acc)/self.unit_cell_volume
        return np.vstack((y[1], -self.chg_density/(self.epsilon_0*self.epsilon_r)))

    def get_scl_charge_per_area(self):
        self.set_scl_array()
        return 2*(self.scl_length/self.n_scl_elements)*np.sum(self.chg_density)


def plot_potential(poissonObj, label=None):
    x = 1e9*np.linspace(0, poissonObj.scl_length, len(poissonObj.phi_scl_list))
    plt.plot(x, poissonObj.phi_scl_list, label=label)
    plt.xlabel("$x$ / nm")
    plt.ylabel("$V_E$ / V")
    plt.legend(loc="best")
=======
import numpy as np
from math import sqrt
from matplotlib import pyplot as plt
from scipy.integrate import solve_bvp



class PoissonEquationBVPsolver(object):
    def __init__(self, epsilon_r, c_v_bulk, c_oh_bulk, c_acc, T, unit_cell_volume, phi_core, core_charge, phi_core_iterator = 0.001, scl_length = 5e-9):
        self.epsilon_r = epsilon_r
        self.c_v_bulk = c_v_bulk
        self.c_oh_bulk = c_oh_bulk
        self.c_acc = c_acc
        self.T = T
        self.unit_cell_volume = unit_cell_volume
        self.phi_core = phi_core
        self.core_charge = core_charge
        self.phi_core_iterator = phi_core_iterator
        self.scl_length = scl_length

        self.set_physical_constants()
        self.set_scl_array()
        self.iterate_poisson_until_convergence()


    def set_physical_constants(self):
        self.e = 1.6e-19
        self.kBev = 8.6173325e-5
        self.epsilon_0 = 8.854187817e-12

    def set_scl_array(self):
        self.n_scl_elements = 1000
        self.scl_array = np.linspace(0, self.scl_length, self.n_scl_elements)

    def iterate_poisson_until_convergence(self):
        self.scl_charge = 0 #dummy value
        if self.core_charge > 0:
            while self.core_charge + self.scl_charge > 0:
                self.solve_poisson_equation()
                self.scl_charge = self.get_scl_charge_per_area()
                self.phi_core = self.phi_core + self.phi_core_iterator
        else:
            while self.core_charge + self.scl_charge < 0:
                self.solve_poisson_equation()
                self.scl_charge = self.get_scl_charge_per_area()
                self.phi_core = self.phi_core - self.phi_core_iterator

    def solve_poisson_equation(self):
        y = np.zeros((2, self.scl_array.size))
        y[0][0] = self.phi_core
        y[1][0] = -1
        res = solve_bvp(self.bvp_solver, self.boundary_conditions, self.scl_array, y)
        self.phi_scl_list = res.sol(self.scl_array)[0]

    def boundary_conditions(self, y_start, y_end):
        return np.array([y_start[0] - self.phi_core, y_end[1]])

    def bvp_solver(self, x, y):
        kBT = self.kBev*self.T
        self.c_v_scl = 3*self.c_v_bulk*np.exp(-2*y[0]/kBT)/(3 + self.c_v_bulk*(np.exp(-2*y[0]/kBT) - 1) + self.c_oh_bulk*(np.exp(-y[0]/kBT) - 1))
        self.c_oh_scl = 3*self.c_oh_bulk*np.exp(-y[0]/kBT)/(3 + self.c_v_bulk*(np.exp(-2*y[0]/kBT) - 1) + self.c_oh_bulk*(np.exp(-y[0]/kBT) - 1))
        self.chg_density = self.e*(2*self.c_v_scl + self.c_oh_scl - self.c_acc)/self.unit_cell_volume
        return np.vstack((y[1], -self.chg_density/(self.epsilon_0*self.epsilon_r)))

    def get_scl_charge_per_area(self):
        self.set_scl_array()
        return 2*(self.scl_length/self.n_scl_elements)*np.sum(self.chg_density)


def plot_potential(poissonObj, label=None):
    x = 1e9*np.linspace(0, poissonObj.scl_length, len(poissonObj.phi_scl_list))
    plt.plot(x, poissonObj.phi_scl_list, label=label)
    plt.xlabel("$x$ / nm")
    plt.ylabel("$V_E$ / V")
    plt.legend(loc="best")
>>>>>>> f118bdd0ca2e9d75f39a78c660627660d3e52366
