import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

from cmath import sqrt
#Jacobi solver

#Fundamental constants
epsilon_0 = 8.854187817620e-12
k_B = 1.3806488e-23
k_B_eV = 8.6173325e-5
e = 1.602176e-19
j = sqrt(-1)


#Material specific constants
epsilon_r = 75
a_0 = 4.2e-10

#Experimental params
T = 600
# c_v_bulk = 0.0075043535 #These numbers are from Comparing models article, at T=600K
# c_oh_bulk = 0.0984991292912 #These numbers are from Comparing models article, at T=600K

c_v_bulk = 0.0
c_oh_bulk = 0.1
c_acc = 0.1
resistivity_bulk = 1

#Numerical params
N = 1000
n_iterations = 50
w_start = -4
w_end = 12
w_steps = 5000
w_list = np.logspace(w_start, w_end, w_steps)
convergence_criterion = 1e-3


def solve_poissons_equation(phi, charge_density, n_iterations, region_width):
    N = len(phi)
    epsilon = epsilon_0*epsilon_r
    iteration = 0
    step_length = region_width/len(phi)
    while iteration < n_iterations:
        phi_temp = np.copy(phi)
        phi_temp[0] = 0 #Boundary condition: phi(-zmax)=0
        phi_temp[N-1] = 0 #Boundary condition: phi(zmax)=0
        for i in range (1, N-1): # Avoid updating the boundaries
            phi[i] = 0.5*(phi_temp[i+1] + phi_temp[i-1] + ((step_length**2)/epsilon)*charge_density[i])
        phi[0], phi[N-1] = phi_temp[0], phi_temp[N-1]
        iteration += 1
    return phi

def get_scl_concentrations(phi, e_seg_c_v, e_seg_c_oh, c_v_bulk, c_oh_bulk, T):
    c_v_exp_term = np.exp(-(e_seg_c_v + 2*phi)/(k_B_eV*T))
    c_oh_exp_term = np.exp(-(e_seg_c_oh + phi)/(k_B_eV*T))

    c_v_scl = (3*c_v_bulk*c_v_exp_term)/(3 + c_v_bulk*(c_v_exp_term - 1) + c_oh_bulk*(c_oh_exp_term - 1))
    c_oh_scl = (3*c_oh_bulk*c_oh_exp_term)/(3 + c_v_bulk*(c_v_exp_term - 1) + c_oh_bulk*(c_oh_exp_term - 1))
    scl_concentrations = SclConcentrations(c_v_scl, c_oh_scl)
    return scl_concentrations

def get_charge_density(scl_concentrations, c_acc):
    return (e/a_0**3)*(2*scl_concentrations.c_v_scl + scl_concentrations.c_oh_scl - c_acc)


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def get_gaussian_params(x, y):
    mean = sum(x*y)/len(x)
    sigma = 5e-10
    popt, pcov = curve_fit(gaussian, x, y, p0=[-1, mean, sigma])
    return popt


def get_scl_impedance_params(phi, w_list):
    resistivity_scl = resistivity_bulk*np.exp(phi/(k_B_eV*T))
    step_length = region_width/len(phi)
    z_scl = []
    for w in w_list:
        z_scl_at_w = np.sum((resistivity_scl*step_length)/(1 + j*w*resistivity_scl*epsilon_0*epsilon_r))
        z_scl.append(z_scl_at_w)
    z_scl = np.asarray(z_scl)
    scl_impedance_params = SclImpendanceParams(resistivity_scl, z_scl)
    return scl_impedance_params

def get_char_freq_ratio(z_scl):
    w0_scl_idx = np.argmax(-z_scl.imag)
    w0_scl = w_list[w0_scl_idx]
    w0_bulk = 1/(resistivity_bulk*epsilon_0*epsilon_r)
    return w0_scl/w0_bulk

def get_tot_chg_e_units(charge_density):
    step_length = region_width/len(phi)
    gb_area = sqrt(5)*a_0**2
    step_volume = step_length*gb_area
    return (1/e)*step_volume*np.sum(charge_density)


class SclConcentrations(object):
    def __init__(self, c_v_scl, c_oh_scl):
        self.c_v_scl = c_v_scl
        self.c_oh_scl = c_oh_scl

class SclImpendanceParams(object):
    def __init__(self, resistivity_scl, z_scl):
        self.resistivity_scl = resistivity_scl
        self.z_scl = z_scl


# temperatures = [373, 398, 423, 448, 473, 498, 523, 548, 573, 598, 623, 648, 673, 698]
temperatures = [423]

for T in temperatures:

    phi = np.zeros(N)
    region_width = 20e-9
    x_range = np.linspace(-region_width/2, region_width/2, N)

    from fit_segregation_energies  import *
    popt_c_v, popt_c_oh = get_fit_params(T)

    e_seg_c_v = gaussian(x_range, *popt_c_v)
    # popt_c_oh[2] = 1*popt_c_oh[2] #Adjustment for smearing out segergation energies
    e_seg_c_oh = gaussian(x_range, *popt_c_oh)

    scl_concentrations = get_scl_concentrations(phi, e_seg_c_v, e_seg_c_oh, c_v_bulk, c_oh_bulk, T)
    charge_density = get_charge_density(scl_concentrations, c_acc)
    tot_chg_e_units = get_tot_chg_e_units(charge_density)

    # for i in range(0, 1000):
    i = 0
    while tot_chg_e_units > convergence_criterion:
        phi = solve_poissons_equation(phi, charge_density, n_iterations, region_width)
        scl_concentrations = get_scl_concentrations(phi, e_seg_c_v, e_seg_c_oh, c_v_bulk, c_oh_bulk, T)
        charge_density = get_charge_density(scl_concentrations, c_acc)
        # print("chg pr formula unit = ", np.average(scl_concentrations.c_oh_scl) - c_acc)
        #(1/e)*self.stepVolume*sum(self.getChargeDensityListJacobi())

        tot_chg_e_units = get_tot_chg_e_units(charge_density)

        if i%10 == 0:
            print(i)
            print("chg pr formula unit = ", np.average(scl_concentrations.c_oh_scl) - c_acc)
            print("chg density e units", tot_chg_e_units)

            # plt.plot(x_range, charge_density)
            # plt.show()
            # plt.plot(x_range, scl_concentrations.c_oh_scl, label="c_oh")
            # plt.plot(x_range, scl_concentrations.c_v_scl, label="c_v")
            # plt.legend(loc="best")
            # plt.yscale("log")
            # plt.show()
            # plt.plot(x_range, phi)
            # plt.show()
        i = i + 1


    np.savetxt("T=" + str(T) + "phi.csv", phi, delimiter=',')
    np.savetxt("T=" + str(T) + "x_range.csv", x_range, delimiter=',')
    np.savetxt("T=" + str(T) + "c_oh_scl.csv", scl_concentrations.c_oh_scl, delimiter=',')

    print("Caution: the function get_scl_impedance_params(phi, w_list).z_scl gives the wrong impedance")
    # np.savetxt("w_list.csv", w_list, delimiter=',')
    # np.savetxt("z_scl_imag.csv", -get_scl_impedance_params(phi, w_list).z_scl.imag, delimiter=',')

    char_freq_ratio = get_char_freq_ratio(get_scl_impedance_params(phi, w_list).z_scl)
    print(max(phi))
    print(char_freq_ratio)

# plt.plot(x_range, scl_concentrations.c_oh_scl)
# plt.show()
