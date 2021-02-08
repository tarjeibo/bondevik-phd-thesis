#Script that inputs potential and outputs char. frequency

import numpy as np
from cmath import sqrt
from matplotlib import pyplot as plt

#Fundamental constants
epsilon_0 = 8.854187817620e-12
k_B = 1.3806488e-23
k_B_eV = 8.6173325e-5
e = 1.602176e-19
j = sqrt(-1)

#Numerical params
w_start = -4
w_end = 12
w_steps = 5000
w_list = np.logspace(w_start, w_end, w_steps)
phi_0_incremenent = 0.001

#Experimental params
resistivity_bulk = 1
T = 423
epsilon_r = 75
c_acc = 0.1
c_oh_bulk = 0.1
a_0 = 4.2e-10



def get_scl_impedance_params_from_concentration(x_range, c_oh_scl):
    resistivity_scl = resistivity_bulk*(c_oh_bulk/c_oh_scl)
    step_length = x_range[1] - x_range[0]
    z_scl = []
    for w in w_list:
        z_scl_at_w = np.sum((resistivity_scl*step_length)/(1 + j*w*resistivity_scl*epsilon_0*epsilon_r))
        z_scl.append(z_scl_at_w)
    z_scl = np.asarray(z_scl)
    scl_impedance_params = SclImpendanceParams(resistivity_scl, z_scl)
    return scl_impedance_params


def get_scl_impedance_params_from_phi(x_range, phi):
    resistivity_scl = resistivity_bulk*np.exp(phi/(k_B_eV*T))
    step_length = x_range[1] - x_range[0]
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

def get_scl_length(phi_0):
    return np.sqrt((2*epsilon_0*epsilon_r*phi_0)/(e*(c_acc/a_0**3)))

def get_phi_params_model_M1(R_scl, phi_0_incremenent):
    n_x_steps = 500
    M1_scl_resistance = -1
    phi_0 = 0.0001
    while M1_scl_resistance < R_scl:
        phi_0 = phi_0 + phi_0_incremenent
        scl_length = get_scl_length(phi_0)
        x_range = np.linspace(0, scl_length, n_x_steps)
        phi = phi_0*(x_range/scl_length - 1)**2
        M1_scl_resistance = 2*scl_length*resistivity_bulk*np.average(np.exp(phi/(k_B_eV*T)))
    c_oh_scl = c_oh_bulk*np.exp(-phi/(k_B_eV*T))
    phi_params = PhiParamsModelM1(phi, x_range, c_oh_scl)
    return phi_params

def get_core_charge(x_range, c_oh_scl):
    net_charge = c_oh_scl - c_acc
    positive_idxs = np.argwhere(net_charge > 1e-4)[:,0]
    core_charge_surf = e*(1/a_0**3)*np.trapz(net_charge[positive_idxs], x_range[positive_idxs])
    return core_charge_surf


class SclImpendanceParams(object):
    def __init__(self, resistivity_scl, z_scl):
        self.resistivity_scl = resistivity_scl
        self.z_scl = z_scl

class PhiParamsModelM1(object):
    def __init__(self, phi, x_range, c_oh_scl):
        self.phi = phi
        self.x_range = x_range
        self.c_oh_scl = c_oh_scl




# temperatures = [373, 398, 423, 448, 473, 498, 523, 548, 573, 598, 623, 648, 673, 698]
temperatures = [423]
char_freq_ratios_smooth = []
char_freq_ratios_M1 = []
for T in temperatures:

    # np.savetxt("T=" + str(T) + "phi.csv", phi, delimiter=',')

    #Input data
    x_range = np.genfromtxt("T=" + str(T) + 'x_range.csv', delimiter=',')
    phi = np.genfromtxt("T=" + str(T) + 'phi.csv', delimiter=',')
    c_oh_scl = np.genfromtxt("T=" + str(T) + 'c_oh_scl.csv', delimiter=',')

    z_scl = get_scl_impedance_params_from_concentration(x_range, c_oh_scl).z_scl
    R_scl = z_scl.real[0]
    print("\n")
    print("R_scl", R_scl)
    print("Smooth", get_char_freq_ratio(z_scl), max(phi))
    char_freq_ratios_smooth.append(get_char_freq_ratio(z_scl))


    #Input data M3
    x_range_M3 = np.genfromtxt('x_range_M3.csv', delimiter=',')
    phi_M3 = np.genfromtxt('phi_M3.csv', delimiter=',')
    c_oh_scl_M3 = np.genfromtxt('c_oh_scl_M3.csv', delimiter=',')

    z_scl_M3 = get_scl_impedance_params_from_concentration(x_range_M3, c_oh_scl_M3).z_scl
    R_scl_M3 = z_scl_M3.real[0]
    # print("\n")
    # print("R_scl_M3", R_scl_M3)
    # print("Smooth, M3", get_char_freq_ratio(z_scl_M3), max(phi_M3))


    phi_params_M1 = get_phi_params_model_M1(R_scl, phi_0_incremenent)
    z_scl_M1 = 2*get_scl_impedance_params_from_phi(phi_params_M1.x_range, phi_params_M1.phi).z_scl
    print("\n")
    print("M1", get_char_freq_ratio(z_scl_M1), max(phi_params_M1.phi))
    char_freq_ratios_M1.append(get_char_freq_ratio(z_scl_M1))

    # np.savetxt("x_range_M1.csv", np.asarray(phi_params_M1.x_range), delimiter=',')
    # np.savetxt("c_oh_scl_M1.csv", np.asarray(phi_params_M1.c_oh_scl), delimiter=',')

    # print("core charge", get_core_charge(x_range, c_oh_scl))

    # plt.plot(x_range, phi)
    # plt.plot(phi_params_M1.x_range, phi_params_M1.phi)
    # plt.show()

    # plt.plot(x_range, c_oh_scl, label="Smooth")
    # plt.plot(phi_params_M1.x_range, phi_params_M1.c_oh_scl, label="M1")
    # plt.plot(x_range_M3, c_oh_scl_M3, label="M3")
    # plt.legend(loc="best")
    # plt.yscale("log")
    # plt.show()

    # plt.plot(z_scl.real, -z_scl.imag, label="Smooth")
    # plt.plot(z_scl_M1.real, -z_scl_M1.imag, label="M1")
    # plt.legend(loc="best")
    # plt.show()


np.savetxt("z_scl_imag_M3.csv", -z_scl_M3.imag, delimiter=',')
np.savetxt("z_scl_imag_M1.csv", -z_scl_M1.imag, delimiter=',')
np.savetxt("w_list.csv", w_list, delimiter=',')


np.savetxt("char_freq_smooth.csv", np.asarray(char_freq_ratios_smooth), delimiter=',')
np.savetxt("char_freq_M1.csv", np.asarray(char_freq_ratios_M1), delimiter=',')
np.savetxt("temperatures.csv", np.asarray(temperatures), delimiter=',')







#
