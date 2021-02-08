import numpy as np
import csv
import sympy as sp
from cmath import sqrt
from matplotlib import pyplot as plt

#Physical constants
e = 1.60217662e-19
epsilon_0 = 8.85e-12
j = sqrt(-1)
k_B = 1.38064852e-23

#BaZrO3 material specific constants
a_0 = 4.22e-10
N_conf = 1

#Numerical input parameters
x_step_length_approx = 1e-12
phi_0_incremenent = 1e-4
w_list = np.logspace(-4, 12, 5000)
resistivity_bulk = 1


def get_scl_length(phi_0):
    return np.sqrt((2*epsilon_0*epsilon_r*phi_0)/(e*(c_acc/a_0**3)))

def get_phi_params(exp_resistivity_ratio, shrinkage_factor=1, effective_gb_area=1):
    scl_resistivity_ratio = -1
    phi_0 = 0.0001
    while scl_resistivity_ratio < exp_resistivity_ratio*effective_gb_area:
        phi_0 = phi_0 + phi_0_incremenent
        scl_length = get_scl_length(phi_0)

        n_x_steps = int(round(scl_length/x_step_length_approx))
        x_scl_list = np.linspace(0, scl_length, n_x_steps)
        x_step_length = scl_length/n_x_steps
        phi_scl = phi_0*(x_scl_list/scl_length - 1)**2

        if shrinkage_factor != 1:
            n_x_steps = int(round(n_x_steps*shrinkage_factor))
        scl_resistivity_ratio = (1/n_x_steps)*np.sum(np.exp((e*phi_scl[:n_x_steps])/(k_B*T)))
    phi_params = PhiParams(phi_scl, phi_0, scl_length, x_scl_list, x_step_length, n_x_steps, scl_resistivity_ratio)
    return phi_params

def get_scl_impedance_params(phi_params, w_list):
    phi_scl = phi_params.phi_scl[:phi_params.n_x_steps]
    resistivity_scl = resistivity_bulk*np.exp((e*phi_scl)/(k_B*T))
    z_scl = []
    for w in w_list:
        z_scl_at_w = np.sum((resistivity_scl*phi_params.x_step_length)/(1 + j*w*resistivity_scl*epsilon_0*epsilon_r))
        z_scl.append(z_scl_at_w)
    z_scl = np.asarray(z_scl)
    scl_impedance_params = SCLImpendanceParams(resistivity_scl, z_scl)
    return scl_impedance_params

def get_char_freq_ratio(z_scl):
    w0_scl_idx = np.argmax(-z_scl.imag)
    w0_scl = w_list[w0_scl_idx]
    w0_bulk = 1/(resistivity_bulk*epsilon_0*epsilon_r)
    return w0_scl/w0_bulk

def get_resistance_ratio(phi_params, shrinkage_factor=1, effective_gb_area=1):
    return (2*phi_params.scl_length*shrinkage_factor/d_grain)*phi_params.scl_resistivity_ratio*(1/effective_gb_area)

def get_c_acc_with_trapping(trapping_energy):
    K = np.exp((-e*trapping_energy)/(k_B*T))
    c_acc = (-1 + sqrt(1 + 4*N_conf*K*c_acc_no_trapping))/(2*N_conf*K)
    return c_acc


class SCLImpendanceParams(object):
    def __init__(self, resistivity_scl, z_scl):
        self.resistivity_scl = resistivity_scl
        self.z_scl = z_scl

class PhiParams(object):
    def __init__(self, phi_scl, phi_0, scl_length, x_scl_list, x_step_length, n_x_steps, scl_resistivity_ratio):
        self.phi_scl = phi_scl
        self.phi_0 = phi_0
        self.scl_length = scl_length
        self.x_scl_list = x_scl_list
        self.x_step_length = x_step_length
        self.n_x_steps = n_x_steps
        self.scl_resistivity_ratio = scl_resistivity_ratio

#Physical input parameters
c_acc = 0.10
exp_resistivity_ratio = (1/0.00442166162797173)
# exp_resistivity_ratio = (1/0.00100960210518789)
print("exp_resistivity_ratio", exp_resistivity_ratio)
T = 600
epsilon_r = 75
d_grain = 0.4e-6


#Modelling input parameters
shrinkage_factor = 1
effective_gb_area = 1
# trapping_energy = -0.14
# c_acc_no_trapping = c_acc
# c_acc = get_c_acc_with_trapping(trapping_energy)

#Run code
char_freq_ratios = []
phi_0s = []
scl_lengths = []
resistance_ratios = []
# effective_gb_areas = np.linspace(1/20, 1, 100)
effective_gb_areas = [1]
for effective_gb_area in effective_gb_areas:
    # T = temperatures[i]

    # exp_resistivity_ratio = resistivity_ratios[i]
    phi_params = get_phi_params(exp_resistivity_ratio, shrinkage_factor, effective_gb_area)
    scl_impedance_params = get_scl_impedance_params(phi_params, w_list)
    char_freq_ratios.append(get_char_freq_ratio(scl_impedance_params.z_scl))
    phi_0s.append(phi_params.phi_0)
    scl_lengths.append(phi_params.scl_length)

    # np.savetxt("phi_0.csv", phi_params.phi_scl.real, delimiter=',')
    # np.savetxt("x_scl_list.csv", phi_params.x_scl_list.real, delimiter=',')

    resistance_ratios.append(get_resistance_ratio(phi_params, shrinkage_factor, effective_gb_area))
    print("T=", T, "RR=", exp_resistivity_ratio, "ratio=", get_char_freq_ratio(scl_impedance_params.z_scl))
    print("phi_0=", phi_params.phi_0)


# data_to_csv = zip(temperatures, resistivity_ratios, char_freq_ratios, phi_0s, resistance_ratios, scl_lengths)
data_to_csv = zip(effective_gb_areas, char_freq_ratios, phi_0s, resistance_ratios, scl_lengths)
with open("char_freq_data.csv", "w") as f:
    writer = csv.writer(f)
    for row in data_to_csv:
        writer.writerow(row)




'''Stored input params with samples from literature'''
# Chen data to test different modifications to SC model





# #Chen, BZY5:
# c_acc = 0.05
# resistivity_ratios = [6.39e5, 3.82e5, 2.29e5, 1.23e5, 7.77e4, 5.42e4, 3.00e4, 2.26e4]
# temperatures = [523, 548, 573, 598, 623, 648, 673, 698]
# epsilon_r = 66
# d_grain = 0.24e-6


# #Chen, BZY8:
# c_acc = 0.08
# resistivity_ratios = [4.05e3, 3.05e3, 1.6e3, 1.06e3]
# temperatures = [473, 498, 523, 548]
# epsilon_r = 85
# d_grain = 0.09e-6


# #Chen, BZY20:
# c_acc = 0.20
# resistivity_ratios = [2.55e3, 1.52e3, 1.04e3, 7.23e2, 5.45e2, 4.01e2]
# temperatures = [373, 398, 423, 448, 473, 498]
# epsilon_r = 58
# d_grain = 0.4e-6


# #My own experimental data
# c_acc = 0.15
# resistivity_ratios = [1/3.63e-6, 1/5.22e-6, 1/3.13e-5, 1/8.20e-5]
# temperatures = [373, 423, 473, 523]
# epsilon_r = 60
# d_grain = 0.50e-6 #Dummy



'''Start for loop:
for i in range(0, len(resistivity_ratios)):
    T = temperatures[i]
'''



#
