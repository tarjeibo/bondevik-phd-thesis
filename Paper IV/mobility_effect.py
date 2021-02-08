import numpy as np
import csv
from cmath import sqrt
from matplotlib import pyplot as plt

#Physical constants
e = 1.60217662e-19
epsilon_0 = 8.85e-12
j = sqrt(-1)
k_B = 1.38064852e-23

#BaZrO3 material specific constants
a_0 = 4.22e-10

#Numerical input parameters
x_step_length_approx = 1e-12
phi_0_incremenent = 1e-4
w_list = np.logspace(-4, 12, 5000)
resistivity_bulk = 1


def get_scl_length(phi_0):
    return np.sqrt((2*epsilon_0*epsilon_r*phi_0)/(e*(c_acc/a_0**3)))

def get_core_resistivity_ratio(scl_length):
    return (delta_core/(delta_core + 2*scl_length))*(1/mobility_ratio)
    # return (delta_core/delta_core)*(1/mobility_ratio)

def get_phi_params(exp_resistance_ratio):
    model_resistance_ratio = -1
    phi_0 = 0.0001
    while model_resistance_ratio < exp_resistance_ratio:
        phi_0 = phi_0 + phi_0_incremenent
        scl_length = get_scl_length(phi_0)

        n_x_steps = int(round(scl_length/x_step_length_approx))
        x_scl_list = np.linspace(0, scl_length, n_x_steps)
        x_step_length = scl_length/n_x_steps
        phi_scl = phi_0*(x_scl_list/scl_length - 1)**2

        core_resistivity_ratio = get_core_resistivity_ratio(scl_length)
        core_resistance_ratio = (delta_core/d_grain)*core_resistivity_ratio
        # print("core_resistance_ratio", core_resistance_ratio)

        scl_resistance_ratio = (2*scl_length/d_grain)*(1/n_x_steps)*np.sum(np.exp((e*phi_scl)/(k_B*T)))
        model_resistance_ratio = scl_resistance_ratio + core_resistance_ratio

    phi_params = PhiParams(phi_scl, core_resistivity_ratio, phi_0, x_step_length, core_resistance_ratio, scl_resistance_ratio)
    return phi_params

def get_model_impedance(phi_params, w_list):
    resistivity_scl = resistivity_bulk*np.exp((e*phi_params.phi_scl)/(k_B*T))
    resistivity_core = resistivity_bulk*phi_params.core_resistivity_ratio
    z_model = []
    for w in w_list:
        z_scl_at_w = np.sum((resistivity_scl*phi_params.x_step_length)/(1 + j*w*resistivity_scl*epsilon_0*epsilon_r))
        z_core_at_w = (resistivity_core*delta_core)/(1 + j*w*resistivity_core*epsilon_0*epsilon_r)
        z_model.append(z_scl_at_w + z_core_at_w)
    z_model = np.asarray(z_model)
    return z_model

def get_char_freq_ratio(z_model):
    w0_model_idx = np.argmax(-z_model.imag)
    w0_model = w_list[w0_model_idx]
    w0_bulk = 1/(resistivity_bulk*epsilon_0*epsilon_r)
    return w0_model/w0_bulk


class PhiParams(object):
    def __init__(self, phi_scl, core_resistivity_ratio, phi_0, x_step_length, core_resistance_ratio, scl_resistance_ratio):
        self.phi_scl = phi_scl
        self.core_resistivity_ratio = core_resistivity_ratio
        self.phi_0 = phi_0
        self.x_step_length = x_step_length
        self.core_resistance_ratio = core_resistance_ratio
        self.scl_resistance_ratio = scl_resistance_ratio


#Physical input parameters
c_acc = 0.20
exp_resistance_ratio = 12.69875
T = 423
epsilon_r = 58
d_grain = 0.4e-6
delta_core = 0.5e-9
activation_energies = np.linspace(0, 0.35, 36)
print(activation_energies)
mobility_ratios = np.exp(-(activation_energies*e)/(k_B*T))
print(mobility_ratios)


char_freq_ratios = []
core_resistance_ratios = []
scl_resistance_ratios = []
for mobility_ratio in mobility_ratios:
    phi_params = get_phi_params(exp_resistance_ratio)
    z_model = get_model_impedance(phi_params, w_list)

    char_freq_ratios.append(get_char_freq_ratio(z_model))
    core_resistance_ratios.append(phi_params.core_resistance_ratio)
    scl_resistance_ratios.append(phi_params.scl_resistance_ratio)

scl_resistance_ratios = np.asarray(scl_resistance_ratios)
core_resistance_ratios = np.asarray(core_resistance_ratios)
char_freq_ratios = np.asarray(char_freq_ratios)


np.savetxt("activation_energies.csv", activation_energies, delimiter=",")
np.savetxt("scl_resistance_ratios.csv", scl_resistance_ratios, delimiter=",")
np.savetxt("core_resistance_ratios.csv", core_resistance_ratios, delimiter=",")
np.savetxt("char_freq_ratios.csv", char_freq_ratios, delimiter=",")

#
#
# def plot_char_freq_ratios():
#     plt.plot(activation_energies, char_freq_ratios, marker='o')
#     plt.xlabel("(Additional) activation energy GB core / eV")
#     plt.ylabel("Characteristic frequency")
#     plt.yscale("log")
#     plt.savefig("char_freq_mobility.png", dpi=400)
#
# def plot_core_and_scl_resistances():
#     plt.plot(activation_energies, scl_resistance_ratios, label="scl_resistance", marker='o')
#     plt.plot(activation_energies, core_resistance_ratios, label="core_resistance", marker='o')
#     plt.plot(activation_energies, scl_resistance_ratios + core_resistance_ratios, label="total_resistance", marker='o')
#     plt.legend(loc="best")
#     plt.xlabel("(Additional) activation energy GB core / eV")
#     plt.ylabel("Resistance / a.u.")
#     plt.ylim(0, 20)
#     plt.savefig("core_scl_resistance_mobility.png", dpi=400)
#
#
#
# # plot_core_and_scl_resistances()
#
# # plot_char_freq_ratios()
#
# print(core_resistance_ratios + scl_resistance_ratios)
#
#
# activation_energies = np.linspace(0.30, 0.31, 2)
# print(activation_energies)
# mobility_ratios = np.exp(-(activation_energies*e)/(k_B*T))
# phi_params = get_phi_params(exp_resistance_ratio)
# z_model = get_model_impedance(phi_params, w_list)
# plt.plot(z_model.real, -z_model.imag)
# x=3
# print("phi_0", phi_params.phi_0)
# print("w", char_freq_ratios)
# print("R_core", core_resistance_ratios)
# print("R_scl", scl_resistance_ratios)

# plt.plot(z_model.real, -z_model.imag)
# plt.show()


# data_to_csv = zip(mobility_ratios, char_freq_ratios, core_resistance_ratios, scl_resistance_ratios)
# with open("mobility_data.csv", "w") as f:
#     writer = csv.writer(f)
#     for row in data_to_csv:
#         writer.writerow(row)


# print(char_freq_ratios)
# print(core_resistance_ratios)
# print(scl_resistance_ratios)
# print(np.asarray(core_resistance_ratios) + np.asarray(scl_resistance_ratios))
