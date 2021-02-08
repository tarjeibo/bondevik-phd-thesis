from matplotlib import pyplot as plt
import numpy as np
from cmath import sqrt

#Physical constants
e = 1.60217662e-19
epsilon_0 = 8.85e-12
j = sqrt(-1)
k_B = 1.38064852e-23


def get_phi_scl_gaussian(x_list, mu, sigma, phi_0):
    return phi_0*np.exp(-np.power(x_list - mu, 2.) / (2 * np.power(sigma, 2.)))

def get_second_derivative(array):
    array = np.diff(np.diff(array))
    array = np.insert(array, 0, array[0])
    array = np.append(array, array[-1])
    return array/x_step_size**2

# def get_c_oh_scl(charge_density_per_formula_unit):
#     return charge_density_per_formula_unit + c_acc

def get_c_oh_scl(gaussian_phi_scl):
    return c_oh_bulk*(1/np.exp((e*gaussian_phi_scl)/(k_B*T)))

def get_z_scl(resistivity_scl):
    z_scl = []
    for w in w_list:
        z_scl_at_w = np.sum((resistivity_scl*x_step_size)/(1 + j*w*resistivity_scl*epsilon_0*epsilon_r))
        z_scl.append(z_scl_at_w)
    z_scl = np.asarray(z_scl)
    return z_scl

def get_char_freq_ratio(z_scl):
    w0_scl_idx = np.argmax(-z_scl.imag)
    w0_scl = w_list[w0_scl_idx]
    w0_bulk = 1/(resistivity_bulk*epsilon_0*epsilon_r)
    return w0_scl/w0_bulk

def get_resistance_ratio(gaussian_phi_scl, T):
    return (2*x_range/d_grain)*np.average(np.exp((e*gaussian_phi_scl)/(k_B*T)))


x_range = 40e-9
x_steps = 2001
x_step_size = x_range/x_steps
x_list = np.linspace(-x_range/2, x_range/2, x_steps)
mu = 0
# sigmas = np.linspace(1e-10, 3.5e-9, 35)
sigmas = [2e-9]
phi_0 = 0.2
phi_0_incremenent = 1e-2
T = 423
c_oh_bulk = 0.2
c_acc = 0.2
epsilon_r = 75
a_0 = 4.2e-10
d_grain = 0.4e-6
expt_resistance_ratio = 12.69875

w_start = -4
w_end = 12
w_steps = 5000
w_list = np.logspace(w_start, w_end, w_steps)
resistivity_bulk = 1

# charge_density = -get_second_derivative(gaussian_phi_scl)*(epsilon_0*epsilon_r)
# charge_density_per_formula_unit = (charge_density*a_0**3)/e
char_freq_ratios = []
for sigma in sigmas:
    phi_0 = 0.2
    resistance_ratio = -1 #dummy value
    while resistance_ratio < expt_resistance_ratio:
        gaussian_phi_scl = get_phi_scl_gaussian(x_list, mu, sigma, phi_0)
        c_oh_scl = get_c_oh_scl(gaussian_phi_scl)
        resistivity_scl = resistivity_bulk*(c_oh_bulk/c_oh_scl)
        resistance_ratio = get_resistance_ratio(gaussian_phi_scl, T)
        phi_0 += phi_0_incremenent
        # print("phi_0", phi_0)


    z_scl = get_z_scl(resistivity_scl)
    char_freq_ratios.append(get_char_freq_ratio(z_scl))

# #Plot sigmas vs char freq
# plt.plot(sigmas, char_freq_ratios, label="gaussian space charge model", marker='o')
# plt.yscale("log")
# plt.axhline(y=9.65e-4, color="r", linestyle="dotted", label="expt value")
# plt.xlabel("Standard deviation, GB potential")
# plt.ylabel("$\omega_{0,gb}$ / $\omega_{0,bulk}$")
# plt.legend(loc="best")
# plt.savefig("foo.png", dpi=300)
# plt.show()



plt.plot(x_list/1e-9, gaussian_phi_scl)
plt.xlabel("x / nm")
plt.ylabel("Potential / V")
plt.savefig("foo.png", dpi=300)
plt.show()


# plt.plot(resistivity_scl)
# plt.plot(gaussian_phi_scl)
# plt.show()
#
# print(z_scl)
# print(np.average(charge_density_per_formula_unit))
# plt.plot(z_scl.real, -z_scl.imag)
# plt.show()

# plt.plot(x, gaussian_phi_scl)
# plt.plot(x, charge_density_per_formula_unit)
# plt.show()
