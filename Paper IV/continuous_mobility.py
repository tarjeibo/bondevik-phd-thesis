#Import proton concentration and segregation energy, and output mobility adjusted resistivity and characteristic frequency
import numpy as np
from matplotlib import pyplot as plt
from cmath import sqrt

epsilon_0 = 8.854187817620e-12
k_B_eV = 8.6173325e-5
T = 423
j = sqrt(-1)
epsilon_r = 75
resistivity_bulk = 1

w_start = -4
w_end = 12
w_steps = 5000
w_list = np.logspace(w_start, w_end, w_steps)

#Input data

x_range = np.genfromtxt('MobilityCalcs/' + 'T=' + str(T) + 'x_range.csv', delimiter=',')
c_oh_scl = np.genfromtxt('MobilityCalcs/' + 'T=' + str(T) + 'c_oh_scl.csv', delimiter=',')

x_range_e_seg = np.genfromtxt('MobilityCalcs/' + 'x_for_plot_fit.csv', delimiter=',')
e_seg_c_oh = np.genfromtxt('MobilityCalcs/' + 'e_seg_c_oh_for_plot.csv', delimiter=',')

max_proton_jump_distance = 1.58e-10
perp_proton_jump_distance = max_proton_jump_distance*(2/np.pi)
# perp_proton_jump_distance = max_proton_jump_distance

def get_energy_barrier():
    energy_barrier = []
    for initial_idx, x in enumerate(x_range_e_seg):
        final_idx = (np.abs(x_range_e_seg - (x + perp_proton_jump_distance))).argmin()
        energy_barrier.append(e_seg_c_oh[final_idx] - e_seg_c_oh[initial_idx])
    energy_barrier = np.asarray(energy_barrier).clip(min=0)
    energy_barrier[0], energy_barrier[-1] = 0, 0
    return energy_barrier

def get_local_mobility_ratio():
    local_mobility_ratio_e_seg_idxs = np.exp(energy_barrier/(k_B_eV*T))
    local_mobility_ratio = []
    for x in x_range:
        nearest_e_seg_idx = (np.abs(x_range_e_seg - x)).argmin()
        local_mobility_ratio.append(local_mobility_ratio_e_seg_idxs[nearest_e_seg_idx])
    return np.asarray(local_mobility_ratio)


def get_local_resistivity():
    c_oh_bulk = c_oh_scl[0]
    return (c_oh_bulk/c_oh_scl)*local_mobility_ratio

def get_local_resistivity_excluding_mobility():
    c_oh_bulk = c_oh_scl[0]
    return (c_oh_bulk/c_oh_scl)

def get_z_scl(local_resistivity):
    step_length = x_range[1] - x_range[0]
    z_scl = []
    for w in w_list:
        z_scl_at_w = np.sum((local_resistivity*step_length)/(1 + j*w*local_resistivity*epsilon_0*epsilon_r))
        z_scl.append(z_scl_at_w)
    z_scl = np.asarray(z_scl)
    return z_scl

def get_char_freq_ratio(z_scl):
    w0_scl_idx = np.argmax(-z_scl.imag)
    w0_scl = w_list[w0_scl_idx]
    w0_bulk = 1/(resistivity_bulk*epsilon_0*epsilon_r)
    return w0_scl/w0_bulk


energy_barrier = get_energy_barrier()
local_mobility_ratio = get_local_mobility_ratio()

local_resistivity = get_local_resistivity()
local_resistivity_excluding_mobility = get_local_resistivity_excluding_mobility()

z_scl = get_z_scl(local_resistivity)
w_0 = get_char_freq_ratio(z_scl)
print(w_0)

z_scl_no_mu = get_z_scl(local_resistivity_excluding_mobility)
w_0_no_mu = get_char_freq_ratio(z_scl_no_mu)
print(w_0_no_mu)

# plt.plot(x_range_e_seg, e_seg_c_oh)
# plt.plot(x_range_e_seg, energy_barrier)


plt.plot(x_range, local_resistivity)
plt.plot(x_range, local_resistivity_excluding_mobility)
plt.xlim(-3e-9, 3e-9)
plt.yscale("log")
plt.show()



# plt.plot(x_range, c_oh_scl)
# plt.xlim(-3e-9, 3e-9)
# plt.yscale("log")
# plt.show()













# more space
