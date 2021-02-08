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
N_conf = 6

#Numerical input parameters
x_step_length_approx = 1e-12
phi_0_incremenent = 1e-4
w_start = -4
w_end = 12
w_steps = 5000
w_list = np.logspace(w_start, w_end, w_steps)
resistivity_bulk = 1


def get_scl_length(phi_0):
    return np.sqrt((2*epsilon_0*epsilon_r*phi_0)/(e*(c_acc/a_0**3)))

def get_phi_params(exp_resistance_ratio, shrinkage_factor=1, effective_gb_area=1):
    scl_resistance_ratio = -1
    phi_0 = 0.0001
    while scl_resistance_ratio < exp_resistance_ratio:
        phi_0 = phi_0 + phi_0_incremenent
        scl_length = get_scl_length(phi_0)

        n_x_steps = int(round(scl_length/x_step_length_approx))
        x_scl_list = np.linspace(0, scl_length, n_x_steps)
        x_step_length = scl_length/n_x_steps
        phi_scl = phi_0*(x_scl_list/scl_length - 1)**2

        if shrinkage_factor != 1:
            n_x_steps = int(round(n_x_steps*shrinkage_factor))
        scl_resistivity_ratio = (1/n_x_steps)*np.sum(np.exp((e*phi_scl[:n_x_steps])/(k_B*T)))
        scl_resistance_ratio = get_resistance_ratio(scl_resistivity_ratio, scl_length, shrinkage_factor, effective_gb_area)
    phi_params = PhiParams(phi_scl, phi_0, scl_length, x_scl_list, x_step_length, n_x_steps, scl_resistivity_ratio)
    return phi_params

def get_phi_scl_continuous(phi_0, scl_length):
    n_x_steps = int(round(scl_length/x_step_length_approx))
    x_scl_list = np.linspace(0, scl_length, n_x_steps)
    phi_scl_continuous = phi_0*(x_scl_list/scl_length - 1)**2
    return x_scl_list, phi_scl_continuous

def get_phi_params_discrete(exp_resistance_ratio, steps_per_a_0=2):
    scl_resistance_ratio = -1
    # phi_0 = 0.0001
    phi_0 = 0.01
    while scl_resistance_ratio < exp_resistance_ratio:
        phi_0 = phi_0 + phi_0_incremenent
        scl_length = get_scl_length(phi_0)

        step_length = a_0/steps_per_a_0
        n_x_steps = int(scl_length/step_length)
        x_scl_list = np.linspace(0, n_x_steps*step_length, n_x_steps)
        phi_scl = phi_0*(x_scl_list/scl_length - 1)**2

        #Quick fix for super short SCLs in start of while loop, to avoid division by zero
        if n_x_steps < 2:
            n_x_steps = 2

        scl_resistivity_ratio = (1/(n_x_steps-1))*np.sum(np.exp((e*phi_scl[1:])/(k_B*T)))
        scl_resistance_ratio = get_resistance_ratio(scl_resistivity_ratio, scl_length, shrinkage_factor, effective_gb_area)
        x_step_length = scl_length/n_x_steps
        x_scl_list_continuous, phi_scl_continuous = get_phi_scl_continuous(phi_0, scl_length)

    phi_params = PhiParams(phi_scl, phi_0, scl_length, x_scl_list, x_step_length, n_x_steps, scl_resistivity_ratio, x_scl_list_continuous, phi_scl_continuous)
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

def get_scl_impedance_params_discrete(phi_params, w_list):
    resistivity_scl = resistivity_bulk*np.exp((e*phi_params.phi_scl[1:])/(k_B*T))
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

def get_resistance_ratio(scl_resistivity_ratio, scl_length, shrinkage_factor=1, effective_gb_area=1):
    return (2*scl_length*shrinkage_factor/d_grain)*scl_resistivity_ratio*(1/effective_gb_area)

def get_c_acc_with_trapping(trapping_energy):
    #Wolfram alpha: https://www.wolframalpha.com/input/?i=K+%3D+3*(c-x)%2F(N*x%5E2)+for+x
    K = np.exp((-e*trapping_energy)/(k_B*T))
    # c_acc = (-1 + sqrt(1 + 4*N_conf*K*c_acc_no_trapping))/(2*N_conf*K) #Old, error
    c_acc = (2*sqrt(3)*c_acc_no_trapping)/(sqrt(4*c_acc_no_trapping*K*N_conf + 3) + sqrt(3))
    return c_acc.real

def get_impedance_from_constant_resistance():
    return resistance_gb/(1 + j*w_list*resistance_gb*capacitance_gb)

def get_c_oh_scl(phi_scl):
    return c_oh_bulk*(1/np.exp((e*phi_scl)/(k_B*T)))

class SCLImpendanceParams(object):
    def __init__(self, resistivity_scl, z_scl):
        self.resistivity_scl = resistivity_scl
        self.z_scl = z_scl

class PhiParams(object):
    def __init__(self, phi_scl, phi_0, scl_length, x_scl_list, x_step_length, n_x_steps, scl_resistivity_ratio, x_scl_list_continuous=0, phi_scl_continuous=0):
        self.phi_scl = phi_scl
        self.phi_0 = phi_0
        self.scl_length = scl_length
        self.x_scl_list = x_scl_list
        self.x_step_length = x_step_length
        self.n_x_steps = n_x_steps
        self.scl_resistivity_ratio = scl_resistivity_ratio
        self.x_scl_list_continuous = x_scl_list_continuous
        self.phi_scl_continuous = phi_scl_continuous

#Physical input parameters
c_acc = 0.20
# resistance_ratios = [31.21664, 18.67118, 12.69875, 8.86158, 6.67949, 4.90698]
# temperatures = [373, 398, 423, 448, 473, 498]
resistance_ratios = [12.69875]
temperatures = [423]
epsilon_r = 58
d_grain = 0.40e-6

# resistivity_bulk = 5.5e1
# resistance_gb = 1.9e4
# capacitance_gb = 1.9e-8

#Modelling input parameters
shrinkage_factor = 1
effective_gb_area = 1
c_oh_bulk = c_acc

# trapping_energy = 0
# c_acc_no_trapping = c_acc
# c_acc = get_c_acc_with_trapping(trapping_energy)
# print("c_acc", c_acc)
# c_accs.append(c_acc.real)

#Run code
char_freq_ratios = []
phi_0s = []
scl_lengths = []
z = 0


# exp_resistance_ratios = [1e-1, 1, 5, 10, 15, 20, 25]
# exp_resistance_ratios = np.random.normal(12.7, 1, 100)
# exp_resistance_ratios = np.linspace(1, 500, 200)

# c_accs = []
# trapping_energies = np.linspace(-0.6, 0.2, 41)
# trapping_energies = [-0.30]



for i, T in enumerate(temperatures):
    exp_resistance_ratio = resistance_ratios[i]

    # c_acc = get_c_acc_with_trapping(trapping_energy)
    # c_accs.append(c_acc)

    # phi_params = get_phi_params_discrete(exp_resistance_ratio, 3)
    # scl_impedance_params = get_scl_impedance_params_discrete(phi_params, w_list)

    phi_params = get_phi_params(exp_resistance_ratio)
    scl_impedance_params = get_scl_impedance_params(phi_params, w_list)

    char_freq_ratios.append(get_char_freq_ratio(scl_impedance_params.z_scl))
    phi_0s.append(phi_params.phi_0)
    scl_lengths.append(phi_params.scl_length)
    print("T=", T, "resistance ratio=", exp_resistance_ratio, "w ratio=", get_char_freq_ratio(scl_impedance_params.z_scl))
    print("phi_0", phi_params.phi_0)
    c_oh_scl = get_c_oh_scl(phi_params.phi_scl)

    # plt.plot(phi_params.x_scl_list, c_oh_scl)
    # plt.yscale("log")
    # plt.show()


    np.savetxt("z_scl.csv", -scl_impedance_params.z_scl.imag, delimiter=',')
    np.savetxt("w_list.csv", w_list, delimiter=',')

    # np.savetxt("phi_scl.csv", phi_params.phi_scl_continuous.real, delimiter=',')
    # np.savetxt("x_scl_list_continuous.csv", phi_params.x_scl_list_continuous.real, delimiter=',')
    # np.savetxt("x_scl_list.csv", phi_params.x_scl_list.real, delimiter=',')
    # np.savetxt("c_oh_scl.csv", c_oh_scl, delimiter=',')

    # '''Uncomment to add impedances, in case we want to introduce inhomogeneity'''
    # z = z + get_impedance_from_constant_resistance()
    # # z = z + scl_impedance_params.z_scl
    #
    # # Normalize resistance for comparison with expt. data
    # z = z*(resistance_gb/z.real[0])



# np.savetxt("phi_0s.csv", phi_0s)
# np.savetxt("trapping_energies.csv", trapping_energies)
# np.savetxt("char_freq_ratios.csv", char_freq_ratios)
# np.savetxt("c_accs.csv", np.asarray(c_accs))


# # z = scl_impedance_params.z_scl
# np.savetxt("real_impedance.csv", z.real)
# np.savetxt("imag_impedance.csv", z.imag)
# np.savetxt("w_list.csv", w_list)
# print("z.real[0]", z.real[0])
#
#
# z_derivative = np.diff(np.log10(1/z.real))*w_steps/(w_end - w_start)
# z_derivative = np.insert(z_derivative, 0, 1e-9)
#
# fig, ax1 = plt.subplots()
# ax1.plot(np.log10(w_list), np.log10(1/z.real), 'b-')
# ax1.set_ylabel("lg 1/Z'", color='b')
# ax1.tick_params('y', colors='b')
# ax1.set_xlabel('lg $\omega$')
#
# ax2 = ax1.twinx()
# ax2.plot(np.log10(w_list), z_derivative, 'r-', label="Slope")
# ax2.set_ylabel("Slope, lg 1/Z'", color='r')
# ax2.tick_params('y', colors='r')
# plt.savefig("impedance_slope.png", dpi=250)
# plt.show()


'''Save data to csv'''
# data_to_csv = zip(char_freq_ratios, phi_0s)
# with open("char_freq_data.csv", "w") as f:
#     writer = csv.writer(f)
#     for row in data_to_csv:
#         writer.writerow(row)




'''Stored input params with samples from literature'''
# Chen data to test different modifications to SC model



# #Chen, BZY5:
# c_acc = 0.05
# resistance_ratios = [33038.02, 19760.58, 11819.13, 6378.64, 4016.39, 2802.76, 1551.99, 1169.83]
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
# resistance_ratios = [31.21664, 18.67118, 12.69875, 8.86158, 6.67949, 4.90698]
# temperatures = [373, 398, 423, 448, 473, 498]
# epsilon_r = 58
# d_grain = 0.4e-6


# #My own experimental data
# c_acc = 0.15
# resistivity_ratios = [1/3.63e-6, 1/5.22e-6, 1/3.13e-5, 1/8.20e-5]
# temperatures = [373, 423, 473, 523]
# epsilon_r = 60
# d_grain = 0.50e-6 #Dummy

#In the case when we want model to match expt data ish-exactly
# c_acc = 0.15
# exp_resistance_ratios = [2.7e1]
# T = 523
# epsilon_r = 60
# d_grain = 0.50e-6 #Dummy
# resistivity_bulk = 5.5e1
# resistance_gb = 1.9e4
# capacitance_gb = 1.9e-8



'''Start for loop:
for i in range(0, len(resistivity_ratios)):
    T = temperatures[i]
'''



#
