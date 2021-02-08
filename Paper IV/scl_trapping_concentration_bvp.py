import numpy as np
import csv
import sympy as sp
from cmath import sqrt
from matplotlib import pyplot as plt
from scipy.integrate import solve_bvp

#Physical constants
e = 1.60217662e-19
epsilon_0 = 8.85e-12
j = sqrt(-1)
k_B = 1.38064852e-23

#BaZrO3 material specific constants
a_0 = 4.22e-10
N_conf = 6
resistivity_bulk = 1 #Dummy

#Numerical params
w_start = -4
w_end = 12
w_steps = 5000
w_list = np.logspace(w_start, w_end, w_steps)
scl_length = 5e-9 #End of integration zone
n_scl_elements = 1000
x_step_length = scl_length/n_scl_elements
scl_array = np.linspace(0, scl_length, n_scl_elements)


def get_c_acc_bulk_with_trapping(c_acc_no_trapping, K_trapping):
    # c_acc_bulk = (2*sqrt(3)*c_acc_no_trapping)/(sqrt(4*c_acc_no_trapping*K_trapping*N_conf + 3) + sqrt(3)) #Old
    c_acc_bulk =  (-3 + sqrt(9 + 12*K_trapping*c_acc_no_trapping))/(2*K_trapping)
    return c_acc_bulk.real

def get_trapping_equilibrium_constant(trapping_energy, T):
    #Wolfram alpha: https://www.wolframalpha.com/input/?i=K+%3D+3*(c-x)%2F(N*x%5E2)+for+x
    return np.exp((-e*trapping_energy)/(k_B*T))

def get_scl_length(phi_0, c_acc):
    return np.sqrt((2*epsilon_0*epsilon_r*phi_0)/(e*(c_acc/a_0**3)))

def get_c_acc_scl(c_acc_no_trapping, c_oh_scl, K_trapping):
    #https://www.wolframalpha.com/input/?i=K+%3D+3*(D-y)%2F(y*h)+for+y
    return 3*c_acc_no_trapping/(c_oh_scl*K_trapping + 3)

def get_c_oh_scl_from_c_acc_scl(c_acc_no_trapping, c_acc_scl, K_trapping):
    return 3*(c_acc_no_trapping - c_acc_scl)/(K_trapping*c_acc_scl)



def bvp_solver(x, y):
    # print("len(y[0])", len(y[0]))
    return np.vstack((y[1], - B*np.exp(-(e*y[0])/(k_B*T)) + C[0:len(y[0])]))

def boundary_conditions(y_start, y_end):
    return np.array([y_start[0] - phi_0, y_end[1]])

def solve_poisson_equation_bvp():
    #Setting initial guess to y = phiGuess, and y' = 0
    y = np.zeros((2, scl_array.size))
    y[0][0] = phi_0
    y[1][0] = 0
    res = solve_bvp(bvp_solver, boundary_conditions, scl_array, y)
    phi_scl = res.sol(scl_array)[0]
    phi_scl = phi_scl - phi_scl[-1]
    return phi_scl

def get_c_oh_scl(phi_scl, T):
    return c_oh_bulk*(1/np.exp((e*phi_scl)/(k_B*T)))

def get_resistance_ratio(phi_scl, T):
    return (2*scl_length/d_grain)*np.average(np.exp((e*phi_scl)/(k_B*T)))


def get_scl_impedance_params(resistivity_scl, w_list):
    z_scl = []
    for w in w_list:
        z_scl_at_w = np.sum((resistivity_scl*x_step_length)/(1 + j*w*resistivity_scl*epsilon_0*epsilon_r))
        z_scl.append(z_scl_at_w)
    z_scl = np.asarray(z_scl)
    return z_scl

def get_char_freq_ratio(z_scl):
    w0_scl_idx = np.argmax(-z_scl.imag)
    w0_scl = w_list[w0_scl_idx]
    w0_bulk = 1/(resistivity_bulk*epsilon_0*epsilon_r)
    return w0_scl/w0_bulk

#Experimental params
epsilon_r = 58
c_acc_no_trapping = 0.20
trapping_energy = -0.30
T = 423
d_grain = 0.4e-6
exp_resistance_ratio = 12.69875

K_trapping = get_trapping_equilibrium_constant(trapping_energy, T)
c_acc_bulk = get_c_acc_bulk_with_trapping(c_acc_no_trapping, K_trapping)
print("c_acc_bulk", c_acc_bulk)
c_oh_bulk = c_acc_bulk
B = e*c_oh_bulk/(epsilon_0*epsilon_r)/a_0**3

print("c_oh_bulk", c_oh_bulk)
#Initialize model
phi_0 = 0.36
phi_0_incremenent = 1e-2
C = np.zeros(n_scl_elements)



#Uncomment to calculate space charge with BVP-method with no trapping
# C.fill(e*c_acc_no_trapping/(epsilon_0*epsilon_r)/a_0**3)
# B = e*c_acc_no_trapping/(epsilon_0*epsilon_r)/a_0**3
# scl_resistance_ratio = -1 #Dummy value
# while scl_resistance_ratio < exp_resistance_ratio:
#     phi_scl = solve_poisson_equation_bvp()
#     scl_resistance_ratio = get_resistance_ratio(phi_scl, T)
#     c_oh_scl = get_c_oh_scl(phi_scl, T)
#     print("scl_resistance_ratio", scl_resistance_ratio)
#     print("phi_0", phi_0)
#     phi_0 += phi_0_incremenent


#BVP trapping
C.fill(e*c_acc_bulk/(epsilon_0*epsilon_r)/a_0**3)
phi_scl = solve_poisson_equation_bvp()
c_oh_scl = get_c_oh_scl(phi_scl, T)
scl_resistance_ratio = -1 #Dummy value
while scl_resistance_ratio < exp_resistance_ratio:
    for i in range(5):
        c_acc_scl = get_c_acc_scl(c_acc_no_trapping, c_oh_scl, K_trapping)
        C = e*c_acc_scl/(epsilon_0*epsilon_r)/a_0**3
        phi_scl = solve_poisson_equation_bvp()
        c_oh_scl = get_c_oh_scl_from_c_acc_scl(c_acc_no_trapping, c_acc_scl, K_trapping)

    scl_resistance_ratio = get_resistance_ratio(phi_scl, T)
    print("scl_resistance_ratio", scl_resistance_ratio)
    print("phi_0", phi_0)
    phi_0 += phi_0_incremenent



resistivity_scl = resistivity_bulk*(c_oh_bulk/c_oh_scl)
z_scl = get_scl_impedance_params(resistivity_scl, w_list)
char_freq_ratio = get_char_freq_ratio(z_scl)
np.savetxt("z_scl_imag.csv", -z_scl.imag, delimiter=',')
np.savetxt("w_list.csv", w_list, delimiter=',')

# plt.plot(w_list, -z_scl.imag)
# plt.xscale("log")
# plt.show()

print("char_freq_ratio", char_freq_ratio)
print("z_scl.real[0]", z_scl.real[0])

# np.savetxt("scl_array_le_chat_trapping.csv", scl_array, delimiter=',')
# np.savetxt("phi_scl_le_chat_trapping.csv", phi_scl, delimiter=',')
# np.savetxt("c_oh_scl_le_chat_trapping.csv", c_oh_scl/c_oh_bulk, delimiter=',')
# np.savetxt("c_acc_scl_le_chat_trapping.csv", c_acc_scl, delimiter=',')

# plt.plot(scl_array, c_oh_scl/c_oh_bulk)
# plt.yscale("log")
# plt.show()



# phi_scl_2 = solve_poisson_equation_bvp()
# c_oh_scl_2 = get_c_oh_scl(phi_scl_2, T)
#
# plt.plot(c_oh_scl)
# plt.plot(c_oh_scl_2)
# plt.show()


# phi_scl = solve_poisson_equation_bvp()
# c_oh_scl_2 = get_c_oh_scl(phi_scl, T)


# plt.plot(c_oh_scl)
# plt.plot(c_oh_scl_2)
# plt.yscale("log")
# plt.show()


# c_oh_scl = get_c_oh_scl(phi_0, c_acc_bulk, T, x_step_length_approx)
#
# for i in range(0, 3):
#     print("sum(c_oh_scl)", sum(c_oh_scl))
#     c_acc_scl = get_c_acc_scl(c_acc_no_trapping, c_oh_scl, K_trapping)
#     c_oh_scl = get_c_oh_scl_from_c_acc_scl(c_acc_no_trapping, c_acc_scl, K_trapping)
