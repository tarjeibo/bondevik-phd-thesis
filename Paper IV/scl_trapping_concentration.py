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

#Numerical params
x_step_length_approx = 1e-12


def get_c_acc_bulk_with_trapping(c_acc_no_trapping, K_trapping):
    c_acc = (2*sqrt(3)*c_acc_no_trapping)/(sqrt(4*c_acc_no_trapping*K_trapping*N_conf + 3) + sqrt(3))
    return c_acc.real

def get_trapping_equilibrium_constant(trapping_energy, T):
    #Wolfram alpha: https://www.wolframalpha.com/input/?i=K+%3D+3*(c-x)%2F(N*x%5E2)+for+x
    return np.exp((-e*trapping_energy)/(k_B*T))

def get_scl_length(phi_0, c_acc):
    return np.sqrt((2*epsilon_0*epsilon_r*phi_0)/(e*(c_acc/a_0**3)))

def get_c_oh_scl(phi_0, c_acc, T, x_step_length_approx):
    scl_length = get_scl_length(phi_0, c_acc)
    n_x_steps = int(round(scl_length/x_step_length_approx))
    x_scl_list = np.linspace(0, scl_length, n_x_steps)

    phi_scl = phi_0*(x_scl_list/scl_length - 1)**2
    c_oh_scl = c_acc_no_trapping/np.exp((phi_scl*e)/(k_B*T))
    return c_oh_scl

def get_c_acc_scl(c_acc_no_trapping, c_oh_scl, K_trapping):
    #https://www.wolframalpha.com/input/?i=K+%3D+3*(D-y)%2F(y*h)+for+y
    return 3*c_acc_no_trapping/(c_oh_scl*K_trapping + 3)

def get_c_oh_scl_from_c_acc_scl(c_acc_no_trapping, c_acc_scl, K_trapping):
    return 3*(c_acc_no_trapping - c_acc_scl)/(K_trapping*c_acc_scl)


#Experimental params
epsilon_r = 58
c_acc_no_trapping = 0.20
trapping_energy = 0
T = 423

K_trapping = get_trapping_equilibrium_constant(trapping_energy, T)
c_acc = get_c_acc_bulk_with_trapping(c_acc_no_trapping, K_trapping)

print("c_acc", c_acc)
phi_0 = 0.5

c_oh_scl = get_c_oh_scl(phi_0, c_acc, T, x_step_length_approx)

for i in range(0, 3):
    print("sum(c_oh_scl)", sum(c_oh_scl))
    c_acc_scl = get_c_acc_scl(c_acc_no_trapping, c_oh_scl, K_trapping)
    c_oh_scl = get_c_oh_scl_from_c_acc_scl(c_acc_no_trapping, c_acc_scl, K_trapping)
