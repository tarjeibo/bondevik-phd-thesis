#Script to fit segregation energies, for 210 grain boundary
from math import sqrt
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

a_0 = 4.2e-10
interplanar_distance = a_0/(2*sqrt(5))

e_seg_c_v_input = [0.00, 0.16, 0.01, 0.15, -0.04, -0.14, -0.08, 0.08, -0.16, 0.10, -1.41, -1.02, -1.02, 0.63, -1.41, -1.83, -0.42, -0.52, 1.79, -0.33, -0.11, 0.05, -0.29, -0.07, 0.03, 0.30, 0.10, 0.21, 0.09, 0.01, 0.00, 0.16]
e_seg_c_oh_input = [0.09, 0.00, -0.03, 0.12, -0.03, -0.07, -0.02, 0.12, -0.18, 0.14, 0.08, -0.30, -1.30, 0.52, -0.71, -0.71, -0.77, 0.13, 0.21, -0.12, -0.13, 0.10, -0.10, -0.08, 0.01, 0.16, -0.01, 0.10, 0.00, -0.03, 0.09, 0.0]

e_seg_c_v_input = np.asarray(e_seg_c_v_input)
e_seg_c_oh_input = np.asarray(e_seg_c_oh_input)


def get_x_values_vector():
    x_values = np.zeros(len(e_seg_c_v_input))
    current_x = -11*interplanar_distance
    for i in range(0, len(e_seg_c_v_input)):
        if (i-1)%3!=0:
            current_x = current_x + interplanar_distance
        x_values[i] = current_x
    return x_values


def get_gaussian_params(x, y):
    mean = np.sum(x*y)/len(x)
    sigma = 5e-10
    p0 = np.asarray([-1, mean, sigma])
    popt, pcov = curve_fit(gaussian, x, y, p0=[-1, mean, sigma])
    return popt

def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def get_popt_c_oh():
    return popt_c_oh

def get_popt_c_v():
    return popt_c_v

def set_modified_gaussian_params(popt_c_oh, std_dev_multiplier, x_for_plot, T):
    # T = 423
    k_B_eV = 8.6173325e-5

    # initial_integral = np.trapz(gaussian(x_for_plot, *popt_c_oh)**3, x_for_plot)
    # current_integral = np.trapz(gaussian(x_for_plot, *popt_c_oh)**3, x_for_plot)

    initial_integral = np.trapz(np.exp(-gaussian(x_for_plot, *popt_c_oh)/(k_B_eV*T)), x_for_plot)

    popt_c_oh[2] = popt_c_oh[2]*std_dev_multiplier #Use double std.deviation

    current_integral = np.trapz(np.exp(-gaussian(x_for_plot, *popt_c_oh)/(k_B_eV*T)), x_for_plot)

    # print(initial_integral)
    # print(popt_c_oh[0])
    # print(current_integral)

    while current_integral > initial_integral:
        popt_c_oh[0] = popt_c_oh[0] + 1e-3
        current_integral = np.trapz(np.exp(-gaussian(x_for_plot, *popt_c_oh)/(k_B_eV*T)), x_for_plot)


def get_fit_params(T):

    x_values = get_x_values_vector()
    use_only_negatives_energies = False
    x_for_plot = np.linspace(-10*interplanar_distance, 10*interplanar_distance, 1000)

    if use_only_negatives_energies:
        negative_idxs_v = np.argwhere(e_seg_c_v_input < 0)[:,0]
        x_values_v = x_values[negative_idxs_v]
        e_seg_c_v_negative = e_seg_c_v_input[negative_idxs_v]

        negative_idxs_oh = np.argwhere(e_seg_c_oh_input < 0)[:,0]
        x_values_oh = x_values[negative_idxs_oh]
        e_seg_c_oh_negative = e_seg_c_oh_input[negative_idxs_oh]

        popt_c_v = get_gaussian_params(x_values_v, e_seg_c_v_negative)
        popt_c_oh = get_gaussian_params(x_values_oh, e_seg_c_oh_negative)

        '''Plot fit'''
        # e_seg_c_oh_for_plot = gaussian(x_for_plot, *popt_c_oh)
        # plt.scatter(x_values_oh, e_seg_c_oh_negative, label="post")
        # plt.xlim(-10.5*interplanar_distance, 10.5*interplanar_distance)
        # plt.plot(x_for_plot, e_seg_c_oh_for_plot)
        # plt.legend(loc="best")
        # plt.show()

    else:
        popt_c_v = get_gaussian_params(x_values, e_seg_c_v_input)
        popt_c_oh = get_gaussian_params(x_values, e_seg_c_oh_input)

        std_dev_multiplier = 2
        set_modified_gaussian_params(popt_c_oh, std_dev_multiplier, x_for_plot, T)



        # '''Plot fit'''
        # e_seg_c_oh_for_plot = gaussian(x_for_plot, *popt_c_oh)
        # plt.scatter(x_values, e_seg_c_oh_input, label="post")
        # plt.xlim(-10.5*interplanar_distance, 10.5*interplanar_distance)
        # plt.plot(x_for_plot, e_seg_c_oh_for_plot)
        # plt.legend(loc="best")
        #
        # # '''Save plots'''
        # # np.savetxt("x_values_input.csv", x_values)
        # np.savetxt("e_seg_c_oh_input.csv", e_seg_c_oh_input)
        # np.savetxt("x_for_plot_fit.csv", x_for_plot)
        # np.savetxt("e_seg_c_oh_for_plot.csv", e_seg_c_oh_for_plot)
        #
        # plt.show()



    return popt_c_v, popt_c_oh
#
#
# '''Save plots'''
# np.savetxt("x_values_input.csv", x_values)
# np.savetxt("e_seg_c_oh_input.csv", e_seg_c_oh_input)
# np.savetxt("x_for_plot_fit.csv", x_for_plot)
# np.savetxt("e_seg_c_oh_for_plot.csv", e_seg_c_oh_for_plot)
