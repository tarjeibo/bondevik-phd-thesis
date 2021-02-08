from functions import *
from librariesAndConstants import *

#Variable input params
d_core = 0.5e-9
u_bulk_gb_ratio = 1


#Constant input params
T = 373
phi_0 = 0.37
rho_bulk = 1
c_bulk = c_acc
n_steps = 200
w_list = np.logspace(-4, 12, 2000)

#SCL impedance
scl_length = get_scl_length(phi_0)
x_list = np.linspace(0, scl_length, n_steps)
x_step = scl_length/n_steps
phi_list = get_phi_list(x_list, phi_0, scl_length)
rho_scl = get_rho_scl(rho_bulk, phi_list, T)
z_scl = get_scl_impedance(rho_scl, x_step, w_list)

# print("scl_length", scl_length)
# print("tot resistance", np.sum(z_scl.real))

# print(get_w0_ratio(z_scl, w_list, rho_bulk))


d_cores = np.linspace(1e-10, 40e-10, 70)
u_bulk_gb_ratios = np.logspace(3, 7, 70)

w0_ratios = np.zeros(shape=(len(u_bulk_gb_ratios), len(d_cores)))

for i, d_core in enumerate(d_cores):
    for j, u_bulk_gb_ratio in enumerate(u_bulk_gb_ratios):

        #GB core impedance
        c_core = get_c_core(d_core, scl_length)
        rho_core = get_rho_core(rho_bulk, c_bulk, c_core, u_bulk_gb_ratio)
        z_core = get_core_impedance(d_core, rho_core, w_list)

        #Total impedance
        z_gb = z_scl + z_core

        # plt.plot(w_list, -z_scl.imag, label='scl')
        # plt.plot(w_list, -z_core.imag, label='core')
        # plt.plot(w_list, -z_gb.imag, label='gb')
        # plt.legend(loc="best")
        # plt.xscale("log")
        # plt.xlim(1e1, 1e9)
        # plt.show()

        w0_ratios[j][i] = get_w0_ratio(z_gb, w_list, rho_bulk)

        # print(get_w0_ratio(z_gb, w_list, rho_bulk))

        # print(get_w0_ratio(z_gb, w_list, rho_bulk))
#
#
#
#

plot_data(d_cores, u_bulk_gb_ratios, w0_ratios)
plt.savefig("foo.png", dpi=350)
plt.show()



# print(get_w0_ratio(z_gb, w_list, rho_bulk))

# print(w0_max_gb)
# print(w0_max_bulk)
# print(w0_max_gb/w0_max_bulk)









#
