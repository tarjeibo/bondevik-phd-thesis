from librariesAndConstants import *

def get_scl_length(phi_0):
    return sqrt((2*epsilon_0*epsilon_r*phi_0)/(e*(c_acc/a_0**3)))

def get_c_core(d_core, scl_length):
    q_scl = 2*scl_length*c_acc
    return q_scl/d_core

def get_phi_list(x_list, phi_0, scl_length):
    return phi_0*(x_list/scl_length - 1)**2

def get_rho_scl(rho_bulk, phi_list, T):
    return rho_bulk*np.exp((e*phi_list)/(k_B*T))

def get_scl_impedance(rho_scl, x_step, w_list):
    z_scl = []
    for w in w_list:
        z_scl_at_w = np.sum((rho_scl*x_step)/(1 + j*w*rho_scl*epsilon_0*epsilon_r))
        z_scl.append(z_scl_at_w)
    return np.asarray(z_scl)

def get_rho_core(rho_bulk, c_bulk, c_core, u_bulk_gb_ratio):
    return rho_bulk*(c_bulk/c_core)*u_bulk_gb_ratio

def get_core_impedance(d_core, rho_core, w_list):
    z_core = []
    for w in w_list:
        z_core_at_w = (d_core*rho_core)/(1 + j*w*rho_core*epsilon_0*epsilon_r)
        z_core.append(z_core_at_w)
    return np.asarray(z_core)

def get_w0_max(z, w_list):
    return w_list[np.argmin(z.imag)]

def get_w0_ratio(z_gb, w_list, rho_bulk):
    w0_max_gb = get_w0_max(z_gb, w_list)
    w0_max_bulk = (1/(rho_bulk*epsilon_0*epsilon_r))
    return w0_max_gb/w0_max_bulk




def plot_data(x, y, data):
    fontsize=16

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # cax = ax.matshow(data, interpolation='nearest', vmin=0, vmax=1)

    extent = [min(x)*1e10, max(x)*1e10, log10(max(y)).real, log10(min(y)).real]
    aspect = max(x)*1e10/(log10(max(y)).real + 1)
    print(aspect)
    cax = ax.matshow(data, interpolation='nearest', norm=LogNorm(), extent=extent, aspect=aspect)

    fig.colorbar(cax)

    # put text on matrix elements
    for i, x_val in enumerate(np.arange(len(x))):
        for j, y_val in enumerate(np.arange(len(y))):
            # c = "${0:.1f}\\%$".format( 100*data[j,i])
            c = "${0:.1f}$".format(data[j,i])

            # c = data[j, i]
            c = '%.2E' % Decimal(data[j,i])
            #'%.2E' % Decimal('40800000000.00000000000000')

            # ax.text(x_val, y_val, c, va='center', ha='center')

    ax.set_xlabel('$\\mathrm{\delta_{core}}$ / Å',fontsize=fontsize)
    ax.set_ylabel('log $\\mathrm{\mu_{bulk} / \mu_{gb}}$',fontsize=fontsize)

    plt.tight_layout()
