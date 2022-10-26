from ast import For
from cProfile import label
from sre_parse import fix_flags
import numpy as np
import matplotlib.pyplot as plt


# Plot one particle in xy-plane
def plot_single_xy(filename):
    x, y = np.loadtxt("data/" + filename + ".txt", usecols=(1,3), unpack=True)

    fig_aspect = 1
    fig_width = 3.5
    fig_height = fig_width * fig_aspect

    xy_lim = 60

    fig, ax = plt.subplots()

    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)

    ax.set_xlim(-xy_lim, xy_lim)
    ax.set_ylim(-xy_lim, xy_lim)
    ax.set_aspect('equal')

    ax.plot(x, y, 'k', lw=1.0)

    ax.set_xlabel(r'$x$ (\textmu m)')
    ax.set_ylabel(r'$y$ (\textmu m)')

    fig.savefig("imgs/" + filename + "_xy.pdf")



# Plot z(t) for one particle
def plot_single_tz(filename):
    t, z = np.loadtxt("data/" + filename + ".txt", usecols=(0,5), unpack=True)

    z_height = 25
    t_len = 50

    fig_aspect = 1./2.
    fig_width = 3.5
    fig_height = fig_width * fig_aspect

    fig, ax = plt.subplots()

    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)

    ax.set_xlim(0, t_len)
    ax.set_ylim(-z_height, z_height)

    ax.plot(t, z, 'k', lw=1.0)
    # Mark analytic top 
    # ax.plot(5*2*np.pi*np.sqrt(40*500*500/(2*2.41e6)), 20, 'ro', ms=2.0)

    ax.set_xlabel(r'$t$ (\textmu s)')
    ax.set_ylabel(r'$z$ (\textmu m)')

    ax.grid()
    fig.savefig("imgs/" + filename + "_tz.pdf")


# Plot two particles in xy-plane
def plot_two_xy(filename):
    x1, y1 = np.loadtxt("data/" + filename + "1.txt", usecols=(1,3), unpack=True)
    x2, y2 = np.loadtxt("data/" + filename + "2.txt", usecols=(1,3), unpack=True) 

    xy_width = 80

    fig_aspect = 1.
    fig_width = 3.5
    fig_height = fig_width * fig_aspect

    fig, ax = plt.subplots()

    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)

    ax.set_xlim(-xy_width, xy_width)
    ax.set_ylim(-xy_width, xy_width)
    ax.set_aspect('equal')

    ax.plot(x1, y1, 'k', lw=1.)
    ax.plot(x2, y2, 'r', lw=1.)

    ax.set_xlabel(r'$x$ (\textmu m)')
    ax.set_ylabel(r'$y$ (\textmu m)')

    fig.savefig("imgs/"+ filename + "_xy.pdf")


def relative_error(r_num, r_an):
    rel_err = np.linalg.norm(r_an - r_num, axis=1)/np.linalg.norm(r_an, axis=1)
    return rel_err


def max_error(r_num, r_an):
    diff = np.linalg.norm(r_an - r_num, axis=1)
    return np.max(diff)


def plot_rel_err(n_steps, filenames):
    fig_aspect = 3./2.
    fig_width = 3.5
    fig_height = fig_width*fig_aspect

    fig, ax = plt.subplots(2, 1, figsize=(fig_width,fig_height))

    for j in range(2):
        max_err = np.zeros(len(n_steps))
        h = np.zeros(len(n_steps))

        for i in range(len(n_steps)):
            t = np.loadtxt("data/" + filenames[j] + "_" + str(n_steps[i]) 
                            + "_1.txt", usecols=(0))
            r_num = np.loadtxt("data/" + filenames[j] + "_" + str(n_steps[i]) 
                                + "_1.txt", usecols=(1,3,5))
            r_an = np.loadtxt("data/" + filenames[2] + str(n_steps[i]) + ".txt", 
                                usecols=(1,3,5))

            rel_err = relative_error(r_num, r_an)

            max_err[i] = max_error(r_num, r_an)
            h[i] = t[-1]/n_steps[i]

            ax[j].semilogy(t, rel_err, label=f'$n = {n_steps[i]}$')

        ax[j].set_xlabel(r'$t$ (\textmu s)')
        ax[j].set_ylabel(r'Relative error')
        ax[j].legend()

        converg_rate = 0.

        for i in range(3):
            converg_rate += np.log(max_err[i+1]/max_err[i])/(3 * 
                                    np.log(h[i+1]/h[i]))

        print("Convergence rate: ", converg_rate)

    fig.savefig("imgs/" + filenames[0] + "_Euler.pdf")

    

def plot_two_xvx(filename):
    x1, vx1 = np.loadtxt("data/" + filename + "1.txt", usecols=(1,2), unpack=True)
    x2, vx2 = np.loadtxt("data/" + filename + "2.txt", usecols=(1,2), unpack=True)

    xy_width = 80
    v_width = 80

    plt.figure(figsize=(4,3))
    plt.xlim(-xy_width, xy_width)
    plt.ylim(-v_width, v_width)
    plt.plot(x1, vx1, 'k', lw=1.)
    plt.plot(x2, vx2, 'r', lw=1.)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$v_x$')
    plt.savefig("imgs/"+ filename + "_xvx.pdf")


def plot_two_rvr(filenames):
    xy_width = 80
    v_low = 20
    v_high = 55

    fig_aspect = 3./2.
    fig_width = 3.5
    fig_height = fig_width * fig_aspect

    fig, ax = plt.subplots(2, 1, figsize=(fig_width, fig_height))

    for i in range(2):
        x1, vx1, y1, vy1 = np.loadtxt("data/" + filenames[i] + "1.txt", 
                                        usecols=(1,2,3,4), unpack=True)
        x2, vx2, y2, vy2 = np.loadtxt("data/" + filenames[i] + "2.txt", 
                                        usecols=(1,2,3,4), unpack=True)

        r1 = np.sqrt(x1**2 + y1**2)
        vr1 = np.sqrt(vx1**2 + vy1**2)

        r2 = np.sqrt(x2**2 + y2**2)
        vr2 = np.sqrt(vx2**2 + vy2**2)


        ax[i].set_xlim(0, xy_width)
        ax[i].set_ylim(v_low, v_high)

        ax[i].plot(r1, vr1, 'k', label=r'$p_1$')
        ax[i].plot(r2, vr2, 'r', label=r'$p_2$')

        ax[i].set_xlabel(r'$\rho$ (\textmu m)')
        ax[i].set_ylabel(r'$v_\rho$ (\textmu m/\textmu s)')

        ax[i].legend()

    fig.savefig("imgs/"+ filenames[0] + "_rvr.pdf")


def plot_two_zvz(filenames):
    z_height = 30
    v_width = 20

    fig_aspect = 3./2.
    fig_width = 3.5
    fig_height = fig_width * fig_aspect

    fig, ax = plt.subplots(2, 1, figsize=(fig_width, fig_height))

    for i in range(2):
        z1, vz1 = np.loadtxt("data/" + filenames[i] + "1.txt", usecols=(5,6), 
                                unpack=True)
        z2, vz2 = np.loadtxt("data/" + filenames[i] + "2.txt", usecols=(5,6), 
                                unpack=True) 

        ax[i].set_xlim(-z_height, z_height)
        ax[i].set_ylim(-v_width, v_width)

        ax[i].plot(z1, vz1, 'k', label=r'$p_1$')
        ax[i].plot(z2, vz2, 'r', label=r'$p_2$')

        ax[i].set_xlabel(r'$z$ (\textmu m)')
        ax[i].set_ylabel(r'$v_z$ (\textmu m/\textmu s)')

        ax[i].legend()

    fig.savefig("imgs/"+ filenames[0] + "_zvz.pdf")



def plot_two_xy(filenames):
    xy_width = 80

    fig_aspect = 1.9
    fig_width = 3.5
    fig_height = fig_width * fig_aspect

    fig, ax = plt.subplots(2, 1, figsize=(fig_width, fig_height))

    for i in range(2):
        x1, y1 = np.loadtxt("data/" + filenames[i] + "1.txt", usecols=(1,3), unpack=True)
        x2, y2 = np.loadtxt("data/" + filenames[i] + "2.txt", usecols=(1,3), unpack=True) 

        ax[i].set_xlim(-xy_width, xy_width)
        ax[i].set_ylim(-xy_width, xy_width)
        ax[i].set_aspect('equal')

        ax[i].plot(x1, y1, 'k', label=r'$p_1$')
        ax[i].plot(x2, y2, 'r', label=r'$p_2$')

        ax[i].set_xlabel(r'$x$ (\textmu m)')
        ax[i].set_ylabel(r'$y$ (\textmu m)')

        ax[i].legend()


    fig.savefig("imgs/"+ filenames[0] + "_xy.pdf")


def plot_xyz(filename):
    # fig_aspect = 3./4.
    # fig_width = 3.5
    # fig_height = fig_width * fig_aspect

    # fig = plt.figure(figsize=(fig_width, fig_height))
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    x, y, z = np.loadtxt("data/" + filename + "1.txt", usecols=(1,3,5), unpack=True)

    ax.plot(x, y, z, 'k', label=r'$p_1$')

    ax.set_xlabel(r'$x$ (\textmu m)')
    ax.set_ylabel(r'$y$ (\textmu m)')
    ax.set_zlabel(r'$z$ (\textmu m)')

    ax.legend()

    fig.savefig("imgs/" + filename + "_xyz.pdf")


def plot_two_xyz(filenames):
    fig_aspect = 1.9
    fig_width = 3.5
    fig_height = fig_width * fig_aspect

    fig, ax = plt.subplots(2, 1, subplot_kw={'projection': '3d'}, 
                            figsize=(fig_width, fig_height))

    for i in range(2):
        x1, y1, z1 = np.loadtxt("data/" + filenames[i] + "1.txt", usecols=(1,3,5), 
                                    unpack=True)
        x2, y2, z2 = np.loadtxt("data/" + filenames[i] + "2.txt", usecols=(1,3,5), 
                                    unpack=True) 

        # ax[i].set_xlim(-xy_width, xy_width)
        # ax[i].set_ylim(-xy_width, xy_width)
        # ax[i].set_aspect('equal')

        ax[i].plot(x1, y1, z1, 'k', label=r'$p_1$')
        ax[i].plot(x2, y2, z2, 'r', label=r'$p_2$')

        # ax[i].set_xlabel(r'$x$ (\textmu m)')
        # ax[i].set_ylabel(r'$y$ (\textmu m)')
        # ax[i].set_zlabel(r'$z$ (\textmu m)')

        ax[i].legend()

    fig.savefig("imgs/" + filenames[0] + "_xyz.pdf")


def plot_resonance(filename):
    omega_V = np.loadtxt("data/" + filename + ".txt", usecols=(0))
    p = np.loadtxt("data/" + filename + ".txt", usecols=(1,2,3))

    f = np.array([.1, .4, .7])

    fig, ax = plt.subplots()

    # Fix range
    for i in range(3):
        ax.plot(omega_V, p[:,i], label=f'$f = {f[i]}$')
    
    ax.set_xlabel(r'$\omega_V$ (MHz)')
    ax.set_ylabel(r'$N/N_0$')
    # ax.set_ylabel(r'Fraction of remaining particles')

    ax.grid()
    ax.legend()

    fig.savefig("imgs/" + filename + ".pdf")


def main():

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "Computer Modern Roman"
    plt.rcParams.update({'lines.linewidth': 1., 'axes.grid': True, 'grid.linewidth': 0.5})
    plt.rcParams.update({'figure.autolayout': True})

    # plot_single_xy("singlepart1")

    # # First plot
    # plot_single_tz("singlepart1")

    # plot_two_xy("twoparts_int")
    # plot_two_xy("twoparts_noint")

    # # Second plot
    # plot_two_xy(["twoparts_noint", "twoparts_int"])

    # plot_two_xvx("twoparts_noint")
    # plot_two_xvx("twoparts_int")

    # # Third plot
    # plot_two_zvz(["twoparts_noint", "twoparts_int"])

    # # Forth plot
    # plot_two_rvr(["twoparts_noint", "twoparts_int"])

    # # Fifth plot
    # plot_two_xyz(["twoparts_noint", "twoparts_int"])

    # plot_xyz("twoparts_noint")

    # plot_two_xy_int("twoparts_noint", "twoparts_int")

    # plot_single_xy("exact_4000")

    # n_steps = np.array([2000, 5000, 8000, 10000, 20000, 40000, 50000])
    n_steps = np.array([4000, 8000, 16000, 32000])
    plot_rel_err(n_steps, ["errorRK4", "errorEuler", "exact_"])

    # plot_resonance("resonance_broad_500_5000")

    # plot_resonance("resonance_analysis_50_500_0")

    return 0


main()