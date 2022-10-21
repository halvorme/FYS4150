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


def plot_two_zvz(filename):
    z1, vz1 = np.loadtxt("data/" + filename + "1.txt", usecols=(5,6), unpack=True)
    z2, vz2 = np.loadtxt("data/" + filename + "2.txt", usecols=(5,6), unpack=True)

    height = 30
    v_width = 20

    plt.figure(figsize=(4,3))
    plt.xlim(-height, height)
    plt.ylim(-v_width, v_width)
    plt.plot(z1, vz1, 'k', lw=1.)
    plt.plot(z2, vz2, 'r', lw=1.)

    plt.plot(5*np.sqrt(2),0, 'bo')

    plt.xlabel(r'$z$')
    plt.ylabel(r'$v_z$')
    plt.savefig("imgs/"+ filename + "_zvz.pdf")


def plot_two_rvr(filename):
    x1, vx1, y1, vy1 = np.loadtxt("data/" + filename + "1.txt", usecols=(1,2,3,4), unpack=True)
    x2, vx2, y2, vy2 = np.loadtxt("data/" + filename + "2.txt", usecols=(1,2,3,4), unpack=True)

    r1 = np.sqrt(x1**2 + y1**2)
    vr1 = np.sqrt(vx1**2 + vy1**2)

    r2 = np.sqrt(x2**2 + y2**2)
    vr2 = np.sqrt(vx2**2 + vy2**2)

    xy_width = 80
    v_width = 55

    plt.figure(figsize=(4,3))
    plt.xlim(0, xy_width)
    plt.ylim(20, v_width)
    plt.plot(r1, vr1, 'k', lw=1.)
    plt.plot(r2, vr2, 'r', lw=1.)
    plt.xlabel(r'$r$')
    plt.ylabel(r'$v_r$')
    plt.savefig("imgs/"+ filename + "_rvr.pdf")


def plot_two_xy_int(filename1, filename2):
    x1, y1 = np.loadtxt("data/" + filename1 + "1.txt", usecols=(1,3), unpack=True)
    x2, y2 = np.loadtxt("data/" + filename1 + "2.txt", usecols=(1,3), unpack=True) 

    xy_width = 80

    plt.figure(figsize=(3.5,7))


    plt.subplot(211)
    plt.xlim(-xy_width, xy_width)
    plt.ylim(-xy_width, xy_width)
    plt.plot(x1, y1, 'k', lw=1.)
    plt.plot(x2, y2, 'r', lw=1.)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')

    x1, y1 = np.loadtxt("data/" + filename2 + "1.txt", usecols=(1,3), unpack=True)
    x2, y2 = np.loadtxt("data/" + filename2 + "2.txt", usecols=(1,3), unpack=True) 

    plt.subplot(212)
    plt.xlim(-xy_width, xy_width)
    plt.ylim(-xy_width, xy_width)
    plt.plot(x1, y1, 'k', lw=1.)
    plt.plot(x2, y2, 'r', lw=1.)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')

    plt.savefig("imgs/"+ filename1 + "_xy.pdf")


def plot_two_xy_int_new(filenames):
    xy_width = 80

    fig_aspect = 2.
    fig_width = 3.5
    fig_height = fig_width * fig_aspect

    fig, ax = plt.subplots(2, 1, figsize=(fig_width, fig_height))
    for i in range(2):
        x1, y1 = np.loadtxt("data/" + filenames[i] + "1.txt", usecols=(1,3), unpack=True)
        x2, y2 = np.loadtxt("data/" + filenames[i] + "2.txt", usecols=(1,3), unpack=True) 

        ax[i].set_xlim(-xy_width, xy_width)
        ax[i].set_ylim(-xy_width, xy_width)
        ax[i].set_aspect('equal')

        ax[i].plot(x1, y1, 'k', lw=1.)
        ax[i].plot(x2, y2, 'r', lw=1.)

        ax[i].set_xlabel(r'$x$ (\textmu m)')
        ax[i].set_ylabel(r'$y$ (\textmu m)')

    fig.savefig("imgs/"+ filenames[0] + "_xy_new.pdf")

def main():

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "Computer Modern Roman"
    plt.rcParams.update({'figure.autolayout': True})

    # plot_single_xy("singlepart1")

    # # First plot
    # plot_single_tz("singlepart1")

    # plot_two_xy("twoparts_int")
    # plot_two_xy("twoparts_noint")

    # plot_two_xy_int_new(["twoparts_noint", "twoparts_int"])

    plot_two_xvx("twoparts_noint")
    plot_two_xvx("twoparts_int")
    plot_two_zvz("twoparts_noint")
    plot_two_zvz("twoparts_int")
    plot_two_rvr("twoparts_noint")
    plot_two_rvr("twoparts_int")

    # plot_two_xy_int("twoparts_noint", "twoparts_int")

    return 0


main()