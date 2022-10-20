import numpy as np
import matplotlib.pyplot as plt


# Plot one particle in xy-plane
def plot_single_xy(filename):
    x, y = np.loadtxt("data/" + filename + ".txt", usecols=(1,3), unpack=True)

    xy_width = 60

    plt.figure(figsize=(4,3))
    plt.xlim(-xy_width, xy_width)
    plt.ylim(-xy_width, xy_width)
    plt.plot(x, y, 'k', linewidth=1.0)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.savefig("imgs/" + filename + "_xy.pdf")


# Plot z(t) for one particle
def plot_single_tz(filename):
    t, z = np.loadtxt("data/" + filename + ".txt", usecols=(0,5), unpack=True)

    z_height = 25
    t_len = 50

    plt.figure(figsize=(4,3))
    plt.xlim(0, t_len)
    plt.ylim(-z_height, z_height)
    plt.plot(t, z, 'k', lw=1.0)
    # Mark analytic top 
    plt.plot(5*2*np.pi*np.sqrt(40*500*500/(2*2.41e6)), 20, 'ro')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$z(t)$')
    plt.savefig("imgs/" + filename + "_tz.pdf")


# Plot two particles in xy-plane
def plot_two_xy(filename):
    x1, y1 = np.loadtxt("data/" + filename + "1.txt", usecols=(1,3), unpack=True)
    x2, y2 = np.loadtxt("data/" + filename + "2.txt", usecols=(1,3), unpack=True) 

    xy_width = 80

    plt.figure(figsize=(4,3))
    plt.xlim(-xy_width, xy_width)
    plt.ylim(-xy_width, xy_width)
    plt.plot(x1, y1, 'k', lw=1.)
    plt.plot(x2, y2, 'r', lw=1.)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.savefig("imgs/"+ filename + "_xy.pdf")



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

def main():

    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'figure.autolayout': True})

    plot_single_xy("singlepart1")
    plot_single_tz("singlepart1")
    plot_two_xy("twoparts_int")
    plot_two_xy("twoparts_noint")
    plot_two_xvx("twoparts_noint")
    plot_two_xvx("twoparts_int")
    plot_two_zvz("twoparts_noint")
    plot_two_zvz("twoparts_int")
    plot_two_rvr("twoparts_noint")
    plot_two_rvr("twoparts_int")

    return 0


main()