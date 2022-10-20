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
    t1, x1, y1, z1 = np.loadtxt("data/" + filename + "1.txt", usecols=(0,1,3,5), unpack=True)
    t2, x2, y2, z2 = np.loadtxt("data/" + filename + "2.txt", usecols=(0,1,3,5), unpack=True) 

    xy_width = 80


    plt.figure(figsize=(4,3))
    plt.xlim(-xy_width, xy_width)
    plt.ylim(-xy_width, xy_width)
    plt.plot(x1, y1, 'k', lw=1.)
    plt.plot(x2, y2, 'r', lw=1.)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.savefig("imgs/"+ filename + "_xy.pdf")

def main():

    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'figure.autolayout': True})

    plot_single_xy("singlepart1")
    plot_single_tz("singlepart1")
    plot_two_xy("twoparts")

main()