import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'figure.autolayout': True})


xy_width = 60
z_height = 25
t_len = 50

# Plot one particle in xy-plane

t, x, y, z = np.loadtxt("data/singlepart.txt",usecols=(0,1,3,5),unpack=True)

plt.figure(1, figsize=(4,3))
plt.xlim(-xy_width, xy_width)
plt.ylim(-xy_width, xy_width)
plt.plot(x, y, 'k', linewidth=1.0)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig("imgs/singlepart_xy.pdf")

# Plot z(t) for one particle
plt.figure(2, figsize=(4,3))
plt.xlim(0, t_len)
plt.ylim(-z_height, z_height)
plt.plot(t, z, 'k', lw=1.0)
# Mark analytic top 
plt.plot(5*2*np.pi*np.sqrt(40*500*500/(2*2.41e6)), 20, 'ro')
plt.xlabel(r'$t$')
plt.ylabel(r'$z(t)$')
plt.savefig("imgs/singlepart_tz.pdf")

# Plot two particles in xy-plane

t1, x1, y1, z1 = np.loadtxt("data/twoparts1.txt", usecols=(0,1,3,5), unpack=True)
t2, x2, y2, z2 = np.loadtxt("data/twoparts2.txt", usecols=(0,1,3,5), unpack=True) 

plt.figure(3, figsize=(4,3))
#plt.xlim(-xy_width, xy_width)
#plt.ylim(-xy_width, xy_width)
plt.plot(x1, y1, 'k', lw=1.)
plt.plot(x2, y2, 'r', lw=1.)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig("imgs/twoparts_xy.pdf")

