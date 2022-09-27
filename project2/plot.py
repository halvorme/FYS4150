import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'figure.autolayout': True})


# Plot for problem 5
n, iter = np.loadtxt("data/tri_iter.txt",usecols=(0,1),unpack=True)

plt.figure(1, figsize=(4,3))
plt.scatter(n, iter, color='k', marker='.')
plt.xlabel(r'$n$')
plt.ylabel(r'$N$')
plt.savefig("imgs/n_iter.pdf")


plt.figure(2, figsize=(4,3))
plt.scatter(n, np.sqrt(iter), color='k', marker='.')
plt.xlabel(r'$n$')
plt.ylabel(r'$\sqrt{N}$')
plt.savefig("imgs/n_sqrt_iter.pdf")


# Plot for problem 6

x, v1, v2, v3 = np.loadtxt("data/x_v10.txt",usecols=(0,1,2,3),unpack=True)
u1, u2, u3 = np.loadtxt("data/x_u10.txt",usecols=(1,2,3),unpack=True)

plt.figure(3, figsize=(4,3))
plt.scatter(x, v1, color='k', marker='.')
plt.plot(x, u1)
plt.xlabel(r'$x$')
plt.ylabel(r'$v_1$')
plt.savefig("imgs/v1_10.pdf")

plt.figure(4, figsize=(4,3))
plt.scatter(x, v2, color='k', marker='.')
plt.plot(x, u2)
plt.xlabel(r'$x$')
plt.ylabel(r'$v_2$')
plt.savefig("imgs/v2_10.pdf")

plt.figure(5, figsize=(4,3))
plt.scatter(x, v3, color='k', marker='.')
plt.plot(x, u3)
plt.xlabel(r'$x$')
plt.ylabel(r'$v_3$')
plt.savefig("imgs/v3_10.pdf")


x, v1, v2, v3 = np.loadtxt("data/x_v100.txt",usecols=(0,1,2,3),unpack=True)
u1, u2, u3 = np.loadtxt("data/x_u100.txt",usecols=(1,2,3),unpack=True)

plt.figure(6, figsize=(4,3))
plt.scatter(x, v1, color='k', marker='.')
plt.plot(x, u1)
plt.xlabel(r'$x$')
plt.ylabel(r'$v_1$')
plt.savefig("imgs/v1_100.pdf")

plt.figure(7, figsize=(4,3))
plt.scatter(x, v2, color='k', marker='.')
plt.plot(x, u2)
plt.xlabel(r'$x$')
plt.ylabel(r'$v_2$')
plt.savefig("imgs/v2_100.pdf")

plt.figure(8, figsize=(4,3))
plt.scatter(x, v3, color='k', marker='.')
plt.plot(x, u3)
plt.xlabel(r'$x$')
plt.ylabel(r'$v_3$')
plt.savefig("imgs/v3_100.pdf")

