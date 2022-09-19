import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'figure.autolayout': True})

# Plot for problem 2
x, u = np.loadtxt("data/u_1000.txt",usecols=(0,1),unpack=True)

plt.figure(figsize=(4,3))
plt.plot(x, u, color='k', linewidth=1.0)
plt.xlabel(r'$x$')
plt.ylabel(r'$u(x)$')
plt.savefig("imgs/u_exact.pdf")


# Plot for problem 7
plt.figure(figsize=(4,3))
plt.plot(x, u, linewidth=1., label=r'$u(x)$')

x, v = np.loadtxt("data/x_v_10.txt", usecols=(0,1), unpack=True)
plt.plot(x, v, linewidth=1., linestyle='-', label=r'$v(x)$, $n_{steps}=10$')
x, v = np.loadtxt("data/x_v_100.txt", usecols=(0,1), unpack=True)
plt.plot(x, v, linewidth=1., linestyle='--', label=r'$v(x)$, $n_{steps}=100$')
x, v = np.loadtxt("data/x_v_1000.txt", usecols=(0,1), unpack=True)
plt.plot(x, v, linewidth=1., linestyle=':', label=r'$v(x)$, $n_{steps}=1000$')



plt.xlabel(r'$x$')
plt.ylabel(r'$v(x)$')
plt.legend()
plt.savefig("imgs/v.pdf")
