import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'figure.autolayout': True})

x, y = np.loadtxt("data/singlepart.txt",usecols=(1,3),unpack=True)

plt.figure(1, figsize=(4,3))
plt.plot(x, y, color='k', linewidth=1.0)
plt.xlabel(r'$x$')
plt.ylabel(r'$u(x)$')
plt.savefig("imgs/singlepart.pdf")
