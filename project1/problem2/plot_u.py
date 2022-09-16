# from turtle import color
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'figure.autolayout': True})

x, y = np.loadtxt("x_u_exact.txt",usecols=(0,1),unpack=True)

plt.figure(figsize=(4,3))
plt.plot(x, y, color='k', linewidth=1.0)
plt.xlabel(r'$x$')
plt.ylabel(r'$u(x)$')
plt.savefig("../imgs/u_exact.pdf")

