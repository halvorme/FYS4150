import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'figure.autolayout': True})

# Set maximal number of steps on the x-axis (10^{max_steps+1})
# Must be lower than max_steps set in 'main.cpp'
max_steps7 = 3
max_steps8 = 5


# Set proportions of figures
width = 5
height = 3/4*width


# Plot for problem 2
x, u = np.loadtxt("data/u_1000.txt",usecols=(0,1),unpack=True)

plt.figure(1, figsize=(4,3))
plt.plot(x, u, color='k', linewidth=1.0)
plt.xlabel(r'$x$')
plt.ylabel(r'$u(x)$')
plt.savefig("imgs/u_exact.pdf")


# Plot for problem 7

plt.figure(2, figsize=(4,3))
plt.plot(x, u, linewidth=1., label=r'$u(x)$')

styles = ['-','--',':']

for i in range(max_steps7):
	x, v = np.loadtxt("data/x_v_" + str(i+1) + ".txt", usecols=(0,1), unpack=True)
	plt.plot(x, v, linewidth=1., linestyle=styles[i], label=f'$v(x)$, $n=10^{i+1}$')

plt.xlabel(r'$x$')
plt.ylabel(r'$v(x)$')
plt.legend()
plt.savefig("imgs/v.pdf")


# Plot for problem 8
plt.rcParams.update({'legend.loc': 'upper center'})


plt.figure(3, figsize=(width,height))
plt.figure(4, figsize=(width,height))

for i in range(max_steps8):
	x, ab, rel = np.loadtxt("data/err_" + str(i+1) + ".txt", usecols=(0,1,2), unpack=True)
	plt.figure(3)
	plt.plot(x[1:-1], ab[1:-1], label=f'$n=10^{i+1}$')

	plt.figure(4)
	plt.plot(x[1:-1], rel[1:-1], label=f'$n=10^{i+1}$')

plt.figure(3)
plt.xlabel(r"$x$")
plt.ylabel(r'$\log_{10} \vert u-v\vert$')
plt.legend()
plt.savefig("imgs/abs_err.pdf")

plt.figure(4)
plt.xlabel(r"$x$")
plt.ylabel(r'$\log_{10} \vert (u-v)/u\vert$')
plt.legend()
plt.savefig("imgs/rel_err.pdf")



