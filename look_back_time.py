# -*- coding: utf-8 -*-
#@author: Juho-Petteri Lesonen
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

# Defining the function for the lookback time
def g(z):
	return 1 / ((1+z)*H0*Gyr*sqrt(o_r*(1+z)**4+o_m*(1+z)**3+o_l+(1-o_z)*(1+z)**2))

# Initial values
a0 = 1E-10 # start value of z
af = 25 # end value of z
da = 0.01 # delta_z
steps = int((af-a0)/da)
o_r = 0*9.2364E-5 # radiation density
o_m = 0.3089 # matter density
o_l = 0.6911 # dark energy density
o_z = o_r + o_m + o_l # curvature term
H = 67.74 # H-L parameter
Mpc = 3.085677581E19 #
H0 = (H/Mpc)
Gyr = 3.1536E16
t = np.zeros(steps,float)
z = np.zeros(steps,float)

# Calculating the integration points and weighs as 
# Gauss-Legendre approximation of Nth order
N = 100
x,w = np.polynomial.legendre.leggauss(N)

# Integration for every z
b = a0
for i in range(steps):
    s = 0.0
    xp = 0.5*(b-a0)*x+0.5*(b+a0)
    wp = 0.5*(b+a0)*w
    for k in range(N):
        s += wp[k]*g(xp[k])
    t[i] = s
    z[i] = b
    b += da

# Plotting the graph
plt.plot(z,t, color='blue', linewidth=1.5, label='Benchmark Model')
plt.title('Lookback time(t0-te) as a function of redshift(z)\n Four different')
plt.legend(loc=4)
plt.axis([0,20,0,17.5])
plt.xlabel('z')
plt.ylabel('t0-te (Gyr)')
plt.grid(True)
plt.show()
