# -*- coding: utf-8 -*-
#@author: Juho-Petteri Lesonen
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

# Defining the function for the evolution
def f(a):
    return 1 / (sqrt(o_r/a**2 + o_m/a + o_l*a**2 + (1 - o_z)))
    
# Initial values
a0 = 1E-17 # Start value of a
af = 1000 # End value of a
da = 0.01 # delta_a
steps = int((af - a0) / da)
o_r = 9.2364E-5 # radiation density
o_m = 0.3089 # matter density
o_l = 0.6911 # dark energy denisty
o_z = o_r + o_m + o_l # curvature term

t = np.zeros(steps,float)
a = np.zeros(steps,float)

# Calculating the integration points and weights as 
# Gauss-Legendre approximation of Nth order
N = 50
x,w = np.polynomial.legendre.leggauss(N)

# Integration for every a
b = a0
for i in range(steps):
    s = 0.0
    xp = 0.5*(b - a0)*x + 0.5*(b + a0)
    wp = 0.5*(b + a0)*w
    for k in range(N):
        s += wp[k]*f(xp[k])
    t[i] = s
    a[i] = b
    b += da

# Turning to logarithmic scale
t=np.log10(t)
a=np.log10(a)
# Plotting the graph
plt.plot(t,a, color='blue', linewidth=1.5)
plt.title('Evolution of the scale-factor(a) as a function of time(H0t) \n
           Benchmark model')
plt.axis([-10,1,-6,3])
plt.xlabel('H0t')
plt.ylabel('a')
plt.grid(True)
plt.show()
