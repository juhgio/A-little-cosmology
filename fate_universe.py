# -*- coding: utf-8 -*-
#@author: Juho-Petteri Lesonen
import numpy as np
import matplotlib.pylab as plt
from sympy.solvers import solve
from sympy import Symbol

# Initial lists
O_l = (list(np.arange(-1,3.125,0.125)))*33
O_m = []
a = 0
da = 0.125
for i in range(33):
    for k in range(33):
        O_m.append(a)
    a = a + da

# Solving all roots for every point of the grid
for i in range(len(O_l)):
    O_l1=O_l[i]
    O_m1=O_m[i]
    a = Symbol('a', real=True)
    sc = solve((O_m1/a**3) + O_l1 + ((1 - (O_m1 + O_l1))/a**2), a)

    # Eliminating all the un-physical results
    if len(sc) == 3:
        sc1 = sc[2]
    elif len(sc) == 2:
        sc1 = sc[1]
    elif len(sc) == 1:
        sc1 = sc[0]
    else:
        sc1 = sc
    # Plotting different fates in a different color
    if sc1 == []:
        plt.plot(O_m1,O_l1, 'ro', color='red')
        plt.axis([0,2.7,-1,3])
        plt.title('Different fates of the universe')
    elif sc1 > 1:
        plt.plot(O_m1,O_l1, 'ro', color='green')
    elif sc1 < 1 and sc1 > 0:
        plt.plot(O_m1,O_l1, 'ro', color='blue')
    elif sc1 < 0:
        plt.plot(O_m1,O_l1, 'ro', color='red')