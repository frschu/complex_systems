"""
Just for temporary use: plots calculated vectors from "henon_heiles_integrator"
F. Schuessler, 29.07.2014
"""

import math
import pylab as plt 
import numpy as np
import os

def U(Q, q):
    return(0.5 * (Q ** 2 + q ** 2) + (Q ** 2) * q - q ** 3 / 3)

def H(x):
    Q = x[0]
    q = x[1]
    P = x[2]
    p = x[3]
    return(0.5 * (P ** 2 + p ** 2) + U(Q, q))

# choose Energy:
j = 0
E = [1/6, 1/8, 1/12][j]
output_dir = 'E_%.3f' % E
dt = 0.1

fig1 = plt.figure()
fig1.suptitle("Henon-Heiles system: Q(t)")
ax1 = fig1.add_subplot(111)
#ax2 = fig1.add_subplot(222)
#ax3 = fig1.add_subplot(223)
#ax4 = fig1.add_subplot(224)
i = 0
for file in os.listdir(output_dir):
    if file.startswith('henon_heiles_trajec'):
        if i == 1:
            xs = np.load(output_dir + '/' + file)
            t_max = int(len(xs) * dt)
            t = np.arange(0, t_max, dt)
            xs = np.transpose(xs)
            E = H(xs)
            dE = E - E[0]
            ax1.plot(t, xs[0], '-o')
            ax1.set_xlabel('t')
            ax1.set_ylabel('Q(t)')
            """
            ax2.plot(t, xs[1])
            ax2.set_xlabel('t')
            ax2.set_ylabel('q(t)')
            ax3.plot(t, dE)
            ax3.set_xlabel('t')
            ax3.set_ylabel('dE(t)')
            ax4.plot(t, xs[3])
            ax4.set_xlabel('t')
            ax4.set_ylabel('p(t)')
            """
        i += 1
        if i > 1: break
fig1.show()
