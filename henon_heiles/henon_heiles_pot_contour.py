"""
This program plots a contour plot of the potential of the Henon-Heiles system.
F. Schuessler, 20.07.2014
"""

import math
import pylab as plt 
import numpy as np

def U(Q, q):
    return(0.5 * (Q ** 2 + q ** 2) + (Q ** 2) * q - q ** 3 / 3)

# Contour plot of potential
Q_lim = [-1.5, 1.5]
q_lim = [-1.5, 2.]
q_steps = 1001
Qs = np.linspace(Q_lim[0], Q_lim[1], q_steps)
qs = np.linspace(q_lim[0], q_lim[1], q_steps)
Q, q = np.meshgrid(Qs, qs)
outer = np.arange(5., 1., -0.5)
inner = [1 /(math.factorial(i)) for i in range(1, 10)]
levels = np.concatenate((outer, inner))


fig1 = plt.figure(figsize=plt.figaspect(0.5))
fig1.suptitle("Henon-Heiles system: Contour plot of potential")
ax = fig1.add_subplot(111)
cont = ax.contour(Q, q, U(Q, q), levels)
ax.clabel(cont, inline=1, fontsize=10)
ax.set_xlabel('Q')
ax.set_ylabel('q')
fig1_name = 'henon_heiles_pot_contour.png'
fig1.savefig(fig1_name)
fig1.show()
