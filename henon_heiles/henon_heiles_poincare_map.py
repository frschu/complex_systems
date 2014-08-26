"""
This program plots the poincare map of the Henon-Heiles system for given points. 
F. Schuessler, 20.07.2014
"""

import math
import pylab as plt 
import numpy as np
import os

plt.close('all')

# choose Energy:
j = 0
E = [1/6, 1/8, 1/12][j]
dir = 'E_%.3f' % E

fig1_name = 'henon_heiles_poincare_map_' + dir + '.png'
fig1 = plt.figure()
fig1.suptitle("Henon-Heiles system: Poincare-section, E = %.3f" % E)
ax = fig1.add_subplot(111)
colors = [['b', 'r', 'g', 'c', 'm', 'lightgrey', 'white', 'k'], \
        ['k', 'lightgrey', 'b', 'r', 'g', 'c', 'm', 'orange', 'lightgreen', 'y','purple'], \
        ['b', 'k', 'r', 'g', 'c', 'm', 'orange', 'lightgreen', 'darkred','purple', 'k', 'grey']][j]

# read data from file
i = 0
for file in os.listdir(dir):
    if file.startswith("henon_heiles_zeros"):
        name2 = file.split('zeros_')[1]
        print('reading data from ', file )
        f = open(dir + '/' + file, 'r')
        t_zeros = []
        q_zeros = []
        p_zeros = []
        for line in f:
            a = line.split()
            if len(a) == 3:
                t_zeros.append(float(a[0]))
                q_zeros.append(float(a[1]))
                p_zeros.append(float(a[2]))
            if len(a) == 2:
                q_zeros.append(float(a[0]))
                p_zeros.append(float(a[1]))
        f.close()
        ax.plot(q_zeros, p_zeros, '.', markersize=2, color=colors[i])
        i += 1

ax.set_xlabel('q')
ax.set_ylabel('p')
#ax.set_xlim(-0.6, 1.1)
#ax.set_ylim(-0.6, 0.6)
fig1.savefig(fig1_name)
fig1.show()

