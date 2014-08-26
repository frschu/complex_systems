"""
This program integrates the Henon-Heiles system numerically.
F. Schuessler, 20.07.2014
"""

import math
import pylab as plt 
import numpy as np
import random
import os

def U(Q, q):
    return(0.5 * (Q ** 2 + q ** 2) + (Q ** 2) * q - q ** 3 / 3)

def P(Q, q, p, E):
    return(np.sqrt(2 * E - p ** 2 - 2 * U(Q, q)))

def H(x):
    Q = x[0]
    q = x[1]
    P = x[2]
    p = x[3]
    return(0.5 * (P ** 2 + p ** 2) + U(Q, q))

# Ruth 3rd order
def Ruth_3(x, dt):
    return(symp3(symp3(symp3(x, 2), 1), 0))

def symp3(x, i):
    c = [1, -2/3, 2/3]
    d = [-1/24, 3/4, 7/24]
    P = x[2] - d[i] * dt * x[0] * (2 * x[1] + 1)
    p = x[3] - d[i] * dt * (x[0] ** 2 - x[1] ** 2 + x[1])
    Q = x[0] + c[i] * dt * P
    q = x[1] + c[i] * dt * p
    return([Q, q, P, p])

def progress(i, n, step = 0.1):
    if i == 0:
        print('progress')
    if i % (n * step) < (n * step + 1) % (n * step):
        print('%.2f' %(i / n))

############################################################
# INTEGRATION
###########################################################

# initial conditions
# choose Energy:
j = 0
E = [1/6, 1/8, 1/12][j]
if j == 0:
    dq0 = 10 ** (-6)
    q0_reg = [-0.25, 0.02, 0.04, 0.07, 0.1]
    q0_chaos = [0.3, -0.030000, -0.030000 + dq0, 0.030000]
    q0 = np.array(q0_reg + q0_chaos)
    p0_reg = [-0.4, 0, 0, 0, 0]
    p0_chaos = [0, 0, 0, 0. + dq0]
    p0 = np.array(p0_reg + p0_chaos)
if j == 1:
    q0_reg = [-0.32, -0.28, 0.07, 0.16, 0.58, 0, 0]
    q0_chaos = [-0.01, -0.10]
    q0 = np.array(q0_reg + q0_chaos)
    p0_reg = [0, 0, 0, 0, 0, -0.2, 0.2]
    p0_chaos = [0.03, 0.]
    p0 = np.array(p0_reg + p0_chaos)
if j == 2:
    q0_reg = [-0.25, -0.20, -0.18, -0.08, 0., 0.1, 0.08, 0.08]
    q0_chaos = [-0.1175]
    q0 = np.array(q0_reg + q0_chaos)
    p0_reg = [0, 0, 0, 0, 0, 0, 0.2, -0.2]
    p0_chaos = [0]
    p0 = np.array(p0_reg + p0_chaos)
n_runs = len(q0)                                            # total number of runs
n_reg = len(q0_reg)                                    # number of non chaotic runs
Q0 = np.zeros(n_runs)
P0 = P(Q0, q0, p0, E)
x0 = np.dstack((Q0, q0, P0, p0))[0]
output_dir = 'E_%.3f' % E
if not os.path.exists(output_dir): os.makedirs(output_dir)

# symplectic integration, Ruth 3rd order algorithm
dt = 0.01
for i in range(n_runs):
    x = x0[i]
    filename = output_dir + '/henon_heiles_trajec_' + output_dir + '_q0_%.6f_p0_%.6f'%(x[1], x[3])
    if i < n_reg:
        t_max = 1000                   
    else:
        t_max = 10000              
        filename += '_chaos'
    if os.path.isfile(filename + '.npy'):       # check for existing file
        print('filename: ' , filename)
        choice = input('File exists, do you really want to continue? (y / n):  ')
        if choice != 'y': continue
    t_steps = int(t_max / dt)
    xs = np.zeros([t_steps, 4])
    print('i = %i, integration for t_max = %i' %(i, t_max))
    for j in range(t_steps):
        progress(j, t_steps)
        xs[j] = x
        x = (Ruth_3(x, dt))
    print('Initial position: x0 = ', x0[i])
    print('Initial energy E0 = ', E)
    print('dE = E(t_max) - E0 = ', H(x) - E)
    np.save(filename, xs)
