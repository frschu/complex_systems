"""
This program integrates the Henon-Heiles system numerically and calculates the lyapunov spectrum. 
The program uses an algorithm suggested by Benettin (1978, 1980), calculating the evolution of k-th volume elemtents. 
After tau intersections of the Q = 0 plane, the vectors in q and p direction are renormalized. 
F. Schuessler, 20.07.2014
"""

import math
import pylab as plt 
import numpy as np
from numpy.linalg import norm
import time

def U(q, Q):        # Potential
    return(0.5 * (q ** 2 + Q ** 2) + (Q ** 2) * q - q ** 3 / 3)

def P(q, p, Q, E):  # Impulse on surface of energy = E
    return(np.sqrt(2 * E - p ** 2 - 2 * U(q, Q)))

def H(x):           # Hamiltonian
    q = x[0]
    p = x[1]
    Q = x[2]
    P = x[3]
    return(0.5 * (p ** 2 + P ** 2) + U(q, Q))

def Ruth_3(x, dt):          # Ruth 3rd order symplectic integration
    c = [2/3, -2/3, 1]
    d = [7/24, 3/4, -1/24]
    q = x[0]
    Q = x[2]
    p = x[1] - d[0] * dt * (Q ** 2 - q ** 2 + q)
    P = x[3] - d[0] * dt * Q * (2 * q + 1)
    q += c[0] * dt * p
    Q += c[0] * dt * P
    p += - d[1] * dt * (Q ** 2 - q ** 2 + q)
    P += - d[1] * dt * Q * (2 * q + 1)
    q += c[1] * dt * p
    Q += c[1] * dt * P
    p += - d[2] * dt * (Q ** 2 - q ** 2 + q)
    P += - d[2] * dt * Q * (2 * q + 1)
    q += c[2] * dt * p
    Q += c[2] * dt * P
    return(np.array([q, p, Q, P]))

def Vol(x):             # Calculation of volume in 1 and 2 dimensions
    d1 = x[1][:2] - x[0][:2]
    d2 = x[2][:2] - x[0][:2]
    return(np.array([norm(d1), np.cross(d1, d2)]))

def progress(i, n, steps=10):
    if i == 0:
        print('progress')
        str0 = '0.' + math.log(steps, 10) * '0'
        print(str0)
    else:
        j = np.floor(i / n * steps)
        if (i - 1) / n * steps < j:
            str1 = '%.' + '%i'%math.log(steps, 10) + 'f'
            print(str1 %(i / n))

############################################################
# INTEGRATION
###########################################################

# initial conditions
j = 2                       # choose energy
E = [1/6, 1/8, 1/12][j]
if j == 0:
    q_ref = -0.03
    p_ref = 0.0
if j == 1:
    q_ref = -0.10
    p_ref = 0.
 j == 2:
    q_ref = -0.1175
    p_ref = 0.0
dq0 = 10 ** (-10)
q0 = np.array([q_ref, q_ref + dq0, q_ref])
p0 = np.array([p_ref, p_ref, p_ref + dq0])
n_runs = len(q0)            # total number of runs
Q0 = np.zeros(n_runs)       # zero plane
P0 = P(q0, p0, Q0, E)
x0 = np.dstack((q0, p0, Q0, P0))[0] # array([[q0, p0, Q0, P0], [...], [...]])

# initialising
dt = 0.01                   # time step of integration
n_max = 10                  # number of renormalizations
tau = 15                    # Q = 0 intersects before renormalization
dim = 2                     # dimension of the problem
x = x0
a = np.ones([n_max, dim])   # volume contraction factors in each dimension
d = np.zeros([dim + 1, dim])# renormalization factors
Vol_0 = Vol(x0)             # reference volume

# symplectic integration, Ruth 3rd order algorithm
qzs = [[],[],[]]            # coordinates for zero plane intersects
start = time.time()
for n in range(n_max):
    print('n =', n)
    for i in range(3):
        j = 0
        xnew = x[i]
        while j < tau:      # integration up to t = t0 + tau
            xold = xnew
            xnew = Ruth_3(xold, dt)
            if xold[2] < 0: 
                if xnew[2] > 0:                 # condition for Q = 0, P > 0
                    j += 1
                    qzs[i].append(xnew[0])
        frac = -xold[2] / (xnew[2] - xold[2])   # interpolation factor
        x[i][:2] = xold[:2] + (xnew[:2] - xold[:2]) * frac  # interpolation
        x[i][2] = 0
        x[i][3] = P(x[i][0], x[i][1], x[i][2], E)
    Vol_1 = Vol(x)
    a[n] = abs(Vol_1 / Vol_0)                   # saving volumes for evalution of l
    d[1] = (x[1][:2] - x[0][:2]) / a[n][0]      # renormalization
    d[2][0] = -Vol_0[1] * d[1][1] / Vol_0[0] ** 2
    d[2][1] = -d[1][0] * d[2][0] / d[1][1]
    for i in [1, 2]:
        x[i][:2] = x[0][:2] + d[i]              # adapting the pertubed trajectories
        x[i][3] = P(x[i][0], x[i][1], 0, E)
    for i in range(3):
        qzs[i].append(x[i][0])
elapsed = (time.time() - start)
print('elapsed time = ', elapsed)

# calculation lyapunov exponents
a = np.transpose(a)
l = [0] * dim
for i in range(dim):
    l[i] = sum(np.log(a[i])) / (n_max * tau) - sum(l[:i])
print('l1 = l_max = ', l[0])
print('l2 = ', l[1])

# Plotting trajactories
plt.close('all')
plt.figure()
plt.plot(np.log(abs(np.array(qzs[0]) - np.array(qzs[1]))/ dq0))
for i in range(n_max):
    plt.plot([i * (tau + 1), i * (tau + 1)], [0, 15], 'k')
plt.show()
