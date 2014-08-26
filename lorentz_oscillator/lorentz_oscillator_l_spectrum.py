"""
This program intergrates the loretz model of convection numerically and calculates the lyapunov spectrum for a given initial condition.
F. Schuessler, 18.07.2014
"""
import pylab as plt
import numpy as np
import math
from scipy.integrate import odeint
from numpy.linalg import norm
import time

def f(q, t):
    x = q[0]
    y = q[1]
    z = q[2]
    f_x = s * (-x + y)
    f_y = x * (-z + r) - y
    f_z = x * y - b * z
    return([f_x, f_y, f_z])

def Dfun(q, t):
    x = q[0]
    y = q[1]
    z = q[2]
    return([[-s, s, 0], [-z + r, -1, -x], [y, x, -b]])

def Vol(x):                     # Volume in 1, 2, and 3 dimensions, absolute values
    dx1 = (x[1] - x[0])
    dx2 = (x[2] - x[0])
    dx3 = (x[3] - x[0])
    Vol1 = norm(dx1)
    Vol2 = np.cross(dx1, dx2)
    Vol3 = np.dot(dx3, Vol2)
    return(np.array([Vol1, norm(Vol2), abs(Vol3)]))

def proj(u, v):                 # projection of v onto u
    return(np.dot(u, v) / np.dot(u, u) * u)

def Gram_Schmidt(x, V_0, dim):   # Gram Schmidt orthonormalization
    l0 = V_0[0]                 # length of output vectors
    dx = [x[i] - x[0] for i in range(dim + 1)]
    dx[2] -= proj(dx[1], dx[2])
    dx[3] -= (proj(dx[1], dx[3]) + proj(dx[2], dx[3]))
    for i in range(1, dim + 1):
        dx[i] *= l0 / norm(dx[i])
        x[i] = x[0] + dx[i]
    return(x[1:])

def progress(i, n, steps=10):
    if i == 0:
        print('progress')
        str0 = '0.0'
        print(str0)
    else:
        j = np.floor(i / n * steps)
        if (i - 1) / n * steps < j:
            str1 = '%.' + '%i'%math.log(steps, 10) + 'f'
            print(str1 %(i / n))

##########################################
#  INTEGRATION
##########################################

# physical control parameters
s = 10
b = 8 / 3
r = 28
# r > r_c = 24.74       # r_c = value above which irregular convection takes place
c = np.sqrt(b * (r - 1))
R_p = [ c,  c, r - 1]   # unstable fixed points
R_m = [-c, -c, r - 1]   

# initial conditions

x_ref = 10
y_ref = 10
z_ref = 30
dx0 = 10 ** (-2)
x0 = np.array([x_ref, x_ref + dx0, x_ref, x_ref])
y0 = np.array([y_ref, y_ref, y_ref + dx0, y_ref])
z0 = np.array([z_ref, z_ref, z_ref, z_ref + dx0])
X0 = np.dstack((x0, y0, z0))[0]
n_runs = len(x0)    # total number of trajectories

# initialising
dt = 0.01       # time step of integration
n_max = 1000    # number of renormalizations
tau = 0.2      # Q = 0 intersects before renormalization
dim = 3         # dimension of the problem
x = X0.copy()
a = np.ones([n_max, 3])
V_0 = Vol(X0)

start = time.time()
for n in range(n_max):
    progress(n, n_max)
    for i in range(dim + 1):
        x[i] = odeint(f, x[i], [0, tau], Dfun=Dfun)[1]
    a[n] = Vol(x) / V_0
    x[1:] = Gram_Schmidt(x, V_0, dim)
elapsed = (time.time() - start)
print('elapsed time = ', elapsed)

# calculation lyapunov exponents
l = np.zeros(dim)
a = np.transpose(a)
for i in range(dim):
    l[i] = sum(np.log(a[i])) / ((n_max) * tau) - sum(l[:i])
print('l1 = l_max = ', l[0])
print('l2 = ', l[1])
print('l3 = ', l[2])
