"""
This program integrates the Henon-Heiles system numerically.
In this case, two different methods (scipy odeint and 3rd and 4th order symplect Ruth algorithm) are confronted. The latter two conserve energy.
F. Schuessler, 20.07.2014
"""

import math
import pylab as plt 
import numpy as np
from scipy.integrate import odeint

def U(Q, q):
    return(0.5 * (Q ** 2 + q ** 2) + (Q ** 2) * q - q ** 3 / 3)

def H(x):
    Q = x[0]
    q = x[1]
    P = x[2]
    p = x[3]
    return(0.5 * (P ** 2 + p ** 2) + U(Q, q))

def P(Q, q, p, E):
    return(np.sqrt(2 * E - p ** 2 - 2 * U(Q, q)))

def f(x, t):
    Q = x[0]
    q = x[1]
    P = x[2]
    p = x[3]
    return([P, p, -Q * (2 * q + 1), q ** 2 - q - Q ** 2])

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

# Ruth 4rd order
def Ruth_4(x, dt):
    return(symp4(symp4(symp4(symp4(x, 3), 2), 1), 0))

def symp4(x, i):
    c1 = 1 / (2 * (2 - 2 ** (1 / 3)))
    c2 = (1 - 2 ** (1 / 3)) / (2 * (2 - 2 ** (1 / 3)))
    d1 = 2 * c1
    d2 = - 2 ** (4/3) * c1
    d4 = 0
    c = [c1, c2, c2, c1]
    d = [d1, d2, d1, d4]
    P = x[2] - d[i] * dt * x[0] * (2 * x[1] + 1)
    p = x[3] - d[i] * dt * (x[0] ** 2 - x[1] ** 2 + x[1])
    Q = x[0] + c[i] * dt * P
    q = x[1] + c[i] * dt * p
    return([Q, q, P, p])

def Dfun(x, t):
    Q = x[0]
    q = x[1]
    Df1 = [0, 0, 1, 0]
    Df2 = [0, 0, 0, 1]
    Df3 = [-(2 * q + 1), -2 * Q, 0, 0]
    Df4 = [-2 * Q, (2 * q - 1), 0, 0]
    return([Df1, Df2, Df3, Df4])

def progress(i, n, step):
    if i == 0:
        print('progress')
    if i % (n * step) < (n * step + 1) % (n * step):
        print('%.2f' %(i / n))

plt.close('all')

#initial conditions
E = 1 / 24
Q0 = 0
q0 = 0
p0 = 0
P0 = P(Q0, q0, p0, E)
x0 = [Q0, q0, P0, p0]
R = np.sqrt(2 * (E - U(Q0, q0)))

# integration
n_runs = 1
t_max = 1
t_steps = 100000
dt = t_max / t_steps
t = np.linspace(0, t_max, t_steps)   # time grid
xs = []
for i in range(n_runs):
    if n_runs > 1:
        x0 = direct_surface(R)         #sampling a point with d(q0, position) < delta_max
    x = np.array(odeint(f, x0, t, Dfun=Dfun))
    E = [H(xi) for xi in x]
    x = x.transpose()
    xs.append(x)

# symplectic integration, Ruth 3rd order
x = x0
xR3 = [x0]
for i in range(t_steps -1):
    x = (Ruth_3(x, dt))
    xR3.append(x)
ER3 = [H(xi) for xi in xR3]
xR3 = np.array(xR3).transpose()

# symplectic integration, Ruth 4th order
x = x0
xR4 = [x0]
for i in range(t_steps -1):
    x = (Ruth_4(x, dt))
    xR4.append(x)
ER4 = [H(xi) for xi in xR4]
xR4 = np.array(xR4).transpose()


fig2 = plt.figure()
fig2.suptitle("Henon-Heiles system: Q(t)")
ax = fig2.add_subplot(221)
ax.plot(t, xs[0][0])
ax.plot(t, xR3[0], 'r')
ax.plot(t, xR4[0], 'g')
ax.set_xlabel('t')
ax.set_ylabel('Q(t)')
ax = fig2.add_subplot(222)
ax.plot(t, xs[0][1], '--')
ax.plot(t, xR3[1], 'r')
ax.plot(t, xR4[1], 'g')
ax.set_xlabel('t')
ax.set_ylabel('q(t)')
ax = fig2.add_subplot(223)
ax.plot(t, E)
ax.plot(t, ER3, 'r')
ax.plot(t, ER4, 'g')
ax.set_xlabel('t')
ax.set_ylabel('E(t)')
ax = fig2.add_subplot(224)
ax.plot(t, xs[0][3])
ax.plot(t, xR3[3], 'r')
ax.plot(t, xR4[3], 'g')
ax.set_xlabel('t')
ax.set_ylabel('p(t)')

fig2.show()
