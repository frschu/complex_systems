"""
This program intergrates the loretz model of convection numerically, plots the trajectories and plots the return map.
F. Schuessler, 18.07.2014
"""
import pylab as plt
import numpy as np
import random
from  mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
from scipy.signal import argrelextrema

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

def direct_sphere(position, radius, dim = 3):
    norm = 0
    x = np.zeros(dim)
    for k in range(dim):
        x_k = (random.gauss(0,1))
        norm += x_k ** 2
        x[k] = x_k
    upsilon = radius * random.uniform(0, 1) ** (1 / dim)
    for k in range(dim):
        x[k] = x[k] * upsilon / np.sqrt(norm) + position[k]
    return x

def progress(i, n, step):
    if i == 0:
        print('progress')
    if i % (n * step) < (n * step + 1) % (n * step):
        print('%.2f' %(i / n))

# physical contral parameters
s = 10
b = 8 / 3
r = 28
# r > r_c = 24.74
a = np.sqrt(b * (r - 1))
R_p = [ a,  a, r - 1]
R_m = [-a, -a, r - 1]

# conditions
n_runs = 1
t_max = 100
t_steps = 10001
dt = t_max / t_steps
t = np.linspace(0, t_max, t_steps)   # time grid
position = R_p
position = [10,10,30]
#position = [random.uniform(-20, 20), random.uniform(-30, 30), random.uniform(0, 50)]
delta_max = 10 ** (-6)

# PLOTTING
plt.close('all')

# SUBPLOT 1: Trajectories and separation in x
fig1 = plt.figure(figsize=plt.figaspect(0.75))
fig1.suptitle('Lorentz model')
# trajectories
ax = fig1.add_subplot(111, projection = '3d')
qs = []
for i in range(n_runs):
    if n_runs == 1:
        q0 = position                                   # use 'position' as point
    else:
        q0 = direct_sphere(position, delta_max)        #sampling a point with d(q0, position) < delta_max
    q = odeint(f, q0, t, Dfun=Dfun).transpose()
    qs.append(q)
    ax.plot(qs[i][0], qs[i][1], qs[i][2])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Trajectories, $\Delta (t = 0) = %.1e$'%delta_max)

# calculating the separation between trajectories in x
if n_runs >= 2:
    x = qs[0][0]
    x1 = qs[1][0]
    delta_x = abs(x - x1)
    ax = fig1.add_subplot(122)
    ax.semilogy(t, delta_x)
    ax.set_ylabel('$\Delta (t)$')
    ax.set_xlabel('t')
    ax.set_title('separation between trajectories (in x)')

# SUBPLOT 2: Time evolution of x and z
fig2 = plt.figure(figsize=plt.figaspect(1.3))
fig2.suptitle('Lorentz model: x(t) and z(t)')
ax = fig2.add_subplot(211)
x = qs[0][0]
ax.plot(t, x)
ax.set_title('x(t)')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax = fig2.add_subplot(212)
z = qs[0][2]
ax.plot(t, z)
ax.set_title('z(t)')
ax.set_xlabel('t')
ax.set_ylabel('z')

# SUBPLOT 3: Z maxima and return map
fig3 = plt.figure(figsize=plt.figaspect(0.5))
fig3.suptitle('Lorentz model: return map')
j = argrelextrema(z, np.greater)
ax = fig3.add_subplot(121)
ax.plot(t, z)
ax.plot(t[j], z[j], 'o')
ax.set_title('z(t)')
ax.set_xlabel('t')
ax.set_ylabel('z')

# return map
ax = fig3.add_subplot(122)
ax.plot(z[j][0:-2], z[j][1:-1], '.')
ax.set_title('z(t)')
ax.set_xlabel('t')
ax.set_ylabel('z')
min = np.ceil(min(z[j])) 
max = np.floor(max(z[j])) 
ax.plot([min, max], [min, max], 'k-')
ax.set_xlim((min, max))
ax.set_ylim((min, max))

fig3.show()
fig1.show()
