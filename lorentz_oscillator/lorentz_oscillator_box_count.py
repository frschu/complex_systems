"""
This program calculates the fractal dimension and it's heigher moments for the lorentz oscillator. 
F. Schuessler, 18.07.2014
"""
import pylab as plt
import numpy as np
import random
from  mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import os
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

def progress(i, n, steps=1000):
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
# r > r_c = 24.74
a = np.sqrt(b * (r - 1))
R_p = np.array([ a,  a, r - 1])
R_m = np.array([-a, -a, r - 1])

# conditions
n_runs = 1000            # number of trajectories
t_min = 100             # time for "warm up"
t_max = t_min + 1000     # runtime for each trajectory
dt = 0.01               # time step
n_t = int(t_max / dt)
t = np.linspace(0, t_max, n_t)     # time grid
l0 = 2 ** 6             # largest length scale
dn = 1                 # steps of powers
n_max = 0 + 4
n_l = int(n_max / dn)
ls = 1 * 2 ** np.linspace(-n_max, -dn, n_l) # length scales
boxes = [{} for l in range(n_l)]          # dictionary for the boxes of lowest length scale
n_total = 0

start = time.time()
print('integration')
for i in range(n_runs):
    progress(i, n_runs)
    q0 = [random.uniform(-20, 20), random.uniform(-30, 30), random.uniform(0, 50)]
    q1 = odeint(f, q0, t, Dfun=Dfun)[int(t_min / dt):]
    for pos in q1:                          # box counting for smalles length scale
        for j  in range(n_l):
            l = ls[j]
            xyz = tuple(pos - (pos % l))
            if xyz in boxes[j]:
                boxes[j][xyz] += 1
            else:
                boxes[j][xyz] = 1
    n_total += len(q1)

# Calculate number of boxes according to number of dots inside
N_l = np.ones(n_l)
n_boxes = dict()
for j in range(n_l):
    l = ls[j]
    N_l[j] = len(boxes[j])
    n_boxes[l] = dict()
    for n in boxes[j].values():
        if n in n_boxes[l]:
            n_boxes[l][n] += 1
        else:
            n_boxes[l][n] = 1
del boxes
elapsed = (time.time() - start)
print('elapsed time = ', elapsed)

pow10 = round(np.log(n_total) / np.log(10))    # 10 ** pow10 points
print('total number of points n_total = 10 ** %i ' % pow10)
filename = 'lorentz_oscillator_n_boxes_10**%i.npy'%(pow10)
np.save(filename, np.array([n_boxes], dtype=object))                      # save n_boxes to file!

# Calculate fractal dimension
x = - np.log(ls / l0) / np.log(2)
y = np.log(N_l) / np.log(2)
A = np.vstack([x, np.ones(n_l)]).T
reg_range = np.where((N_l > 100) * ((N_l / n_total) < 0.5))      # range for regression, 100 < N_l < 0.5 * n_points
D, offset_D = np.linalg.lstsq(A[reg_range], y[reg_range])[0]
reg_label = 'y = %.2f * x + %.2f' %(D, offset_D)
print('fractal dimension D = ', D)

# PLOTTING
plt.close('all')
fig1 = plt.figure(figsize=plt.figaspect(0.75))
fig1.suptitle('Lorentz model')
ax = fig1.add_subplot(111)
ax.plot(x, y)
ax.plot(x[reg_range], y[reg_range], 'x')
ax.plot(x, D * x + offset_D, label=reg_label)
ax.set_xlabel('$-log_2(l / l_0)$')
ax.set_ylabel('$log(N_l)$')
ax.legend(loc=4)
ax.set_title('Fractal dimension')

fig1.show()


# PLOTTING number of boxes with n points
fig2 = plt.figure(figsize=plt.figaspect(0.75))
fig2.suptitle('Lorentz model')
for i in range(1, 5):
    l = ls[i - 1]
    n = np.array(list(n_boxes[l].keys()))
    n.sort()
    alpha = np.log(n / n_total) / np.log(l / l0)
    m = np.array([n_boxes[l][key] for key in n])
    f = - np.log(m) / np.log(l / l0)
    q = 0
    ax = fig2.add_subplot(2,2,i)
    ax.plot(alpha, q * alpha - f)
    ax.set_xlabel('$\\alpha$')
    ax.set_ylabel('$%.1f \cdot \\alpha - f(\\alpha)$' %q)
    ax.set_title('l = 2 ** %i' %(int(np.log(l) / np.log(2))))
fig2.show()
