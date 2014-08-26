"""
The program calculates the attractor of the discrete henon and lozi maps:

    x1(t + 1) = 1 - a * abs(x1(t)) ** m + x2(t)
    x2(t + 1) = b * x1(t)
    m = 1, 2 (lozi and henon map, respectively)

F. Schuessler, 20.08.2014
"""
import pylab as plt
import numpy as np
import math
import random
import os
import time

def direct_sphere(radius, dim = 2):
    norm = 0
    x = np.zeros(dim)
    for k in range(dim):
        x_k = (random.gauss(0,1))
        norm += x_k ** 2
        x[k] = x_k
    upsilon = radius * random.uniform(0, 1) ** (1 / dim)
    for k in range(dim):
        x[k] = x[k] * upsilon / np.sqrt(norm)
    return x

def direct_sphere2d(radius):
    x = random.gauss(0, 1)
    y = random.gauss(0, 1)
    upsilon = radius * math.sqrt(random.uniform(0, 1) / (x ** 2 + y ** 2))
    x *= upsilon
    y *= upsilon
    return (x, y)

def progress(i, n, steps=10):
    if i == 0:
        print('progress')
        str0 = '%.' + '%i'%math.log(steps, 10) + 'f'
        print(str0 % 0.)
    else:
        j = np.floor(i / n * steps)
        if (i - 1) / n * steps < j:
            str1 = '%.' + '%i'%math.log(steps, 10) + 'f'
            print(str1 %(i / n))

##########################################
#  MAIN
##########################################

# initial conditions
print('sampling')
position = [0, 0]
radius = 10 ** -2
n = 10 ** 5
x = np.zeros(n)
y = np.zeros(n)
for i in range(n):
    progress(i, n)
    x[i], y[i] = direct_sphere2d(radius)

# iteration
m = 2
a = 1.4
b = 0.3
if m == 1:
    title = 'Lozi map'
if m == 2:
    title = 'Henon map'
t_max = 100              # number of iterations
t_plot = np.arange(t_max - 3, t_max)  # times at which configuration is plotted
plt.close('all')
fig1 = plt.figure(figsize=plt.figaspect(0.75))
fig1.suptitle(title)
print('iteration')
for t in range(t_max):
    progress(t, t_max)
    if t == 0:
        ax = fig1.add_subplot(2, 2, 1)
    if t in t_plot:
        ax = fig1.add_subplot(2, 2, np.where(t_plot == t)[0] + 2)
        ax.plot(x, y, ',')
        ax.set_xlabel('x1')
        ax.set_ylabel('x2')
        ax.set_title('t = %i' %t)
    xold = x
    x = 1 - a * abs(x) ** m + y
    y = b * xold


# box counting
dn_l = 1                  # steps of powers
n_l_max = 5             # minimum length scale = 2 ** (-n_l_max)
n_l = int(n_l_max / dn_l)
ls = 1 * 2 ** np.linspace(-n_l_max, -dn_l, n_l) # length scales
boxes = [{} for l in range(n_l)]          # dictionary for the boxes of lowest length scale

start = time.time()
print('box counting')
for i in range(n):                          # box counting for smalles length scale
    xi = x[i]
    yi = y[i]
    for j  in range(n_l):
        l = ls[j]
        x_l = xi - xi % l
        y_l = yi - yi % l
        xy_l = tuple([x_l, y_l])
        if xy_l in boxes[j]:
            boxes[j][xy_l] += 1
        else:
            boxes[j][xy_l] = 1
elapsed = (time.time() - start)
print('elapsed time = ', elapsed)

# Calculate number of boxes according to number of dots inside
start = time.time()
N_l = np.ones(n_l)
n_boxes = dict()
print('analyse boxes')
for j in range(n_l):
    l = ls[j]
    N_l[j] = len(boxes[j])
    n_boxes[l] = dict()
    for n_box in boxes[j].values():
        if n_box in n_boxes[l]:
            n_boxes[l][n_box] += 1
        else:
            n_boxes[l][n_box] = 1
del boxes
elapsed = (time.time() - start)
print('elapsed time = ', elapsed)
 
pow10 = round(np.log(n) / np.log(10))    # 10 ** pow10 points
print('total number of points n = 10 ** %i ' % pow10)
filename = 'henon_map_n_boxes_10**%i.npy'%(pow10)
np.save(filename, np.array([n_boxes], dtype=object))                      # save n_boxes to file!

# Calculate fractal dimension
l0 = 1.54             # largest length scale
x = - np.log(ls / l0) / np.log(2)
y = np.log(N_l) / np.log(2)
A = np.vstack([x, np.ones(n_l)]).T
reg_range = np.where((N_l > 100) * ((N_l / n) < 0.5))      # range for regression, 100 < N_l < 0.5 * n_points
D, offset_D = np.linalg.lstsq(A[reg_range], y[reg_range])[0]
reg_label = 'y = %.2f * x + %.2f' %(D, offset_D)
print('fractal dimension D = ', D)


# PLOTTING fractal dimension regression
fig2 = plt.figure(figsize=plt.figaspect(0.75))
fig2.suptitle('Lorentz model')
ax = fig2.add_subplot(111)
ax.plot(x, y)
ax.plot(x[reg_range], y[reg_range], 'x')
ax.plot(x, D * x + offset_D, label=reg_label)
ax.set_xlabel('$-log_2(l / l_0)$')
ax.set_ylabel('$log(N_l)$')
ax.legend(loc=4)
ax.set_title('Fractal dimension')

fig1.show()
fig2.show()
