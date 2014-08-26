"""
Analysis of the nonlinear autonomous ODE
    dx/dt = y * z
    dy/dt = -2 * x * z
    dz/dt = x * y

fixpoints:  

Friedrich Schuessler, 29.04.2014
"""      

import math
import pylab as plt
import numpy as np
import random
from scipy.integrate import odeint
pi = math.pi

def f(x, t):
    return x[0] * (3 - x[0] - x[1]), x[1] * (x[0] - 1)

def df(x, t):
    return 3 - 2 * x[0] - x[1], -x[0], x[1], x[0] - 1

def progress(i, n, step):
    if i % (n * step) < (n * step + 1) % (n * step):
        print(i / n)

show_trajectory = False
show_phase_space = True
lines = True
title1 = 'autonomous nonlinear ODE'

# time grid
t0 = 0
t_max = 5.
dt = 0.01
if lines:
    ts = np.arange(t0, t_max + dt, dt)
else:
    ts = np.arange(t0, t_max + dt, t_max)

# setting initial conditions
initial = 1
n = 20 ** 2
if initial == 0:
    # single point, arbitrary initial condition
    n = 1
    x0 = [10]
    y0 = [3]
    title2 = '\ninitial condition: single point (%.2f, %.2f)' %(x0[0], y0[0])
if initial == 1:
    # uniform random dots
    xlim = 0,5
    ylim = 0,5
    x0 = [random.uniform(xlim[0], xlim[1]) for i in range(n)]
    y0 = [random.uniform(ylim[0], ylim[1]) for i in range(n)]
    title2 = '\ninitial condition: %i uniform random numbers' %n
elif initial == 2:
    # fixpoint x_i = y_i = 0
    # initial points within the square x_i +/- delta_x, y_i +/- delta_y
    x_i = 0.0
    y_i = 0.0
    delta_x = 0.1
    delta_y = delta_x
    n_1 = int(math.sqrt(n))
    x0 = [(x_i - delta_x * (2 / n_1 * a - 1)) % xlim for a in range(n_1) for b in range(n_1)]
    y0 = [y_i - delta_y * (2 / n_1 * b - 1) for a in range(n_1) for b in range(n_1)]
    title2 = '\ninitial condition: square around fixpoint (%.1f, %.1f), dx = %.3f' %(x_i, y_i, delta_x)
elif initial == 3:
    # fixpoint x_i = 1, y_i = 2
    # initial points within the square x_i +/- delta_x, y_i +/- delta_y
    x_i = 1
    y_i = 2
    delta_x = 0.1
    delta_y = delta_x
    n_1 = int(math.sqrt(n))
    x0 = [(x_i - delta_x * (2 / n_1 * a - 1)) % xlim for a in range(n_1) for b in range(n_1)]
    y0 = [y_i - delta_y * (2 / n_1 * b - 1) for a in range(n_1) for b in range(n_1)]
    title2 = '\ninitial condition: square around fixpoint (%.1f, %.1f), dx = %.3f' %(x_i, y_i, delta_x)

xs = []
for i in range(n):
    progress(i, n, 0.1)
    xs.append(odeint(f, [x0[i], y0[i]], ts, Dfun=df))

if show_trajectory:
    plt.figure()
    for i in range(n):
        plt.plot(ts, xs[i][:, 0])
    plt.title(title1 + ', x(t)' + title2)
    plt.xlabel('t')
    plt.ylabel('x(t)')
    plt.show()

if show_phase_space:
    plt.figure()
    plt.plot(x0, y0, 'b.')
    for i in range(n):
        if lines :
            plt.plot(xs[i][:, 0], xs[i][:,1],  '-')
        else:
            plt.plot(xs[i][1:, 0], xs[i][1:,1],  'g.')
    plt.title(title1 + ', y(x)' + title2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
