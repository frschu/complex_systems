"""
Analysis of the driven, damped nonlinear pendulum:
    dx/dt = y
    dy/dt = -omega_0 * x + F * cos(omega_F * t)

    omega_0 = 8
    omega_F = 2
    F = 10

Friedrich Schuessler, 29.04.2014
"""      

import math
import pylab as plt
import numpy as np
pi = math.pi

def runge_kutta_4(x, v, t, dt):
    # (k, l) <--> (dx/dt, dv/dt)
    k1 = f1(x, v, t)
    l1 = f2(x, v, t)
    k2 = v + 0.5 * l1 * dt      # not the general form!
    l2 = f2(x + 0.5 * k1 * dt, v + 0.5 * l1 * dt, t + 0.5 * dt)
    k3 = v + 0.5 * l2 * dt      # "-"
    l3 = f2(x + 0.5 * k2 * dt, v + 0.5 * l2 * dt, t + 0.5 * dt)
    k4 = v + l3 * dt            # "-"
    l4 = f2(x + k3 * dt, v + l3 * dt, t + dt)
    x  = x + dt * (1/6 * (k1 + k4) + 1/3 * (k2 + k3))
    v  = v + dt * (1/6 * (l1 + l4) + 1/3 * (l2 + l3))
    return x, v

def progress(i, n_steps):
    if i==0: print('progress')
    if i % (n_steps / 100) == 0: 
        print('%.2f'%(i/n_steps))
    return 0

def f1(x1, x2, t):
    return x2

def f2(x1, x2, t):
    omega_0 = 8
    omega_F = 2
    F = 10
    return -omega_0 * x1 + F * math.cos(omega_F * t)

initial = 0
if initial == 0:
    # arbitrary initial conditions
    x0 = 0.0
    v0 = 0.0
    t0 = 0.0
    t_max = 100
else:
    print('bad input')

dt = 2 ** (-7)
n_steps = int((t_max - t0) / dt)

x = x0
v = v0
t = t0
xs = [x0]
vs = [v0]
ts = [t0]
for i in range(n_steps):
    progress(i, n_steps)
    x, v = runge_kutta_4(x, v, t, dt)
    xs.append(x)
    vs.append(v)
    ts.append(t)
    t += dt

show_trajectory = False
show_phase_space = True

if show_trajectory:
    plt.figure()
    plt.plot(ts, xs)
    plt.title('forced harmonic oscillator: trajectory phi(t)')
    plt.xlabel('t')
    plt.ylabel('$\phi(t)$')
    plt.show()

if show_phase_space:
    plt.figure()
    plt.plot(xs, vs, '-')
    plt.title('forced harmonic oscillator: dx/dt over x')
    plt.xlabel('x')
    plt.ylabel('v = dx/dt')
#    plt.xlim(0,2*pi)
#    plt.ylim(0,2*pi)
    plt.show()

