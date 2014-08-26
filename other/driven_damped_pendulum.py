"""
Analysis of the driven, damped nonlinear pendulum:
    d^2x / dt^2 + gamma_p * dx /dt + (alpha - beta * cos(t)) * sin(x) = 0
with driving force 
    F_drive = h0 * cos(t)
and gamma_p = gamma / omega
    alpha   = g / (L omega^2)
    beta    = h0 / L

The diff eq. is split up into due 1st order diff eq:
    v = dx/dt
    dv/dt = f(x, v, t)

For different beta, there is different behavior:
    0 < beta < 0.55, all orbits collaps into the same periodic orbit with T = 2pi
    0.55 <= beta < beta_d: period doubling: T = k * 2pi
    where the critical beta_d = 0.64018
    beta_d <= beta: no regularities
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
    gamma_p = 0.03
    alpha = 0.5
    beta = 0.63
    return -(gamma_p * x2 + (alpha - beta * math.cos(t)) * math.sin(x1))

initial = 0
if initial == 0:
    # arbitrary initial conditions
    x0 = 0.0
    v0 = 0.01
    t0 = 0.0
    t_max = 15000
    t_equi = 13100         # approx. time at which equilibrium is reached
elif initial == 1:
    # initial conditions for equilibrium for beta = 0.63, dt = 2 ** (-7)
    x0 = 96.92113251780207860975
    v0 = -1.43109660831644247558
    t0 = 13100.000
    t_max = 16000
    t_equi = 13700          # approx. time at which equilibrium is reached
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
xp = []
vp = []
xp1 = []
vp1 = []
ts1 = []
strobo = 2 *  pi        # Time intervall, at which xp and vp are taken ("strobo phase space")
#strobo = 2**(-5)*pi     # Time intervall, at which xp and vp are taken ("strobo phase space")
counter = 0
for i in range(n_steps):
    progress(i, n_steps)
    x, v = runge_kutta_4(x, v, t, dt)
    xs.append(x)
    vs.append(v)
    ts.append(t)
    t += dt
    if (t-dt) == t_equi:
        print('equilibrium initial values:')
        print('x0 = %.20f\nv0 = %.20f\nt0 = %.3f'%(x, v, t - dt))
    if t > t_equi:
        if (t - dt) % strobo < dt:
            xp.append(x % (2 * pi))
            vp.append(v % (2 * pi))
            if counter % 4 == 0:
                xp1.append(x % (2 * pi))
                vp1.append(v % (2 * pi))
                ts1.append(t)
            counter += 1
# np.var
plt.figure()
plt.plot(ts, xs)
plt.title('nonlinear driven pendulum: trajectory phi(t)')
plt.xlabel('t')
plt.ylabel('$\phi(t)$')
plt.show()

plt.figure()
plt.plot(xp, vp, '.')
plt.title('nonlinear driven pendulum: dx/dt over x')
plt.xlabel('x')
plt.ylabel('v = dx/dt')
plt.xlim(0,2*pi)
plt.ylim(0,2*pi)
plt.show()

"""
plt.figure()
plt.plot(ts1, xp1)
plt.title('phase space: dx/dt over x')
plt.xlabel('t')
plt.ylabel('x')
plt.ylim(0,2*pi)
plt.show()
"""
