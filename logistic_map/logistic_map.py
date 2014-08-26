"""
Analysis of the logistic map
    x(n + 1) = r * x(n) * (1 - x(n))
for x in [0,1], 
    r in [0,4].

r_1 = 3
r_2 = 3.448
r_inf = 3.569945
Friedrich Schuessler, 30.04.2014
"""

import math
import pylab as plt
import numpy as np

def f(x):
    return r * x * (1 - x)

def f2(x):
    return f(f(x))

def f4(x):
    return f2(f2(x))

def f8(x):
    return f4(f4(x))

def lambda1(x):
    return(abs(r * (1 - 2 * x)))

def lambda2(x):
    return lambda1(f(x))*lambda1(x)

r_1 = 3.
r_2 = 3.448
r_inf = 3.569945
r = 3.6
n_runs = 80
x_grid = np.arange(0., 1., 0.001)
fs = [f(x) for x in x_grid]


x0 = 0.1
x = x0
xn = [x]
xn_plot = [x0]
for i in range(10):
    x = f(x)
    xn.append(x)
    xn_plot += [x,x]


xlim = (0,1)
ylim = xlim
figname = 'logistic_map.png'
title = 'logistic map: x(n+1) over x(n)'
plt.figure()
plt.plot(x_grid, x_grid, 'b--')
plt.plot(x_grid, f(x_grid))
if r_1 < r and r < r_2:
    plt.plot(x_grid, f2(x_grid))
elif r < r_inf:
    plt.plot(x_grid, f4(x_grid))
plt.plot([x0] + xn_plot[:-1], [0] +  xn_plot[1:], 'k')
plt.title(title)
plt.xlabel('x(n)')
plt.ylabel('x(n + 1)')
plt.savefig(figname)
plt.show()

