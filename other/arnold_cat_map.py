"""
The Arnold cat map is a 2 dimension nonlinear symplectic map,
acting on points on the torus [0,1] x [0,1] 
    x(n + 1) =  (x(n) + y(n)) % 1
    y(n + 1) =  (x(n) + 2 * y(n)) % 1

fixed point: (x, y) = (0, 0)

eigenvalues: a_pm = (3 +- math.sqrt(5)) / 2
normed eigenvector (x_pm, y_pm) = ((math.sqrt(2 / (5 +- math.sqrt(5)))), (1 +- math.sqrt(5) / math.sqrt(2 * (5 +- math.sqrt(5)))))
Numerical values:   a_m = 0.382, (x_m, y_m) = (0.851, 0.049)
                    a_p = 2.618, (x_p, y_p) = (0.526, 1.588)

choose between the following initial conditions:
    (1) eigenvector for a_m, and two neighboring lines (contraction)
    (2) uniform random dots
    (3) square around intersection of two eigenvectors
    (4) ellipse
"""

import math
import pylab as plt
import numpy as np
import random

def f1(x, y):
    return (x + y) % 1

def f2(x, y):
    return (x + 2 * y) % 1

n = 100 ** 2
n_runs = 12
initial = 3
if initial == 0:
    # initial conditions: eigenvector for a_m / a_p
    a_m = (3 - math.sqrt(5)) / 2
    y_m = a_m - 1
    a_p = (3 + math.sqrt(5)) / 2
    y_p = a_p - 1
    xs = [1 / n * a for a in range(n)] 
    ys = [x * y_m % 1 for x in xs]
    delta_y = 0.01
    ys1 = [(x * y_m  + delta_y) % 1 for x in xs]
    ys2 = [(x * y_m  - delta_y) % 1 for x in xs]
elif initial == 1:
    # uniform random dots
    xs = []
    ys = []
    for i in range(n):
        xs.append(random.uniform(0,1))
        ys.append(random.uniform(0,1))
elif initial == 2:
    # intersection at x_i = 0.4472, y_i = 0.7236
    # initial points within the square x_i +/- delta_x, y_i +/- delta_y
    # check results for n_runs = 12, n = 100 ** 2!
    x_i = 0.4472
    y_i = 0.7236
    delta_x = 0.005
    delta_y = delta_x
    n_1 = int(math.sqrt(n))
    xs = [x_i - delta_x * (2 / n_1 * a - 1) for a in range(n_1) for b in range(n_1)]
    ys = [y_i - delta_y * (2 / n_1 * b - 1) for a in range(n_1) for b in range(n_1)]
elif initial == 3:
    # fixpoint x_i = y_i = 0
    # initial points within the square x_i +/- delta_x, y_i +/- delta_y
    x_i = 0.0
    y_i = 0.0
    delta_x = 0.005
    delta_y = delta_x
    n_1 = int(math.sqrt(n))
    xs = [(x_i - delta_x * (2 / n_1 * a - 1)) % 1 for a in range(n_1) for b in range(n_1)]
    ys = [(y_i - delta_y * (2 / n_1 * b - 1)) % 1 for a in range(n_1) for b in range(n_1)]
else:
    # initial conditions: ellipse
    n_1 = n / 2
    xs = [1 / n_1 * a for a in range(n_1)] 
    ys = [0.5*math.sqrt(x * (1 - x)) for x in xs]
    xs *= 2
    ys += [-y for y in ys]
    ys = [y + 1/2 for y in ys]

plt.figure()
plt.plot(xs,ys, ',')
if initial == 0:
    plt.plot(xs,ys1, ',')
    plt.plot(xs,ys2, ',')
for i in range(n_runs):
    xsnew = [f1(xs[j],ys[j]) for j in range(n)]
    ys = [f2(xs[j],ys[j]) for j in range(n)]
    if initial == 0:
        ys1 = [f2(xs[j],ys1[j]) for j in range(n)]
        ys2 = [f2(xs[j],ys2[j]) for j in range(n)]
    xs = xsnew
    #plt.plot(xs,ys, ',')
plt.plot(xs,ys, '.')
if initial == 0:
    plt.plot(xs,ys1, '.')
    plt.plot(xs,ys2, '.')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0,1)
plt.ylim(0,1)
plt.title('Arnold cat map')
plt.show()
