"""
The standard map is a 2 dimension nonlinear symplectic map,
acting on points on the cylinder (q, p) in [0, 2 * pi] x R.

    q(n + 1) = q(n) + p(n + 1) % (2 * pi)
    p(n + 1) = p(n) + k * sin(q(n))
    k = const in R

fixed points:
    p(n) = 0, q(n) in {0, pi}
    observations:   (q, p) = (0, 0) not very stable (saddle)
                    (q, p) = (pi, 0) stable
                    mixing for uniform numbers with p in [-1, 1]
                    same behavoir up to p = 20
                    p ~ 31.5 seems to be another fix point

"""

import math
import pylab as plt
mport random
pi = math.pi
k = 0.1
xlim = 2 * pi

def fq(q, p):
    return (q + p) % xlim

def fp(q, p):
    return p + k * math.sin(q)

n = 100 ** 2
n_runs = 5
initial = 0
if initial == 0:
    # uniform random dots
    qs = []
    ps = []
    for i in range(n):
        qs.append(random.uniform(0,xlim))
        ps.append(random.uniform(31.400, 31.401))
elif initial == 1:
    # fixpoint q_i = p_i = 0
    # initial points within the square q_i +/- delta_x, p_i +/- delta_y
    q_i = 0.0
    p_i = 0.0
    delta_x = 0.1
    delta_y = delta_x
    n_1 = int(math.sqrt(n))
    qs = [(q_i - delta_x * (2 / n_1 * a - 1)) % xlim for a in range(n_1) for b in range(n_1)]
    ps = [p_i - delta_y * (2 / n_1 * b - 1) for a in range(n_1) for b in range(n_1)]
elif initial == 2:
    # fixpoint q_i = pi, p_i = 0
    # initial points within the square q_i +/- delta_x, p_i +/- delta_y
    q_i = pi
    p_i = 0.0
    delta_x = 0.1
    delta_y = delta_x
    n_1 = int(math.sqrt(n))
    qs = [(q_i - delta_x * (2 / n_1 * a - 1)) % xlim for a in range(n_1) for b in range(n_1)]
    ps = [p_i - delta_y * (2 / n_1 * b - 1) for a in range(n_1) for b in range(n_1)]
else:
    # initial conditions: ellipse
    n_1 = int(n / 2)
    qs = [1 / n_1 * a for a in range(n_1)]
    ps = [0.5*math.sqrt(x * (1 - x)) for x in qs]
    qs *= 2
    ps += [-y for y in ps]
    ps = [y + 1/2 for y in ps]

plt.figure()
plt.plot(qs,ps, ',')
for i in range(n_runs):
    ps = [fp(qs[j],ps[j]) for j in range(n)]
    qs = [fq(qs[j],ps[j]) for j in range(n)]
    plt.plot(qs,ps, ',')
#plt.plot(qs,ps, '.')
plt.xlabel('q')
plt.ylabel('p')
plt.xlim(0, 2 * pi)
plt.title('standard map, k = %.2f, initial = %i, n_runs = %i' %(k, initial, n_runs))
plt.show()
