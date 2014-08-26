"""
This program calculates the lyaponov exponent of the Henon-Heiles system, approximated by that of the Poincare map.
As of now, this is only done for E = 1/6, where most of the phase space is in the chaotic regime. 
F. Schuessler, 30.07.2014
"""

import math
import pylab as plt 
import numpy as np
import random
import os

def U(Q, q):
    return(0.5 * (Q ** 2 + q ** 2) + (Q ** 2) * q - q ** 3 / 3)

def P(Q, q, p, E):
    return(np.sqrt(2 * E - p ** 2 - 2 * U(Q, q)))

plt.close('all')

#initial conditions
E = 1/6
n_runs = 2
dx = 10 ** (-6)
q0 = -0.03
p0 = 0
t_max = 10000
output_dir = 'E_%.3f' % E
tn_max = [40] # max value of linear behaviour, adjust by eyesight (-> graph)
tn_max0 = 100 # max value for plot

# reading reference file
file_ref = 'henon_heiles_zeros_E_0.167_q0_-0.030000_p0_0.000000_chaos.npy'     # reference file
if file_ref in os.listdir(output_dir):
    print('reference file:')
    print(file_ref)
    q0s, P0s, p0s = np.transpose(np.load(output_dir + '/' + file_ref))[[1, 2, 3]]
    len0 = len(q0s)
    t_bar = t_max / len0
else: print("Reference file does not exist!")

# plotting
fig = plt.figure()
fig.suptitle('Henon-Heiles: lyaponov exponent\n$\\tau = %.1f$'%t_bar)
axq = fig.add_subplot(121)
axp = fig.add_subplot(122)

# reading other files, calculating lyapunov exp
delta_q = []
delta_p = []
i = 0
for file in os.listdir(output_dir):
    if file.startswith('henon_heiles_zeros'):
        if file.endswith('chaos.npy'):
            if file != file_ref:
                file_str = str.rsplit(file, '_')
                q0i = float(file_str[-4])
                p0i = float(file_str[-2])
                dist = np.sqrt((q0 - q0i)**2 + (p0 - p0i)**2) 
                if dist < (1.1 * dx):
                    tn = np.arange(tn_max[i])
                    A = np.vstack([tn, np.ones(tn_max[i])]).T
                    print('file no. %i:'%(i+1))
                    print(file)
                    qs, Ps, ps = np.transpose(np.load(output_dir + '/' + file))[[1, 2, 3]]
                    len_min = min(len(qs), len0)
                    dq0 = abs(q0 - q0i)
                    dp0 = abs(p0 - p0i)
                    if dp0 == 0:
                        dp0 = abs(p0s[1] - ps[1])
                    delta_q.append(np.log(abs(q0s[1:len_min] - qs[1:len_min])/ dq0))
                    delta_p.append(np.log(abs(p0s[1:len_min] - ps[1:len_min])/ dp0))
                    d_q = delta_q[i][:tn_max[i]]
                    l_q, offset_q = np.linalg.lstsq(A, d_q)[0]
                    print('l_q = ', l_q)
                    axq.plot(range(tn_max0), delta_q[i][:tn_max0])
                    axq.plot([0, tn_max[i]],[offset_q, l_q * tn_max[i] + offset_q], label='$\lambda = %.3f$' %l_q)
                    d_p = delta_p[i][:tn_max[i]]
                    l_p, offset_p = np.linalg.lstsq(A, d_p)[0]
                    print('l_p = ', l_p)
                    axp.plot(delta_p[i][:tn_max0])
                    axp.plot([0, tn_max[i]],[offset_p, l_p * tn_max[i] + offset_p], label='$\lambda = %.3f$' %l_p)
                    i += 1

axq.set_xlabel('t_n')
axq.legend(loc=4)
axq.set_ylabel('log(|dq(t)|/dq0)')
axp.set_xlabel('t_n')
axp.legend(loc=4)
axp.set_ylabel('log(|dp(t)|/dp0)')
fig.show()
t_bar = t_max / len_min
print('mean return time of Poincare map t_bar = ', t_bar)
