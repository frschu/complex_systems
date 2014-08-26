"""
This program calculates the intersections of trajectories of the henon heils system with the Q = 0 plane
F. Schuessler, 22.07.2014
"""

import pylab as plt 
import numpy as np
import os

# choose Energy:
j = 0
E = [1/6, 1/8, 1/12][j]
output_dir = 'E_%.3f' % E

for file in os.listdir(output_dir):
    if file.startswith('henon_heiles_trajec'):
        name2 = file.split('jec_')[1]
        filename_zeros = output_dir + '/henon_heiles_zeros_' + name2
        print('searching zeros for ' + name2.split('.npy')[0] )
        xs = np.load(output_dir + '/' + file)           # loading file
        Qs = np.transpose(xs)[0]
        max_d = max(Qs[:-1] - Qs[1:])
        epsilon = max_d * 1.01                          # threshold for search for zeros
        a = np.where((Qs < 0) * (Qs > -epsilon))[0]     # candidates for zeros
        idx_zeros = []                                  # indexes of zeros
        j = 0
        while j < len(a) - 1:
            while a[j + 1] == a[j] + 1:
                j += 1
                if j == len(a) - 1: break
            if a[j] < len(Qs) - 1:
                if Qs[a[j] + 1] > 0:
                    idx_zeros.append(a[j])
            j += 1
        np.save(filename_zeros, xs[idx_zeros])
        print('epsilon  = ', epsilon)
        print('number of zeros: ', len(idx_zeros))
