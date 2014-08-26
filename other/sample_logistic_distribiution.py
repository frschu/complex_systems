import pylab as plt
import math
import random

def f(y):
    return(mu - beta * math.log(y / (1 - y)))

mu = 0
beta = 0.5
n = 10 ** 6
x_hist = []
for i in range(n):
    y = random.uniform(0,1)
    x_hist.append(f(y))

# PLOT HISTOGRAM
# if 'normed=True', the histogram is normed over the given range – NOT over the entire space
xmin = -10
xmax = 10
n_bins = 100
width = (xmax - xmin) / n_bins
plt.figure()
n_hist, bins, patches = plt.hist(x_hist, bins=n_bins, range=[xmin, xmax], facecolor='gray', align='mid', normed=True)
plt.show()
integral = sum(n_hist) * width
print(integral)
