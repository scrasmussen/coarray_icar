from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import copy
import sys
# import pylab as pyl

import scipy as sy
from scipy import fftpack

import pandas as pd
import sys

# ---- read input data ----
f = open(sys.argv[1])
l = f.readline()

header = ['n','time','job', 'machine']

df = pd.read_csv(f, sep='\s+',header=None, names=header)

fig, (ax0, ax1) = plt.subplots(2,1)

machine='delta'
job='put-latency-1-node'
x=df[(df.job == job) & (df.machine == machine)].n
y=df[(df.job == job) & (df.machine == machine)].time
ax0.plot(x, y, label='1 node '+machine)
job='put-latency-2-node'
x=df[(df.job == job) & (df.machine == machine)].n
y=df[(df.job == job) & (df.machine == machine)].time
ax0.plot(x, y, label='2 node '+machine)
machine='cray'
job='put-latency-1-node'
x=df[(df.job == job) & (df.machine == machine)].n
y=df[(df.job == job) & (df.machine == machine)].time
ax0.plot(x, y, label='1 node '+machine)
job='put-latency-2-node'
x=df[(df.job == job) & (df.machine == machine)].n
y=df[(df.job == job) & (df.machine == machine)].time
ax0.plot(x, y, label='2 node '+machine)

ax0.set_xlabel('Size')
ax0.set_ylabel('Put Latency (us)')
ax0.legend()
# ax0.set_title('Put Latency')

# ax.set_xlabel('Frequency in Hertz [Hz]')
# ax.set_ylabel('Frequency Domain (Spectrum) Magnitude')
# ax.set_xlim(-f_s / 2, f_s / 2)
# ax.set_ylim(-5, 110)
# ax1.set_xlim(-0.002, 0.002)

# ax1.plot(freqs, sy.log10(abs(FFT)), '.')

machine='delta'
job='pt2pt-latency-1-node'
x=df[(df.job == job) & (df.machine == machine)].n
y=df[(df.job == job) & (df.machine == machine)].time
ax1.plot(x, y, label='1 node '+machine)
job='pt2pt-latency-2-node'
x=df[(df.job == job) & (df.machine == machine)].n
y=df[(df.job == job) & (df.machine == machine)].time
ax1.plot(x, y, label='2 node '+machine)
machine='cray'
job='pt2pt-latency-1-node'
x=df[(df.job == job) & (df.machine == machine)].n
y=df[(df.job == job) & (df.machine == machine)].time
ax1.plot(x, y, label='1 node '+machine)
job='pt2pt-latency-2-node'
x=df[(df.job == job) & (df.machine == machine)].n
y=df[(df.job == job) & (df.machine == machine)].time
ax1.plot(x, y, label='2 node '+machine)

ax1.set_xlabel('Size')
ax1.set_ylabel('Point2Point Latency (us)')
ax1.legend()

ax0.set_title('MVAPICH MPI Benchmarks')
plt.show()
