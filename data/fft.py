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

from tools.my_setup import *


turn_off_graphs=True
turn_off_graphs=False

original = copy.deepcopy(particles)





for i in range(0,num_particles):
    if any(particles.identifier == i):
        change_i = particles.iloc[i].identifier
        particles.loc[particles.identifier == change_i, 'identifier'] = -i
q = particles[particles.timestep == 0]
q = q.sort_values(by=['temperature'])
for i in range(0,num_particles):
    change_i = q.iloc[i].identifier
    new_i = i * (-1) - 1
    particles.replace({'identifier' : change_i}, new_i, inplace=True)
particles.identifier *= -1

df = particles

q = df[df.identifier == 1]
t = q.timestep.to_numpy()
x = q.temperature.to_numpy()
x = x - x[0]


X = fftpack.fft(x)
freqs = fftpack.fftfreq(len(x))

fig, (ax0, ax1) = plt.subplots(2,1)
ax0.plot(x)
ax0.set_ylabel('Amplitude')
ax0.set_xlabel('Time')
# ax1.plot(fftpack.fftfreq(len(t), np.abs(X)))


# fig, ax = plt.subplots()
ax1.stem(freqs, np.abs(X))
# ax.set_xlabel('Frequency in Hertz [Hz]')
# ax.set_ylabel('Frequency Domain (Spectrum) Magnitude')
# ax.set_xlim(-f_s / 2, f_s / 2)
# ax.set_ylim(-5, 110)
ax1.set_xlim(-0.002, 0.002)


# u = x
# FFT = sy.fft(u)
# freqs = fftpack.fftfreq(len(u), 1)

# ax1.plot(freqs, sy.log10(abs(FFT)), '.')


plt.show()
