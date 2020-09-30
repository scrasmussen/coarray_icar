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

header = ['image','timestep','identifier', 'exists','moved', 'x', 'y', 'z',
          'u', 'v', 'w', 'z_meters', 'z_interface', 'pressure','temperature',
          'potential_temperature', 'velocity', 'water_vapor','cloud_water',
          'relative_humidity']
df = pd.read_csv(f, sep='\s+',header=None, names=header)

fig, (ax0, ax1) = plt.subplots(2,1)


# fig, ax = plt.subplots()
# ax1.stem(freqs, np.abs(X))
# ax.set_xlabel('Frequency in Hertz [Hz]')
# ax.set_ylabel('Frequency Domain (Spectrum) Magnitude')
# ax.set_xlim(-f_s / 2, f_s / 2)
# ax.set_ylim(-5, 110)
# ax1.set_xlim(-0.002, 0.002)

# ax1.plot(freqs, sy.log10(abs(FFT)), '.')


plt.show()
