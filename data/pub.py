from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import copy
import sys
from tools.my_setup import *

# from pubplot import Document
# from pubplot.document_classes import acm_sigconf


# These were added to my_setup
# plt.rc('font', family='serif')
# plt.rcParams['font.size'] = 10
# plt.rcParams['axes.linewidth'] = 2

# --- be sure to add ---
# AND fig = plt.figure(figsize=(4,3))
# AND plt.tight_layout() before plt.show()



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

# --- function that gets run every animation timestep ---
blue_cmap = plt.cm.Blues
discrete_cmap = plt.get_cmap('tab20b')
rh_cmap = plt.get_cmap('gist_gray')


particles.z_meters /= 1000
particles.pressure /= 1000

rows = 1 #3 # 2
cols = 1 # 2
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(rows,cols,1)


filename='elevation_temp_dry'
    # particles.temperature -= 273.15
ax.scatter(particles.temperature, particles.z_meters,
           cmap=discrete_cmap, c=particles.identifier, marker='.')
ax.set_xlabel("temperature (K)")
ax.set_ylabel("elevation (km)")


plt.tight_layout()
plt.show()
# fig.save(filename, pdf=False, pgf=True)
print("Fin!")
