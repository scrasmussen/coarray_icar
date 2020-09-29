from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import math
import copy
import sys

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

# sys.exit()

comparison = particles[particles.timestep == 0].relative_humidity.values == particles[particles.timestep == 1].relative_humidity.values
if (comparison.all()):
    relative_humidity = False
else:
    relative_humidity = True



# --- function that gets run every animation timestep ---
blue_cmap = plt.cm.Blues
discrete_cmap = plt.get_cmap('tab20b')
rh_cmap = plt.get_cmap('gist_gray')
# discrete_cmap = plt.get_cmap('tab10')
# discrete_cmap = plt.get_cmap('Set1')
# discrete_cmap = plt.get_cmap('Paired')


g = 9.81
dtype = type(particles.timestep[0])

bv_all = pd.DataFrame(columns=['timestep','identifier','n'])
for id in particles.identifier.unique():
    bv = pd.DataFrame(columns=['timestep','identifier','n'])
    z = particles[particles.identifier == id].z_meters.values
    # theta = particles[particles.identifier == id].potential_temperature.values
    theta = particles[particles.identifier == id].temperature.values

    dz = np.diff(z)

    # -------- N^2_d = g* (dln(theta) / dz) ---------
    # if (relative_humidity == False):
    dln_dz = np.diff(np.log(theta)) / dz
    n2 = g * dln_dz

    # -------- N^2_m = g/T (dT/Dz + gamma_m)(1+Lq_s/RT) - g/(1-q_w) dq_w/dz ----
    # else:


    n = np.sqrt(n2)
    bv.n = n
    bv.timestep = np.arange(len(n),dtype=dtype)
    bv.identifier = id

    bv_all = bv_all.append(bv)
    # particles = particles.merge(bv, how='outer')
    # print(bv.n)


particles = particles.merge(bv_all, how='outer')
# sys.exit()

particles.z_meters /= 1000
particles.pressure /= 1000

rows=2 # 2
cols = 2
fig = plt.figure()
ax = fig.add_subplot(rows,cols,1)

# -----plot temperature-----
# particles.temperature -= 273.15
ax.scatter(particles.timestep, particles.potential_temperature,
           # ax.scatter(particles.z_meters, particles.temperature,
           cmap=discrete_cmap, c=particles.identifier, marker='.')
# ax.set_xlabel("timesteps")
ax.set_ylabel("potential temperature (K)")
# ax.set_xlabel("elevation (km)")
# ax.set_ylabel("temperature (K)")

# -----plot pressure-----
ax2 = fig.add_subplot(2,2,2)
ax2.scatter(particles.timestep, particles.n,
            # ax2.scatter(particles.z_meters, particles.pressure,
            cmap=discrete_cmap, c=particles.identifier, marker='.')
ax2.set_xlabel("Brunt-Vaisala freq.")
plt.setp(ax2,yticklabels=[])
plt.setp(ax2.get_yticklabels(), visible=False)



if relative_humidity == False:
    # -----plot temp over time-----
    ax3 = fig.add_subplot(2,2,3)
    ax3.scatter(particles.timestep, particles.temperature,
                    cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4 = fig.add_subplot(2,2,4)
    ax4.scatter(particles.timestep, particles.pressure,
                cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax4.set_xlabel("timesteps")
    ax4.set_ylabel("pressure")

if relative_humidity == True:
    rh = particles[particles.relative_humidity >= 1.0]
    no_rh = particles[particles.relative_humidity < 1.0]
    # ------- no rh -------
    # -----plot temp over time-----
    ax3 = fig.add_subplot(2,2,3)
    ax3.scatter(no_rh.timestep, no_rh.temperature,
                cmap=discrete_cmap, c=no_rh.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4 = fig.add_subplot(2,2,4)
    ax4.scatter(no_rh.timestep, no_rh.pressure,
                cmap=discrete_cmap, c=no_rh.identifier, marker='.')
    ax4.set_xlabel("timesteps")
    ax4.set_ylabel("pressure")
    # -------- rh ---------
    # -----plot temp over time-----
    ax3.scatter(rh.timestep, rh.temperature,
                cmap=rh_cmap, c=rh.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4.scatter(rh.timestep, rh.pressure,
                cmap=rh_cmap, c=rh.identifier, marker='.')
    ax4.set_xlabel("timesteps")
    ax4.set_ylabel("pressure")
    fig.legend(['Grey Scale when relative humidity < 1.0'], loc='lower left')




title="Temperature and Pressure of \n"
title += str(num_particles) + " particles, " + str(num_t) + " timesteps"
title += ", " + f.name.split('/')[1]
plt.suptitle(title)

plt.show()
# fig.save(filename, pdf=False, pgf=True)
print("Fin!")
