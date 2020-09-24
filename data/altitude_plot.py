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


turn_off_graphs=True
turn_off_graphs=False

original = copy.deepcopy(particles)

fig = plt.figure()


rows=2 # 2
cols = 2
ax = fig.add_subplot(rows,cols,1)

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


oparticles = particles[['timestep','x','y','z_meters', 'identifier', 'image']]
# particles = particles[['timestep','identifier','temperature','z_meters','pressure']]
# particles = particles[particles.identifier < 600000000]
# particles = particles[particles.identifier < -14]



# --- function that gets run every animation timestep ---
blue_cmap = plt.cm.Blues
discrete_cmap = plt.get_cmap('tab20b')
rh_cmap = plt.get_cmap('gist_gray')
# discrete_cmap = plt.get_cmap('tab10')
# discrete_cmap = plt.get_cmap('Set1')
# discrete_cmap = plt.get_cmap('Paired')



# --- animation ----
# - setup graphs
# plot_image_lines(time)

gif    = True
gif    = False
repeat = False if gif else True

# num_t = 4
# ax.view_init(90, 0)

t = 0


plot_temperature = False # plot_pressure = True
plot_temperature = True
particles.z_meters /= 1000
particles.pressure /= 1000

if (plot_temperature == True):
    filename='elevation_temp_dry'
    # particles.temperature -= 273.15
    ax.scatter(particles.temperature, particles.z_meters,
    # ax.scatter(particles.z_meters, particles.temperature,
               cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("elevation (km)")
    # ax.set_xlabel("elevation (km)")
    # ax.set_ylabel("temperature (K)")

if True:
    # -----plot pressure-----
    ax2 = fig.add_subplot(2,2,2)
    filename='elevation_pressure_dry'
    ax2.scatter(particles.pressure, particles.z_meters,
    # ax2.scatter(particles.z_meters, particles.pressure,
                cmap=discrete_cmap, c=particles.identifier, marker='.')
    # ax2.set_xlabel("elevation (km)")
    # ax2.set_ylabel("pressure")

    ax2.set_xlabel("pressure (kPa)")
    # ax2.set_ylabel("elevation (km)")
    plt.setp(ax2,yticklabels=[])
    plt.setp(ax2.get_yticklabels(), visible=False)


relative_humidity = True


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
