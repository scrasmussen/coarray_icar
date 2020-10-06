from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import copy
import sys
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

# oparticles = particles[['timestep','x','y','z_meters', 'identifier', 'image']]
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


particles.z_meters /= 1000
particles.pressure /= 1000

rows = 3 # 2
cols = 2
fig = plt.figure()
ax = fig.add_subplot(rows,cols,1)

if True: # -----plot temperature-----
    filename='elevation_temp_dry'
    # particles.temperature -= 273.15
    sc = ax.scatter(particles.temperature, particles.z_meters,
    # ax.scatter(particles.z_meters, particles.temperature,
               cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("elevation (km)")
    # ax.set_xlabel("elevation (km)")
    # ax.set_ylabel("temperature (K)")

if True: # -----plot pressure-----
    ax2 = fig.add_subplot(3,2,2)
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


comparison = particles[particles.timestep == 0].relative_humidity.values == particles[particles.timestep == 1].relative_humidity.values
if (comparison.all()):
    relative_humidity = False
else:
    relative_humidity = True


if relative_humidity == False:
    # -----plot temp over time-----
    ax3 = fig.add_subplot(3,2,3)
    ax3.scatter(particles.timestep, particles.temperature,
                    cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4 = fig.add_subplot(3,2,4)
    ax4.scatter(particles.timestep, particles.pressure,
                cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax4.set_xlabel("timesteps")
    ax4.set_ylabel("pressure")

if relative_humidity == True:
    rh = particles[particles.relative_humidity >= 1.0]
    no_rh = particles[particles.relative_humidity < 1.0]
    # ------- no rh -------
    # -----plot temp over time-----
    ax3 = fig.add_subplot(3,2,3)
    ax3.scatter(no_rh.timestep, no_rh.temperature,
                cmap=discrete_cmap, c=no_rh.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4 = fig.add_subplot(3,2,4)
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



# def



# ---- plot table ----
t = particles[particles.identifier == 1].timestep.to_numpy()

ax5 = fig.add_subplot(3,2,6)
t_x = 3
t_y = 1
table_dim = (t_x,t_y)
table_data = np.zeros(table_dim, dtype=float)
for row in range(1,num_particles+1):
    x = particles[particles.identifier == row].temperature.to_numpy()
    x = x[~np.isnan(x)]
    x = x - x[0]  # center data
    period = np.diff(np.where(np.diff(np.sign(x)) < 0))
    ave_period = np.average(np.diff(np.where(np.diff(np.sign(x)) < 0)))
    ave_freq = 1 / ave_period
    X = fftpack.fft(x)
    freqs = fftpack.fftfreq(len(x))
    # print(table_data)

    table_data = np.c_[table_data, ['{:,.5f}'.format(ave_freq), ave_period, '']]
table_data = table_data[:,1:]


col_names = ['ave. frequency','ave. period','particle color']
table_data = table_data.transpose()
table_df = pd.DataFrame(data=table_data, columns=col_names)
table_df.index += 1
c_table = pd.plotting.table(ax5,table_df.transpose(), rowLabels=None, loc='center')
                  # cellColours=table_colors)
ax5.axis("off")

table_colors = np.chararray((t_y + 2,num_particles))
for particle in range(0,num_particles):
    c_table.get_celld()[(3,particle)].set_color(matplotlib.colors.to_hex(sc.to_rgba(particle+1)))






plt.tight_layout()
plt.show()
# fig.save(filename, pdf=False, pgf=True)
print("Fin!")
