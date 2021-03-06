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
import os
from scipy import fftpack


from tools.my_setup import *

# plt.rc('text', usetex=True)

turn_off_graphs=True
turn_off_graphs=False

# original = copy.deepcopy(particles)





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
particles.pressure /= 100

rows = 2 # 2
cols = 3
fig = plt.figure()

# ax = fig.add_subplot(rows,cols,1)
# if True: # -----plot temperature-----
#     # particles.temperature -= 273.15
#     sc = ax.scatter(particles.temperature, particles.z_meters,
#     # ax.scatter(particles.z_meters, particles.temperature,
#                cmap=discrete_cmap, c=particles.identifier, marker='.')
#     ax.set_xlabel("temperature (K)")
#     ax.set_ylabel("elevation (km)")
#     # ax.set_xlabel("elevation (km)")
#     # ax.set_ylabel("temperature (K)")

# if True: # -----plot pressure-----
#     ax2 = fig.add_subplot(rows,cols,2)
#     ax2.scatter(particles.pressure, particles.z_meters,
#     # ax2.scatter(particles.z_meters, particles.pressure,
#                 cmap=discrete_cmap, c=particles.identifier, marker='.')
#     # ax2.set_xlabel("elevation (km)")
#     # ax2.set_ylabel("pressure")

#     ax2.set_xlabel("pressure (kPa)")
#     # ax2.set_ylabel("elevation (km)")
#     plt.setp(ax2,yticklabels=[])
#     plt.setp(ax2.get_yticklabels(), visible=False)


comparison = particles[particles.timestep == 0].relative_humidity.values == particles[particles.timestep == 1].relative_humidity.values
if (comparison.all()):
    relative_humidity = False
else:
    relative_humidity = True



add=1 - 3
if relative_humidity == False:
    # -----plot temp over time-----
    ax3 = fig.add_subplot(rows,cols,3+add)
    ax3.scatter(particles.timestep, particles.temperature,
                    cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4 = fig.add_subplot(rows,cols,4+add)
    ax4.scatter(particles.timestep, particles.pressure,
                cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax4.set_xlabel("timesteps")
    ax4.set_ylabel("pressure (hPa)")

if relative_humidity == True:
    rh = particles[particles.relative_humidity >= 1.0]
    no_rh = particles[particles.relative_humidity < 1.0]
    # ------- no rh -------
    # -----plot temp over time-----
    ax3 = fig.add_subplot(rows,cols,3+add)
    ax3.scatter(no_rh.timestep, no_rh.temperature,
                cmap=discrete_cmap, c=no_rh.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4 = fig.add_subplot(rows,cols,4+add)
    ax4.scatter(no_rh.timestep, no_rh.pressure,
                cmap=discrete_cmap, c=no_rh.identifier, marker='.')
    ax4.set_xlabel("timesteps")
    ax4.set_ylabel("pressure (hPa)")
    # -------- rh ---------
    # -----plot temp over time-----
    ax3.scatter(rh.timestep, rh.temperature,
                cmap=rh_cmap, c=rh.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("Temp (K)")
    # -----plot pressure over time-----
    ax4.scatter(rh.timestep, rh.pressure,
                cmap=rh_cmap, c=rh.identifier, marker='.')
    ax4.set_xlabel("timesteps")
    ax4.set_ylabel("pressure (hPa)")
    # fig.legend(['Grey Scale when relative humidity > 1.0'], loc='lower left')

title="Temperature and Pressure of \n"
title += str(num_particles) + " saturated particles, " + str(num_t) + " timesteps"
# title += ", " + f.name.split('/')[1]
# plt.suptitle(title)




# --- plot potential temperature ---
ax55 = fig.add_subplot(rows,cols,5+add)
ax55.scatter(particles.timestep, particles.potential_temperature,
            cmap=discrete_cmap, c=particles.identifier, marker='.')
ax55.set_xlabel("timesteps")
ax55.set_ylabel("potential temp. (K)")



# plot altitude and the change of potential temp, temp, pressure
ax7 = fig.add_subplot(rows,cols,6+add)
ax7.scatter(particles.temperature, particles.z_meters,
            cmap=discrete_cmap, c=particles.identifier, marker='.')
ax7.set_xlabel("Temp (K)")
ax7.set_ylabel("elevation (km)")

ax8 = fig.add_subplot(rows,cols,7+add)
# ax8.scatter(particles.temperature, particles.pressure,
#             cmap=discrete_cmap, c=particles.identifier, marker='.')
# ax8.set_xlabel("Temp (K)")
# ax8.set_ylabel("pressure (hPa)")
# ax8.invert_yaxis()
ax8.scatter(particles.pressure, particles.z_meters,
            cmap=discrete_cmap, c=particles.identifier, marker='.')
ax8.set_xlabel("pressure (hPa)")
ax8.set_ylabel("elevation (km)")

ax9 = fig.add_subplot(rows,cols,8+add)
ax9.scatter(particles.potential_temperature, particles.z_meters,
            cmap=discrete_cmap, c=particles.identifier, marker='.')
ax9.set_xlabel("potential temp. (K)")
ax9.set_ylabel("elevation (km)")


plt.tight_layout()

if 'dry' in os.path.splitext(f.name)[0]:
    filename='validation_dry_parcels'
else:
    filename='validation_saturated_parcels'

filename+=".png"
# fig.set_size_inches((4,3))
fig.set_size_inches((6,4))
plt.tight_layout()
plt.savefig(filename, pad_inches=0.0, dpi=400)
# fig.save(filename, pdf=False, pgf=True)
print("file "+filename+" outputted")
print(fig.get_size_inches(), " and ",fig.dpi)
print(f.name, "  ", os.path.splitext(f.name)[0])

# should be ~ 3.07
# is 6.4, 4.8
#
plt.show()

print("Fin!")

sys.exit()



# ---- plot table ----
t = particles[particles.identifier == 1].timestep.to_numpy()

# ax5 = fig.add_subplot(3,1,3)
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

    # table_data = np.c_[table_data, ['{:,.5f}'.format(ave_freq), ave_period, '']]
    table_data = np.c_[table_data, ['{:,.5f}'.format(ave_freq),
                                    '{:,.2f}'.format(ave_period),
                                    '']]


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




# fig.save(filename, pdf=False, pgf=True)

# plt.tight_layout()
# plt.show()
print("Fin!")
