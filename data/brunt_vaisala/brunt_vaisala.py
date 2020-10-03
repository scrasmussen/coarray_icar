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

import scipy as sy
from scipy import fftpack


f = open(sys.argv[1])

header = ['image', 'timestep', 'identifier', 'env_potential_temperature']

particles = pd.read_csv(f, sep='\s+',header=None, names=header)
num_t = particles['timestep'].max()
unique_particles = particles.identifier.unique()
num_particles = len(list(particles.identifier.unique()))


f = open(sys.argv[2])
l = f.readline()
dims = [int(n) for n in l.split()]
nx = dims[0]; nz = dims[1]; ny = dims[2]
ximages = dims[3]; yimages = dims[4]
dz_value = 500
time=0

nx += 1
ny += 1
nz += 1
header2 = ['image','timestep','identifier', 'exists','moved', 'x', 'y', 'z',
          'u', 'v', 'w', 'z_meters', 'z_interface', 'pressure','temperature',
          'potential_temperature', 'velocity', 'water_vapor','cloud_water',
          'relative_humidity']
df = pd.read_csv(f, sep='\s+',header=None, names=header2)
df = df[df.timestep > 0]

particles = pd.merge(particles,df, how='inner', on=['image','timestep','identifier'])

i = -1
for p in unique_particles:
    particles.loc[particles.identifier == p, 'identifier'] = i
    i-=1

for i in range(-1,-num_particles-1,-1):
    particles.loc[particles.identifier == i, 'identifier'] = -i
    # particles.replace({'identifier' : int(change_i)}, new_i, inplace=True)


if ('no_rh' in f.name):
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
dtype = type(particles.timestep[1])

bv_all = pd.DataFrame(columns=['timestep','identifier','n'])
for id in particles.identifier.unique():
    bv = pd.DataFrame(columns=['timestep','identifier','n'])
    z = particles[particles.identifier == id].z_meters.values
    theta = particles[particles.identifier == id].env_potential_temperature.values
    # theta = particles[particles.identifier == id].potential_temperature.values
    # theta = particles[particles.identifier == id].temperature.values

    dz = np.diff(z)

    # -------- N^2_d = g* (dln(theta) / dz) ---------
    # if (relative_humidity == False):
    dln = np.diff(np.log(theta))
    dln_dz = dln / dz
    n2 = g * dln_dz

    # -------- N^2_m = g/T (dT/Dz + gamma_m)(1+Lq_s/RT) - g/(1-q_w) dq_w/dz ----
    # else:

    n = n2
    n = np.sqrt(n2)

    bv.n = n
    bv.timestep = np.arange(len(n),dtype=dtype)
    bv.identifier = id

    bv_all = bv_all.append(bv)
    # sys.exit()
    # particles = particles.merge(bv, how='outer')
    # print(bv.n)


particles = particles.merge(bv_all, how='outer')
# sys.exit()

particles.z_meters /= 1000
particles.pressure /= 1000

rows = 3 # 2
cols = 2
fig = plt.figure()
ax = fig.add_subplot(rows,cols,1)

# -----plot temperature-----
# particles.temperature -= 273.15
sc = ax.scatter(particles.timestep, particles.env_potential_temperature,
           # ax.scatter(particles.z_meters, particles.temperature,
           cmap=discrete_cmap, c=particles.identifier, marker='.')
# ax.set_xlabel("timesteps")
ax.set_ylabel("environmental potential temp. (K)")
# ax.set_xlabel("elevation (km)")
# ax.set_ylabel("temperature (K)")

# -----plot pressure-----
ax2 = fig.add_subplot(rows,cols,2)
ax2.scatter(particles.timestep, particles.n,
            # ax2.scatter(particles.z_meters, particles.pressure,
            cmap=discrete_cmap, c=particles.identifier, marker='.')
ax2.set_ylabel("Brunt-Vaisala freq.")
# ax2.set_xlabel("timesteps")
# plt.setp(ax2,yticklabels=[])
# plt.setp(ax2.get_yticklabels(), visible=False)



if relative_humidity == False:
    # -----plot temp over time-----
    ax3 = fig.add_subplot(rows,cols,3)
    ax3.scatter(particles.timestep, particles.temperature,
                    cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4 = fig.add_subplot(rows,cols,4)
    ax4.scatter(particles.timestep, particles.pressure,
                cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax4.set_xlabel("timesteps")
    ax4.set_ylabel("pressure")

if relative_humidity == True:
    rh = particles[particles.relative_humidity >= 1.0]
    no_rh = particles[particles.relative_humidity < 1.0]
    # ------- no rh -------
    # -----plot temp over time-----
    ax3 = fig.add_subplot(rows,cols,3)
    ax3.scatter(no_rh.timestep, no_rh.temperature,
                cmap=discrete_cmap, c=no_rh.identifier, marker='.')
    ax3.set_xlabel("timesteps")
    ax3.set_ylabel("temp (K)")
    # -----plot pressure over time-----
    ax4 = fig.add_subplot(rows,cols,4)
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


# np.diff(np.where(np.diff(np.sign(x)))[0][:-1])

    X = fftpack.fft(x)
    freqs = fftpack.fftfreq(len(x))
    # print(table_data)

    table_data = np.c_[table_data, ['{:,.5f}'.format(ave_freq), ave_period, '']]

# sys.exit()
table_data = table_data[:,1:]



# sys.exit()
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
    # table_colors[2][particle] =
    # table_colors[2][particle] = matplotlib.colors.to_hex(sc.to_rgba(particle+1))
    # # print(matplotlib.colors.to_rgba(matplotlib.colors.to_hex(sc.to_rgba(particle))))
    # table_colors[1][particle] = matplotlib.colors.to_hex(sc.to_rgba(particle+1))
    # table_colors[0][particle] = matplotlib.colors.to_hex(sc.to_rgba(particle+1))




if (relative_humidity == False):
    title_end = "dry air parcels"

title = "Brunt-Vaisala Frequency \n"
title += str(num_particles) + " particles, " + str(num_t) + " timesteps"
# title += ", " + f.name.split('/')[1]
plt.suptitle(title)

plt.show()
# fig.save(filename, pdf=False, pgf=True)
print("Fin!")
