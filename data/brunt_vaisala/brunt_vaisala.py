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


present = False
# present = True




plt.rc('font', family='serif')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2


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

bv_all = pd.DataFrame(columns=['timestep','identifier','n', 'n2'])
for id in particles.identifier.unique():
    bv = pd.DataFrame(columns=['timestep','identifier','n', 'n2'])
    z = particles[particles.identifier == id].z_meters.values
    theta = particles[particles.identifier == id].env_potential_temperature.values
    # theta = particles[particles.identifier == id].potential_temperature.values
    # theta = particles[particles.identifier == id].temperature.values

    # print(theta)
    # print("-------")

    dz = np.diff(z)
    # -------- N^2_d = (g / theta)   (dtheta / dz) ---------
    if (1 == g):
        dth = np.diff(theta)  # new
        dth_dz = dth / dz
        n2 = (g / theta[1:]) * dth_dz # new

    # -------- N^2_d = g* (dln(theta) / dz) ---------
    else:
        dln = np.diff(np.log(theta))
        dln_dz = dln / dz
        n2 = (g ) * dln_dz


    # -------- N^2_m = g/T (dT/Dz + gamma_m)(1+Lq_s/RT) - g/(1-q_w) dq_w/dz ----
    # else:

    # n = n2
    n = np.sqrt(n2)


    bv.n = n
    bv.n2 = n2
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


if (present == True):
    rows = 2
    cols = 2

fig = plt.figure()
ax = fig.add_subplot(rows,cols,1)

# -----plot temperature-----
# particles.temperature -= 273.15
sc = ax.scatter(particles.timestep, particles.env_potential_temperature,
           # ax.scatter(particles.z_meters, particles.temperature,
           cmap=discrete_cmap, c=particles.identifier, marker='.')
# ax.set_xlabel("timesteps")
ax.set_ylabel("env. potential temp. (K)")
# ax.set_xlabel("elevation (km)")
# ax.set_ylabel("temperature (K)")

# -----plot Brunt-Vaisala Freq-----
if (0 == 1):
    ax2 = fig.add_subplot(rows,cols,2)
    ax2.scatter(particles.timestep, particles.n,
                cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax2.set_ylabel("Brunt-Vaisala freq.")
    ax2.set_ylim(0)
else:
    ax2 = fig.add_subplot(rows,cols,2)
    ax2.scatter(particles.timestep, particles.potential_temperature,
                cmap=discrete_cmap, c=particles.identifier, marker='.')
    ax2.set_ylabel("potential temp. (K)")



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
    # fig.legend(['Grey Scale when relative humidity < 1.0'], loc='lower left')


# set limits of x axis
x_upper_limit = particles.timestep.max()
ax.set_xbound(0, x_upper_limit)
ax2.set_xbound(0, x_upper_limit)
ax3.set_xbound(0, x_upper_limit)
ax4.set_xbound(0, x_upper_limit)

ax.set_xbound(0, x_upper_limit)
ax2.set_xbound(0, x_upper_limit)
ax.set_xticklabels([])
ax2.set_xticklabels([])



# ax.vlines(plt.xticks()[0], ax.get_ybound()[0],ax.get_ybound()[1],


# ax.axvline(plt.xticks()[0][1],
#           linestyle='dashed', markevery=100)


where = plt.xticks()[0][1:-1]

for an_axe in fig.get_axes():
    for w in where:
        an_axe.axvline(w, linestyle='--', c='black', zorder=0)



# import matplotlib.gridspec as gridspec
# gs = gridspec.GridSpec(rows,cols)
# gs.update(hspace=0.0)



if (present == True):
    plt.tight_layout()
    plt.subplots_adjust(hspace = 0.05)

    plt.show()
    print("Fin!")
    sys.exit()




# bv_all[bv_all.identifier==1].n.median()

# ---- plot table ----
t = particles[particles.identifier == 1].timestep.to_numpy()

ax5 = fig.add_subplot(3,1,3)
t_x = 5
t_y = 1
table_dim = (t_x,t_y)
table_data = np.zeros(table_dim, dtype=float)
for row in range(1,num_particles+1):
    x = particles[particles.identifier == row].env_potential_temperature.to_numpy()
    x = x[~np.isnan(x)]
    x = x - x[0]  # center data
    period = np.diff(np.where(np.diff(np.sign(x)) < 0))
    ave_period = np.average(np.diff(np.where(np.diff(np.sign(x)) < 0)))
    ave_freq = 1 / ave_period
    # bv radians per second
    bv_freq_rad = bv_all[bv_all.identifier==row].n2.median()
    print("bv = ", bv_freq_rad)
    bv_freq = bv_freq_rad * 0.159155
    # ave_freq = ave_freq**2
    percent_diff = (abs(ave_freq - bv_freq)/ bv_freq) * 100
# np.diff(np.where(np.diff(np.sign(x)))[0][:-1])

    X = fftpack.fft(x)
    freqs = fftpack.fftfreq(len(x))
    # print(table_data)
    print(freqs.min(), freqs.max())
    idx = np.argmax(np.abs(x))
    freq = freqs[idx]
    # freq_in_hertz = abs(freq * fra)

    # table_data = np.c_[table_data, ['{:,.5f}'.format(ave_freq), ave_period, '']]
    table_data = np.c_[table_data, ['{:,.2f}'.format(ave_period),
                                    '{:,.5f}'.format(bv_freq),
                                    '{:,.5f}'.format(ave_freq),
                                    '{:,.5f}'.format(percent_diff),
                                    '']]
                                    # '{:,.5f}'.format(abs(freqs).mean()),



# sys.exit()
table_data = table_data[:,1:]



# sys.exit()
col_names = ['period','b.v. freq. (Hz)','freq. (Hz)','% diff.',
'particle color']

table_data = table_data.transpose()
table_df = pd.DataFrame(data=table_data, columns=col_names)
table_df.index += 1
c_table = pd.plotting.table(ax5,table_df.transpose(), rowLabels=None, loc='center')
                  # cellColours=table_colors)
ax5.axis("off")



table_colors = np.chararray((t_y + 2,num_particles))
for particle in range(0,num_particles):
    c_table.get_celld()[(t_x,particle)].set_color(matplotlib.colors.to_hex(sc.to_rgba(particle+1)))
    # table_colors[2][particle] =
    # table_colors[2][particle] = matplotlib.colors.to_hex(sc.to_rgba(particle+1))
    # # print(matplotlib.colors.to_rgba(matplotlib.colors.to_hex(sc.to_rgba(particle))))
    # table_colors[1][particle] = matplotlib.colors.to_hex(sc.to_rgba(particle+1))
    # table_colors[0][particle] = matplotlib.colors.to_hex(sc.to_rgba(particle+1))




if (relative_humidity == False):
    title_end = "dry air parcels"

title = "Brunt-Vaisala Frequency \n"
title += "Dry air parcels \n"
title += str(num_particles) + " particles, " + str(num_t) + " timesteps"
# title += ", " + f.name.split('/')[1]
plt.suptitle(title)

plt.tight_layout()
plt.show()
# fig.save(filename, pdf=False, pgf=True)
print("Fin!")
