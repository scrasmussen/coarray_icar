from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys

from tools.my_setup import *
# from pubplot import Document
# from pubplot.document_classes import acm_sigconf
# doc = Document(acm_sigconf)
# plt.style.use('./publication.style')


turn_off_graphs=True
# turn_off_graphs=False

# num_to_keep = int(len(unique_particles) * 0.05)
# unique_particles = np.random.choice(unique_particles, num_to_keep,
#                                     replace=False)

start_particles = particles.loc[particles['identifier'].isin(unique_particles)]
norm = matplotlib.colors.Normalize(vmin=particles['temperature'].min()-10,
                                   vmax=particles['temperature'].max())
# ---- Set up colors ----
# cmap = plt.cm.Greens
cmap = plt.cm.Blues
cmap_c = cmap(norm(particles.temperature.values))

fig = plt.figure(figsize=plt.figaspect(0.5))
# fig, ax = doc.subfigures()

if not turn_off_graphs:
    ax = fig.add_subplot(1,2,1,projection='3d')
else:
    # ax = fig.add_subplot(1,1,1,projection='3d')
    ax = fig.add_subplot(1,1,1)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_locator(MaxNLocator(integer=True))


for i in range(0,num_particles):
    if any(start_particles.identifier == i):
        change_i = start_particles.iloc[i].identifier
        start_particles.loc[start_particles.identifier == change_i, 'identifier'] = -i

q = start_particles[start_particles.timestep == 0]
q = q.sort_values(by=['temperature'])

for i in range(0,num_particles):
    change_i = q.iloc[i].identifier
    start_particles.replace({'identifier' : change_i}, i, inplace=True)





oparticles = start_particles[['timestep','x','y','z_meters', 'identifier', 'image']]
# particles = particles[['timestep','identifier','temperature','z_meters','pressure']]


# --- function that gets run every animation timestep ---
blue_cmap = plt.cm.Blues
discrete_cmap = plt.get_cmap('tab20b')



gif    = True
gif    = False
repeat = False if gif else True



plot_temperature = False # plot_pressure = True
plot_temperature = True
particles.z_meters /= 1000

plen = len(particles[particles.timestep == 1])
set1 = set(particles[particles.timestep == 1].identifier)
for i in range(1,max(particles.timestep)):
    # if (plen > len(particles[particles.timestep == i])):
    #     print("DIFF! at "+str(i)+" of "+str(len(particles[particles.timestep == i])))
    #     set2 = set(particles[particles.timestep == i].identifier)
    #     set_diff = set1 - set2
    #     print(set_diff)
    # plen = len(particles[particles.timestep == i])

    ax.scatter(i, len(particles[particles.timestep == i]), color='black')

# sys.exit()

ax.set_xlabel("num timesteps")
ax.set_ylabel("num particles")
title=f.name.split('/')[1]
plt.suptitle(title)

plt.tight_layout()
plt.show()
# fig.save(filename, pdf=False, pgf=True)
print("Fin!")
