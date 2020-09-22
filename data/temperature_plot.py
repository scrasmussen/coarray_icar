from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys

from tools.my_setup import *

frame_delay_ms=10
# frame_delay_ms=25
# frame_delay_ms=50


#
turn_off_graphs=True
# turn_off_graphs=False

hillcolor='darkolivegreen'
hillcolor='black'


particles = particles.loc[particles['identifier'].isin(unique_particles)]
unique_particles = particles.identifier.unique()
num_particles = len(list(particles.identifier.unique()))

# print(len(num_particles))
# sys.exit()

norm = matplotlib.colors.Normalize(vmin=particles['temperature'].min()-10,
                                   vmax=particles['temperature'].max())
# ---- Set up colors ----
# cmap = plt.cm.Greens
cmap = plt.cm.Blues
cmap_c = cmap(norm(particles.temperature.values))

fig = plt.figure(figsize=plt.figaspect(0.5))
if not turn_off_graphs:
    ax = fig.add_subplot(1,2,1,projection='3d')
else:
    ax = fig.add_subplot(1,1,1,projection='3d')

def set_limits():
    ax.set_xlim(0,particles.timestep.max());  # time
    ax.set_ylim(1,num_particles);
    ax.set_zlim(1,particles.temperature.max())
    # ax.set_xlim(1,nx); ax.set_ylim(1,ny); ax.set_zlim(1,nz)
    # ax.set_xlim(1,nx); ax.set_ylim(1,ny); ax.set_zlim(1,particles['z_meters'].max())

set_limits()

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.zaxis.set_major_locator(MaxNLocator(integer=True))

ax.set_xlabel("timesteps")
ax.set_ylabel("particles")
ax.set_zlabel("temperature (K)")
# ax.set_xticklabels([])
# ax.set_yticklabels([])


for i in range(0,num_particles):
    if any(particles.identifier == i):
        change_i = particles.iloc[i].identifier
        particles.loc[particles.identifier == change_i, 'identifier'] = -i

q = particles[particles.timestep == 0]
q = q.sort_values(by=['temperature'])

for i in range(0,num_particles):
    change_i = q.iloc[i].identifier
    particles.replace({'identifier' : change_i}, i, inplace=True)


if not turn_off_graphs:
    # --- set up second graph ---
    second_graph_title  = 'temperature'
    second_graph_entity = second_graph_title.replace(' ', '_')
    second_graph        = fig.add_subplot(2,2,2)
    second_graph_max    = particles[second_graph_entity].max()
    if (particles[second_graph_entity].min() < 0):
        second_graph_min = particles[second_graph_entity].min()
    else:
        second_graph_min = 0
    second_graph_max   *= 1.1
    second_graph_min   *= 1.1


    # --- set up third graph ---
    # choose type
    third_graph_title = "density"
    third_graph_title = "water vapor"
    third_graph_title = "potential temperature"
    # third_graph_title = "mixing ratio"

    third_graph_entity = third_graph_title.replace(' ', '_')
    third_graph = fig.add_subplot(2,2,4)
    third_graph_max = particles[third_graph_entity].max()
    if (particles[third_graph_entity].min() < 0):
        third_graph_min = particles[third_graph_entity].min()
    else:
        third_graph_min = 0
    third_graph_max *= 1.1
    third_graph_min *= 1.1


# --- set graph limits ---
def set_graph_lim():
    return
    if turn_off_graphs:
        return

    second_graph.set_title(second_graph_title)
    second_graph.set_xlim(0,num_t)
    second_graph.set_ylim(second_graph_min,second_graph_max)
    third_graph.set_title(third_graph_title, y=-0.3)
    third_graph.set_xlim(0,num_t)
    third_graph.set_ylim(third_graph_min,third_graph_max)


# print(particles)
# print(particles['z_meters'].max())
# sys.exit()
def set_old(p):
    old = ax.scatter3D(p.timestep, p.identifier, p.temperature, marker='o', edgecolor='black',
                       color=cmap_c[0])
    return old


oparticles = particles[['timestep','x','y','z_meters', 'identifier', 'image']]
particles = particles[['timestep','identifier','temperature']]


# --- function that gets run every animation timestep ---
blue_cmap = plt.cm.Blues
discrete_cmap = plt.get_cmap('tab20b')

def updateFig(*args):
    global t, old, time

    ax.set_title("Particle Movement t="+str(t), y=1.05)
    if (t == num_t):
        ax.cla()
        set_limits()
        # plot_image_lines(time)
        time += 1
        t = 0
    else:
        # old.remove()
        pass

    p = particles[ (particles.timestep == t) ]

    # p.water_vapor = np.ma.masked_where(p.water_vapor == 0.0, p.water_vapor)
    # p.water

    p_color = "black"
    # print(type(p.water_vapor))
    # if (p.water_vapor > 0.0):
    #     p_color = "blue"
    old = set_old(p)
    #                    c=p.cloud_water, cmap=blue_cmap, vmin=0.00000)
                       # c=p.water_vapor, cmap=blue_cmap, vmin=0.0)
                       # color=cmap_c[t]) edgecolor='black')

    t += 1
    return old


# --- animation ----
# - setup graphs
# plot_image_lines(time)
set_graph_lim()

gif    = True
gif    = False
repeat = False if gif else True

# num_t = 4
# ax.view_init(90, 0)

t = 0

# debug
# p = particles[ (particles.timestep == 0) ]
# old = set_old(p)
# ani = animation.FuncAnimation(fig, updateFig, interval=frame_delay_ms,
#                               frames=num_t-1,repeat=repeat)

# ax.scatter()

ax.scatter(particles.timestep, particles.identifier, particles.temperature,
           cmap=discrete_cmap, c=particles.identifier)
             # marker='.', edgecolor='black')
#                    color=cmap_c[0])



i = 0

# print(len(ax.get_yticklabels()))
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# sys.exit()

if (gif):
    print("Gif!")
    # ani.save('test.gif', writer=animation.PillowWriter, fps=None,dpi=20) # fps was 5
    # ani.save('test.gif', writer=animation.ImageMagickWriter, fps=None) # fps was 5
else:
    plt.show()
print("Fin!")
