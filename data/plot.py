from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys


if (len(sys.argv) < 2):
    sys.exit("Error: too few arguments for `plot.py [graph_data.txt]`")


frame_delay_ms=10
frame_delay_ms=25
frame_delay_ms=50
frame_delay_ms=100
frame_delay_ms=200
# frame_delay_ms=300


#
turn_off_graphs=True
# turn_off_graphs=False



# ---- read input data ----
f = open(sys.argv[1])
l = f.readline()
dims = [int(n) for n in l.split()]
nx = dims[0]; nz = dims[1]; ny = dims[2]
ximages = dims[3]; yimages = dims[4]

header = ['image','timestep','identifier', 'exists','moved', 'x', 'y', 'z',
          'u', 'v', 'w', 'k', 'pressure','temperature', 'potential_temperature',
          'velocity', 'water_vapor','mixing_ratio']

if (f.name == 'four-200-time-step.txt'):
    header = ['image','timestep','identifier', 'exists','moved', 'x', 'y', 'z',
          'u', 'v', 'w', 'k', 'pressure','temperature', 'potential_temperature',
          'velocity','water_vapor']

particles = pd.read_csv(f, sep='\s+',header=None, names=header)
# sys.exit()
num_t = particles['timestep'].max()
num_particles = list(particles.identifier.unique())
unique_particles = particles.identifier.unique()

print("ARTLESS: fix normalization")
norm = matplotlib.colors.Normalize(vmin=particles['temperature'].min()-10,
                                   vmax=particles['temperature'].max())
# ---- Set up colors ----
cmap = plt.cm.Greens
cmap_c = cmap(norm(particles.temperature.values))

fig = plt.figure(figsize=plt.figaspect(0.5))
if not turn_off_graphs:
    ax = fig.add_subplot(1,2,1,projection='3d')
else:
    ax = fig.add_subplot(1,1,1,projection='3d')
ax.set_xlim(0,nx); ax.set_ylim(0,ny); ax.set_zlim(0,nz)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.zaxis.set_major_locator(MaxNLocator(integer=True))

ax.set_xlabel="x axis"
ax.set_ylabel="y axis"
ax.set_zlabel="z axis"

# --- plot hill lines ---
def hill_contour(hillX,hillY):
    hillZ = np.full_like(hillX,0.0)
    ids = 1
    jds = 1
    ide = nx
    jde = ny
    dz_value = 500
    surface_z = 0
    n_hills = 1.0
    hill_height = 1000.0

    hillZ = (np.sin((hillX-ids)/((ide-ids)/n_hills)*2*3.14159-3.14159/2)+1)/2 *\
            (np.sin((hillY-jds)/((jde-jds)/n_hills)*2*3.14159-3.14159/2)+1)/2
    hillZ = ((hillZ * hill_height ) + (surface_z + dz_value / 2.0)) / dz_value
    return hillZ

fineness=1
hillx = np.linspace(1,nx,nx*fineness)
hilly = np.linspace(1,ny,ny*fineness)
hillX,hillY = np.meshgrid(hillx,hilly)
hillZ = hill_contour(hillX,hillY)


# --- plot image lines ---
def plot_image_lines():
    ax.set_title("Particle Movement", y=1.05)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.zaxis.set_major_locator(MaxNLocator(integer=True))
    for i in range (1,ximages):
        ax.plot(xs=[(i*nx)/ximages,(i*nx)/ximages], ys=[0,ny], color='black')
    for i in range (1,yimages):
        ax.plot(xs=[0,nx], ys=[(i*ny)/yimages,(i*ny)/yimages], color='black')
    ax.contour3D(hillX,hillY,hillZ, cmap='binary')


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
    if turn_off_graphs:
        return

    second_graph.set_title(second_graph_title)
    second_graph.set_xlim(0,num_t)
    second_graph.set_ylim(second_graph_min,second_graph_max)
    third_graph.set_title(third_graph_title, y=-0.3)
    third_graph.set_xlim(0,num_t)
    third_graph.set_ylim(third_graph_min,third_graph_max)



particles = particles[['timestep','x','y','z']]
p = particles[ (particles.timestep == 0) ]
old = ax.scatter3D(p.x, p.y, p.z, marker='o', edgecolor='black',
                       color=cmap_c[0])


# --- function that gets run every animation timestep ---
def updateFig(*args):
    global t, old

    if (t == num_t):
        ax.cla()
        ax.set_xlim(0,nx); ax.set_ylim(0,ny); ax.set_zlim(0,nz)
        plot_image_lines()
        t = 0
    else:
        old.remove()

    p = particles[ (particles.timestep == t) ]
    old = ax.scatter3D(p.x, p.y, p.z, marker='o', color='green')
                       # color=cmap_c[t]) edgecolor='black')

    t += 1
    return old

    # # old style
    # if not turn_off_graphs:
    #     second_graph.cla()
    #     third_graph.cla()
    #     # set_graph_lim()
    # clear_3D_graph = False
    # if (clear_3D_graph):
    #     ax.cla()
    #     ax.set_xlim(0,nx); ax.set_ylim(0,ny); ax.set_zlim(0,nz)
    #     plot_image_lines()
    # for id in unique_particles:
    #     break
    #     if not turn_off_graphs:
    #         p_less_than = particles[ (particles.identifier == id) &
    #                                  (particles.timestep < t) ]
    #         p_time = p_less_than.timestep

    #         p_second_graph = p_less_than[second_graph_entity]
    #         second_graph.plot(p_time, p_second_graph, color='black')
    #         p_third_graph  = p_less_than[third_graph_entity]
    #         third_graph.plot(p_time, p_third_graph, color='black')

    #     # ax.view_init(10+(t/75), 45+int(t/25))

    #     p = particles[ (particles.identifier == id) &
    #                    (particles.timestep == t) ]
    #     if not p.empty:
    #         old = ax.scatter3D(p.x, p.y, p.z, marker='o', edgecolor='black',
    #                            color=cmap_c[t])
    #     else:
    #         if not turn_off_graphs:
    #             second_graph.plot(p_time[-1:], p_second_graph[-1:],
    #                               color='black', marker='x')
    #             third_graph.plot(p_time[-1:], p_third_graph[-1:], color='black',
    #                              marker='x')

    # t += 1
    # return old


# --- animation ----
# - setup graphs
plot_image_lines()
set_graph_lim()

gif    = True
gif    = False
repeat = False if gif else True

# num_t = 4
# ax.view_init(0, 0)

t = 0
ani = animation.FuncAnimation(fig, updateFig, interval=frame_delay_ms,
                              frames=num_t-1,repeat=repeat)

if (gif):
    print("Gif!")
    ani.save('./test.gif', writer='imagemagick', fps=None) # fps was 5
else:
    plt.show()

print("Fin!")
