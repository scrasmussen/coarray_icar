from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
from matplotlib import rcParams
import numpy as np
import pandas as pd
import sys

gif = False
if (len(sys.argv) < 2):
    sys.exit("Error: too few arguments for `plot.py [graph_data.txt]`")
if (len(sys.argv) == 3):
    if (sys.argv[2] == 'gif'):
        gif = True
        print("Creating video/gif")



frame_delay_ms=10 # 25
frame_delay_ms=10
turn_off_graphs=True



# ------------------------------------------------------------------------------
# read input data
# ------------------------------------------------------------------------------
f = open(sys.argv[1])
l = f.readline()
dims = [int(n) for n in l.split()]
nx = dims[0]; nz = dims[1]; ny = dims[2]
ximages = dims[3]; yimages = dims[4]
dz_value = 500
time=0
nx += 1; ny += 1; nz += 1

header = ['image','timestep','identifier', 'exists','moved', 'x', 'y', 'z',
          'u', 'v', 'w', 'z_meters', 'z_interface', 'pressure','temperature',
          'potential_temperature', 'velocity', 'water_vapor','cloud_water']
particles = pd.read_csv(f, sep='\s+',header=None, names=header)
num_t = particles['timestep'].max()
unique_particles = particles.identifier.unique()
# num_to_keep = int(len(unique_particles) * 0.05)
# unique_particles = np.random.choice(unique_particles, num_to_keep,
#                                     replace=False)
# particles = particles.loc[particles['identifier'].isin(unique_particles)]
unique_particles = particles.identifier.unique()
num_particles = list(particles.identifier.unique())


# ---- Set up colors ----
norm = matplotlib.colors.Normalize(vmin=particles['temperature'].min()-10,
                                   vmax=particles['temperature'].max())
cmap = plt.cm.Blues # cmap = plt.cm.Greens
cmap_c = cmap(norm(particles.temperature.values))
blue_cmap = plt.cm.Blues
blue_cmap.set_under(color='black')


# ---- setup graph and plot if needed ----
fig = plt.figure(figsize=plt.figaspect(0.5))
if not turn_off_graphs:
    ax = fig.add_subplot(1,2,1,projection='3d')
else:
    ax = fig.add_subplot(1,1,1,projection='3d')
ax.set_xlim(1,nx); ax.set_ylim(1,ny); ax.set_zlim(1,particles['z_meters'].max())
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.zaxis.set_major_locator(MaxNLocator(integer=True))

ax.set_zlabel("meters")
ax.set_xticklabels([])
ax.set_yticklabels([])


# --- plot hill contour lines ---
fineness=1
hillx = np.linspace(1,nx,nx*fineness)
hilly = np.linspace(1,ny,ny*fineness)
hillX,hillY = np.meshgrid(hillx,hilly)
hillZ = np.full_like(hillX,0.0)
ids = 1; jds = 1; ide = nx; jde = ny
surface_z = 0
n_hills = 1.0
hill_height = 1000.0
hillZ = (np.sin((hillX-ids)/((ide-ids)/n_hills)*2*3.14159-3.14159/2)+1)/2 *\
        (np.sin((hillY-jds)/((jde-jds)/n_hills)*2*3.14159-3.14159/2)+1)/2
hillZ = ((hillZ * hill_height ) + (surface_z + dz_value / 2.0))
hillcolor='black' # hillcolor='darkolivegreen'


# --- plot image lines ---
def plot_image_lines(time):
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title("Air Parcels", y=1.05)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.zaxis.set_major_locator(MaxNLocator(integer=True))
    z_labels = ax.get_zticklabels()
    ax.set_zlabel("meters")
    ax.set_zticklabels(['0k','1.5k','3k','4.5k','6k','7.5k','9k','10.5k'])
    ax.grid(b=None)
    for i in range (1,ximages):
        ax.plot(xs=[(i*nx)/ximages,(i*nx)/ximages], ys=[0,ny], color='black')
    for i in range (1,yimages):
        ax.plot(xs=[0,nx], ys=[(i*ny)/yimages,(i*ny)/yimages], color='black')
    ax.contour3D(hillX,hillY,hillZ, colors=hillcolor)
    ax.set_xticks([])
    ax.set_yticks([])



# ---- setup additional graphs ----
if not turn_off_graphs:
    # --- setup second graph ---
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


    # --- setup third graph ---
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


# ------------------------------------------------------------------------------
# graph 3d scatter plot, called every time step
# ------------------------------------------------------------------------------
def graph_scatter(p):
    scatter = ax.scatter3D(p.x, p.y, p.z_meters, marker='o', edgecolor='black',
                       color='blue')#cmap_c[0])
    #                  c=p.cloud_water, cmap=blue_cmap, vmin=0.00000)
    #                  c=p.water_vapor, cmap=blue_cmap, vmin=0.0)
    return scatter



# --- plot initial timestep ---
particles = particles[['timestep','x','y','z_meters','water_vapor','cloud_water'
                       ,'temperature', 'identifier']]
p = particles[ (particles.timestep == 0) ]
scatter = graph_scatter(p)


# ------------------------------------------------------------------------------
#  function that gets run every animation timestep
# ------------------------------------------------------------------------------
t_interval = 1
def updateFig(*args):
    global t, scatter, time

    ax.set_title("Particle Movement t="+str(t), y=1.05)
    if (t == num_t):
        ax.cla()
        ax.set_xlim(1,nx); ax.set_ylim(1,ny); ax.set_zlim(1,particles['z_meters'].max())
        plot_image_lines(time)
        t = 0
    else:
        scatter.remove()

    p = particles[ (particles.timestep == t) ]
    scatter = graph_scatter(p)
    t += t_interval
    return scatter


# ------------------------------------------------------------------------------
# setup animation
# ------------------------------------------------------------------------------
# decide whether to create gif or not
repeat = False if gif else True

# - setup graph and start animation
plot_image_lines(time)
set_graph_lim()
t = 0
ani = animation.FuncAnimation(fig, updateFig, interval=frame_delay_ms,
                              frames=int((num_t-1) / t_interval),repeat=repeat,
                              blit=False)



# how
# print(particles[(particles.identifier == 5) & (t > particles.timestep )]\
#       .drop(['x','y','identifier'],axis=1).to_string())

if (gif):
    vid_name = 'vid-' + f.name.replace('.txt','.mp4')
    # vid_name = vid_name.replace('-','_')

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=50, metadata=dict(artist='Me'), bitrate=900)
    ani.save(vid_name, writer=writer)
    # ani.save('test.mp4', writer=writer)

    # --- NOTE: not working on OSX, version too old? ---
    # ani.save('giftest.gif', writer=animation.ImageMagickWriter)
else:
    plt.show()
print("Fin!")
