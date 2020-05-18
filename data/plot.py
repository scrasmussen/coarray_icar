from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys


if (len(sys.argv) < 2):
    sys.exit("Error: too few arguments for `plot.py [graph_data.txt]`")

# ---- read input data ----
f = open(sys.argv[1])
l = f.readline()
dims = [int(n) for n in l.split()]
nx = dims[0]; nz = dims[1]; ny = dims[2]
ximages = dims[3]; yimages = dims[4]

header = ['image','timestep','identifier', 'exists','moved', 'x', 'y', 'z',
          'u', 'v', 'w', 'density', 'temperature', 'velocity']

particles = pd.read_csv(f, sep='\s+',header=None, names=header)
num_t = particles['timestep'].max()
num_particles = list(particles.identifier.unique())
field = np.zeros((nx,nz,ny,num_t))

print("ARTLESS: remove this once real data is produced")
new_temp = [i+x for i,x in enumerate(particles['temperature'])]
particles['temperature'] = new_temp

print("ARTLESS: fix normalization")
norm = matplotlib.colors.Normalize(vmin=particles['temperature'].min()-10,
                                   vmax=particles['temperature'].max())
# ---- Set up colors ----
cmap = plt.cm.Greens
cmap_c = cmap(norm(particles.temperature.values))

fig = plt.figure(figsize=plt.figaspect(0.5))

ax = fig.add_subplot(1,2,1,projection='3d')
ax.set_xlim(0,nx); ax.set_ylim(0,ny); ax.set_zlim(0,nz)
ax.set_xlabel="x axis"
ax.set_ylabel="y axis"
ax.set_zlabel="z axis"

# --- plot image lines ---
def plot_image_lines():
    for i in range (1,ximages):
        ax.plot(xs=[(i*nx)/ximages,(i*nx)/ximages], ys=[0,ny], color='black')
    for i in range (1,yimages):
        ax.plot(xs=[0,nx], ys=[(i*ny)/yimages,(i*ny)/yimages], color='black')

# ---- 3d plots ----
plot_image_lines()
temp = fig.add_subplot(2,2,2)
density = fig.add_subplot(2,2,4)
temp_max = particles.temperature.max()
density_max = particles.density.max()
temp_max += 10

t = 0
up = particles.identifier.unique()
q = particles

def set_graph_lim():
    temp.set_xlim(0,num_t)
    temp.set_ylim(0,temp_max)
    density.set_xlim(0,num_t)
    density.set_ylim(0,density_max)
set_graph_lim()

def updateBarGraphs(ps,i):
    p = ['particle 1']
    temp.bar(p,ps.iloc[i].temperature, width=0.5, label='Temp',
            color='#3CB371')
    density.bar(p,ps.iloc[i].density, width=0.5,  label='Density',
            color='#6495ED') #cornflowerblue
old = []
def updateFig(*args):
    global t, old

    temp.cla()
    density.cla()
    set_graph_lim()
    if (t == num_t):
        ax.cla()
        ax.set_xlim(0,nx); ax.set_ylim(0,ny); ax.set_zlim(0,nz)
        plot_image_lines()
        t = 0

    for id in up:
        p = particles[ (particles.identifier == id) &
                       (particles.timestep == t) ]
        p_less_than = particles[ (particles.identifier == id) &
                                 (particles.timestep <= t) ]
        p_density = p_less_than.density
        p_temp    = p_less_than.temperature
        p_time    = p_less_than.timestep

        temp.plot(p_time, p_temp, color='black')
        density.plot(p_time, p_density, color='black')
        if not p.empty:
            old = ax.scatter3D(p.x, p.y, p.z, marker='o', edgecolor='black',
                               color=cmap_c[t])
        else:
            temp.plot(p_time[-1:], p_temp[-1:], color='black', marker='x')
            density.plot(p_time[-1:], p_density[-1:], color='black',
                         marker='x')


    t += 1
    return old

        # temp.bar(p,ps.iloc[i].temperature, width=0.5, label='Temp',
            #          color='#3CB371')
            # density.bar(p,ps.iloc[i].density, width=0.5,  label='Density',
            #             color='#6495ED') #cornflowerblue

# --- animation ----
frame_delay_ms=100
ani = animation.FuncAnimation(fig, updateFig, interval=frame_delay_ms,
                              frames=num_t-1,repeat=True)

gif = True
gif = False
if (gif):
    ani.save('./test.gif', writer='imagemagick', fps=5)
else:
    plt.show()

print("Fin!")
