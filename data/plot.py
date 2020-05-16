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

header = ['image','time step','exists','moved', 'x', 'y', 'z', 'u', 'v', 'w', \
          'density', 'temperature']
particles = pd.read_csv(f, sep='\s+',header=None, names=header)

num_t = particles['time step'].max()
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




# fig = plt.figure()
# ax = plt.axes(projection='3d')
# fig, (ax,ax2) = plt.subplots(2)

fig = plt.figure(figsize=plt.figaspect(0.5))

# ax = fig.add_subplot(2,1,1,projection='3d')  # THIS WORKS!!!!!!!!
ax = fig.add_subplot(1,2,1,projection='3d')
# ax.plot3D()
ax.set_xlim(0,nx); ax.set_ylim(0,ny); ax.set_zlim(0,nz)
ax.set_xlabel="x axis"
ax.set_ylabel="y axis"
ax.set_zlabel="z axis"




# ---- 3d plots ----
# ax.plot3D(particles['x'],particles['y'],particles['z'],
# normalized = particles['temperature'] / particles['temperature'].max()

# --- plot image lines ---
def plotlines():
    for i in range (1,ximages):
        ax.plot(xs=[(i*nx)/ximages,(i*nx)/ximages], ys=[0,ny], color='black')
    for i in range (1,yimages):
        ax.plot(xs=[0,nx], ys=[(i*ny)/yimages,(i*ny)/yimages], color='black')

# ---- 3d plots ----
plotlines()
t = 0




temp = fig.add_subplot(2,2,2)
density = fig.add_subplot(2,2,4)

temp_min = particles.temperature.min()
temp_max = particles.temperature.max()
density_min = particles.density.min()
density_max = particles.density.max()

if (temp_min > temp_max-10):
    temp_min -= 5
    temp_max += 5
if (density_min > density_max-10):
    density_min -= 5
    density_max += 5
temp_min = 0
density_min = 0

temp.set_ylim(temp_min, temp_max)
density.set_ylim(density_min-10, density_max)

def updateBarGraphs(ps,i):
    p = ['particle 1']
    temp.bar(p,ps.iloc[i].temperature, width=0.5, label='Temp',
            color='#3CB371')
    density.bar(p,ps.iloc[i].density, width=0.5,  label='Density',
            color='#6495ED') #cornflowerblue

def updateFig(*args):
    global t, old
    if (t == num_t):
        plt.cla()
        temp.cla()
        density.cla()
        temp.set_ylim(temp_min, temp_max)
        density.set_ylim(density_min-10, density_max)

        ax.set_xlim(0,nx); ax.set_ylim(0,ny); ax.set_zlim(0,nz)
        plotlines()
        t = 0
    if (t != 0):
        old.remove()
        x,y,z = particles.iloc[t-1][['x','y','z']]
        im = ax.scatter3D(x,y,z, marker='.', color=cmap_c[t-1])

    updateBarGraphs(particles,t)
    x,y,z = particles.iloc[t][['x','y','z']]
    old = ax.scatter3D(x,y,z, marker='o', color=cmap_c[t])

    t += 1
    return old


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
