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

header = ['image','time step','exists','moved', 'x', 'y', 'z', 'u', 'v', 'w', \
          'density', 'temperature']
particles = pd.read_csv(f, sep='\s+',header=None, names=header)


num_t = particles['time step'].max()
field = np.zeros((nx,nz,ny,num_t))

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlim(0,nx); ax.set_ylim(0,ny); ax.set_zlim(0,nz)



print("ARTLESS: remove this once real data is produced")
new_temp = [i+x for i,x in enumerate(particles['temperature'])]
particles['temperature'] = new_temp

print("ARTLESS: fix normalization")
norm = matplotlib.colors.Normalize(vmin=particles['temperature'].min()-10,
                                   vmax=particles['temperature'].max())




# ---- Set up colors ----
cmap = plt.cm.Greens
cmap_c = cmap(norm(particles.temperature.values))

# ---- 3d plots ----
# ax.plot3D(particles['x'],particles['y'],particles['z'],
# normalized = particles['temperature'] / particles['temperature'].max()
old = ax.scatter3D(particles.iloc[0]['x'], particles.iloc[0]['y'],
             particles.iloc[0]['z'], marker='X', color=cmap_c[0])


t = 1
t_max = particles['time step'].max()

def updateFig(*args):
    global t, old
    old.remove()
    x,y,z = particles.iloc[t-1][['x','y','z']]
    im = ax.scatter3D(x,y,z, marker='.', color=cmap_c[t-1])
    x,y,z = particles.iloc[t][['x','y','z']]
    old = ax.scatter3D(x,y,z, marker='o', color=cmap_c[t])

    t += 1
    return im


# --- animation ----
frame_delay_ms=100
ani = animation.FuncAnimation(fig, updateFig, interval=frame_delay_ms,
                              frames=t_max-1,repeat=False)

plt.show()
print("Fin!")
