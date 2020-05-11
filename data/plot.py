from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys


if (len(sys.argv) < 2):
    sys.exit("Error: too few arguments for `plot.py [graph_data.txt]`")


f = open(sys.argv[1])
l = f.readline()
xzy = [int(n) for n in l.split()]
x = xzy[0]; z = xzy[1]; y = xzy[2]

header = ['image','time step','exists','moved', 'x', 'y', 'z', 'u', 'v', 'w', \
          'density', 'temperature']

particles = pd.read_csv(f, sep='\s+',header=None, names=header)
p = particles

num_t = particles['time step'].max()
field = np.zeros((x,z,y,num_t))

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlim(0,x); ax.set_ylim(0,y); ax.set_zlim(0,z)

# ---- 3d plots ----
# ax.plot3D(particles['x'],particles['y'],particles['z'],
normalized = particles['temperature'] / particles['temperature'].max()
ax.scatter3D(particles.iloc[0]['x'], particles.iloc[0]['y'],
             particles.iloc[0]['z'], marker='X')
# ax.scatter3D(1,1,1, marker='X')

          # c=normalized.iloc[0])



t = 1
t_max = particles['time step'].max()

def updateFig(*args):
    global t
    if (t == t_max-5):
        anim_running = False
    x,y,z = particles.iloc[t][['x','y','z']]
    im = ax.scatter3D(x,y,z, marker='o')
    t += 1
    return im


# --- animation ----
frame_delay_ms=50
ani = animation.FuncAnimation(fig, updateFig, interval=frame_delay_ms,
                              frames=t_max-1,repeat=False)


plt.show()
print("Fin!")
