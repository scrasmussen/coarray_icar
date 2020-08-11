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

# ---- read input data ----
f = open(sys.argv[1])
header = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time']
df = pd.read_csv(f, sep='\s+',header=None, names=header)


# --- animation ----
gif    = True
gif    = False
repeat = False if gif else True

# ax.scatter(particles.timestep, df.time,
#            cmap=discrete_cmap, c=particles.identifier)
#            marker='.', edgecolor='black')
#            color=cmap_c[0])

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

# --- plot data ---
for i,size in enumerate(df.nx.unique()):
    l = df.loc[df.nx == size, ['nx','ny','nz']].iloc[0]
    label = l.to_csv(header=False, index=False).replace('\n','x')[:-1]

    data = df[df.nx == size]
    t1 = (data.loc[data.np == 1].time.values[0])

    # speedup = t_1 / t_n
    speedup = t1 / data.time

    plt.plot(df[df.nx == size].np, speedup, marker = '.',
             label=label)
             # color=discrete_cmap(i*4))

# --- plot ideal speedup ---
max_np = df.np.max()
plt.plot([1,max_np],[1,max_np],label="ideal speedup", linestyle='--')
# --- plot strong scaling speedup ---
amdahl_x=[]
amdahl=[]
s = 0.1
p = 1 - s
for n in range(1,max_np+1):
    amdahl_x.append(n)
    amdahl.append(1 / (s + p / n))
plt.plot(amdahl_x,amdahl,label="Amdahl's Law s="+str(s), linestyle=':')


plt.legend(title="Dimensions")
plt.xlabel("number of images")
plt.ylabel("speedup")
plt.title("Speedup")

# plt.xscale('log', basex=2)
# ax.set_xticklabels([])
# ax.set_yticklabels([])


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
