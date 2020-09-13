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
plot_title = ""
if (len(sys.argv) > 2):
    if (sys.argv[2] == '-t' ):
        plot_title=sys.argv[3]
    else:
        plot_title="Halo size communication time"


header = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time']
header2 = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time', 'n_nodes']
header3 = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time',  'halo_depth', 'n_nodes']
header4 = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time',  'halo_depth', 'n_nodes', 'compiler_flags']

# df = pd.read_csv(f, sep='\s+',header=None, names=header)
df = pd.read_csv(f, sep='\s+',header=None)
if (len(df.columns) == 9):
    df.columns = header
elif (len(df.columns) == 10):
    df.columns = header2
elif (len(df.columns) == 11):
    df.columns = header3
elif (len(df.columns) == 12):
    df.columns = header4

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
c_flags = [0,3]
# --- plot data ---
for c_flag in c_flags:
    for i,nodes in enumerate(df.n_nodes.unique()):
        if (nodes == 1):
            label = str(nodes) + " node, optimization -O" + str(c_flag)
        else:
            label = str(nodes) + " nodes, optimization -O" + str(c_flag)
        plt.plot(df[(df.n_nodes == nodes) &
                    (df.compiler_flags == c_flag)].halo_depth,
                 df[(df.n_nodes == nodes) &
                    (df.compiler_flags == c_flag)].time,
                 marker = '.',
                 label=label)
        # color=discrete_cmap(i*4))

plt.legend(title="Number of nodes")
plt.xlabel("halo depth")
plt.ylabel("time (seconds)")
plt.title(plot_title)

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
