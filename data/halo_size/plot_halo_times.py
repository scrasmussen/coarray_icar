from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys

plt.rc('font', family='serif')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2

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

rows = 1 # 2
cols = 1
fig = plt.figure()



# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')
# --- plot data ---
df = df[df.timesteps == 200]
# ax = fig.add_subplot(111)
ax1 = fig.add_subplot(rows,cols,1)
# ax2 = fig.add_subplot(rows,cols,2)

def plot_size(nx, ax_in):
    t_ax = ax_in
    if (nx == 2000):
        color = 'blue'
    else:
        color = 'red'
    for i,nodes in enumerate(df.n_nodes.unique()):
        label = str(nodes) + ' node'
        if (nodes != 1):
            label = label + 's'
        if (nodes == 1):
            marker = '.'
        else:
            marker = 'x'

        label = str(nx)+' problem size, '+label
        t_ax.plot(df[(df.n_nodes == nodes) &
                    (df.nx == nx)].halo_depth,
                 df[(df.n_nodes == nodes) &
                    (df.nx == nx)].time,
                 marker=marker,
                 label=label, color=color)
        t_ax.set_xlabel("halo depth")
        t_ax.set_ylabel("time (seconds)")

        # color=discrete_cmap(i*4))

ax1.set_ylim([0,df.time.max()+10])
plot_size(2000, ax1)
plot_size(500,  ax1)



plt.legend(title="Number of nodes, problem size")
# plt.title(plot_title)

# plt.xscale('log', basex=2)
# ax.set_xticklabels([])
# ax.set_yticklabels([])


# print(len(ax.get_yticklabels()))
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# sys.exit()
# --- be sure to add ---
# AND fig = plt.figure(figsize=(4,3))
# AND plt.tight_layout() before plt.show()
plt.tight_layout()
plt.show()
print("Fin!")
