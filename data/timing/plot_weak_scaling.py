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
        plot_title="Strong Scaling"


header = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time']
header2 = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time', 'n_nodes']
header3 = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time',  'halo_size', 'n_nodes']
header4 = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time',  'halo_depth', 'n_nodes', 'compiler_flags']
header5 = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time',  'halo_depth', 'n_nodes', 'compiler_flags', 'scaling_run']

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
elif (len(df.columns) == 13):
    df.columns = header5

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
if (len(df.columns) != 13):
    for i,size in enumerate(df.nx.unique()):
        l = df.loc[df.nx == size, ['nx','ny','nz']].iloc[0]
        label = l.to_csv(header=False, index=False).replace('\n','x')[:-1]
        plt.plot(df[df.nx == size].np, df[df.nx == size].time, marker = '.',
                 label=label)
             # color=discrete_cmap(i*4))
else:
    c_flags = ['O0','O3']
    c_flags = ['O3']
    for c_flag in c_flags:
        for i,run in enumerate(df.scaling_run.unique()):
            print(run)
            if (run == 1):
                marker = 'x'
                label = 20*20*30 / 1000
            else:
                marker = 's'
                label = 160*160*30 / 1000
            label = str(int(label)) + "k"
            plt.plot(df[(df.scaling_run == run)].np,
                     df[(df.scaling_run == run)].time,
                     marker = marker,
                     label=label)



plt.legend(title="Problem size per image")
plt.xlabel("number of images")
plt.ylabel("time (seconds)")
plt.title(plot_title)

# plt.xscale('log', base=2)
# plt.yscale('log', base=2)
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
