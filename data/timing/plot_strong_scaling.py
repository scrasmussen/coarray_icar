from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import matplotlib.scale as scale
import numpy as np
import pandas as pd
import sys

plt.rc('font', family='serif')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
# --- be sure to add ---
# AND fig = plt.figure(figsize=(4,3))
# AND plt.tight_layout() before plt.show()



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


# df = pd.read_csv(open('cheyenne_results.txt'), sep='\s+',header=None)
# df.columns = header4
# sys.exit()
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

# --- plot data ---
if (len(df.columns) != 12):
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
        for i,size in enumerate(df.nx.unique()):
            if (size == 500):
                marker = 'x'
            else:
                marker = 's'
            l = df.loc[df.nx == size, ['nx','ny','nz']].iloc[0]
            label = l.to_csv(header=False, index=False).replace('\n','x')[:-1]
            # label = label + ", optimization -" + str(c_flag)
            plt.plot(df[(df.nx == size) &
                    (df.compiler_flags == c_flag)].np,
                     df[(df.nx == size) &
                    (df.compiler_flags == c_flag)].time,
                     marker = marker,
                     label=label)



plt.legend(title="Dimension and Optimization")
plt.xlabel("number of images")
plt.ylabel("time (seconds)")
plt.title(plot_title)

plt.yscale('log', base=2)
plt.xscale('log', base=2)

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
    plt.tight_layout()
    plt.show()
print("Fin!")
