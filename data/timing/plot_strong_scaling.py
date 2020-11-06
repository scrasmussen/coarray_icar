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

f_cray = open('cray_strong_scaling.txt')
f_cheyenne = open('cheyenne_strong_scaling.txt')

# ---- read input data ----
plot_title = ""
if (len(sys.argv) > 1):
    if (sys.argv[1] == '-t' ):
        plot_title=sys.argv[2]
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


# df = pd.read_csv(open('cheyenne_results.txt'), sep='\s+',header=None)
# df.columns = header4
# sys.exit()
# df = pd.read_csv(f, sep='\s+',header=None, names=header)
df = pd.read_csv(f_cray, sep='\s+',header=None)
df_c = pd.read_csv(f_cheyenne, sep='\s+',header=None)
if (len(df.columns) == 9):
    df.columns = header
elif (len(df.columns) == 10):
    df.columns = header2
elif (len(df.columns) == 11):
    df.columns = header3
elif (len(df.columns) == 12):
    df.columns = header4
    df_c.columns= header4
elif (len(df.columns) == 13):
    df.columns = header5
    df_c.columns= header5

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

# --- plot data ---
for i,size in enumerate(df.nx.unique()):
    l = df.loc[df.nx == size, ['nx','ny','nz']].iloc[0]
    label = l.to_csv(header=False, index=False).replace('\n','x')[:-1]
    plt.plot(df[df.nx == size].np, df[df.nx == size].time, marker = '.',
             label=label)
    plt.plot(df_c[df_c.nx == size].np, df_c[df_c.nx == size].time,
             marker = 's', label=label)

# if (len(df.columns) != 12):
#     for i,size in enumerate(df.nx.unique()):
#         l = df.loc[df.nx == size, ['nx','ny','nz']].iloc[0]
#         label = l.to_csv(header=False, index=False).replace('\n','x')[:-1]
#         plt.plot(df[df.nx == size].np, df[df.nx == size].time, marker = '.',
#                  label=label)
#         plt.plot(df_c[df_c.nx == size].np, df_c[df_c.nx == size].time,
#                  marker = 'o', label=label)

#              # color=discrete_cmap(i*4))

# else: # not being used right now
#     c_flags = ['O0','O3']
#     c_flags = ['O3']
#     for c_flag in c_flags:
#         for i,size in enumerate(df.nx.unique()):
#             if (size == 500):
#                 marker = 'x'
#             else:
#                 marker = 's'
#             l = df.loc[df.nx == size, ['nx','ny','nz']].iloc[0]
#             label = l.to_csv(header=False, index=False).replace('\n','x')[:-1]
#             # label = label + ", optimization -" + str(c_flag)
#             plt.plot(df[(df.nx == size) &
#                     (df.compiler_flags == c_flag)].np,
#                      df[(df.nx == size) &
#                     (df.compiler_flags == c_flag)].time,
#                      marker = marker,
#                      label=label)





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

plt.tight_layout()
plt.show()
print("Fin!")
