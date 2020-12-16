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
plot_title = "" # "Weak Scaling"
header = ['n_nodes', 'nx','nz','ny','np','x_images','y_images','n_particles',
          'timesteps', 'time',  'is_dry', 'wind_speed', 'scaling_run']

f_cray     = open('cray_weak_scaling.txt')
f_cheyenne = open('cheyenne_weak_scaling.txt')
df   = pd.read_csv(f_cray, sep='\s+',header=None, comment='#')
df_c = pd.read_csv(f_cheyenne, sep='\s+',header=None, comment='#')
df.columns = header
df_c.columns = header

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

# --- plot data ---
def plot_data(data_in, name):
    global label_n
    for i,run in enumerate(df.scaling_run.unique()):
        if (run == 1):
            marker = 'x'
            # continue
            label_n = 20*20*30 / 1000
        else:
            marker = '.'
            continue
            label_n = 160*160*30 / 1000

        label = name + " " + str(int(label_n)) + "k "

        data = data_in[(data_in.scaling_run == run) &
                       (data_in.n_particles == 0)]
        data_p = data_in[(data_in.scaling_run == run) &
                         (data_in.n_particles != 0)]
        plt.plot(data.np,   data.time, marker = '.', label=label)
        plt.plot(data_p.np, data_p.time, marker = 'x',
                 label=label+' w/ particles')

plot_data(df, 'Cray')
plot_data(df_c, 'SGI')

ax = plt.gca()
ax.set_ylim(0.0)
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


filename="weak_scaling_cray_"+str(int(label_n))+"k.png"
fig = plt.gcf()
fig.set_size_inches((4,3))
plt.tight_layout()
plt.savefig(filename, dpi=300)
plt.show()
print("Fin!")
