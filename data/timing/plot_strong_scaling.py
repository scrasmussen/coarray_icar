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
# plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.linewidth'] = 2

f_cray = open('cray_strong_scaling.txt')
f_cheyenne = open('cheyenne_strong_scaling.txt')

# ---- read input data ----
plot_title = ""
# plot_title="Strong Scaling"

header = ['n_nodes', 'nx','nz','ny','np','x_images','y_images','n_particles',
          'timesteps', 'time',  'is_dry', 'wind_speed', 'scaling_run']

df = pd.read_csv(f_cray, sep='\s+',header=None)
df_c = pd.read_csv(f_cheyenne, sep='\s+',header=None, comment='#')

df.columns = header
df_c.columns= header

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

# --- plot data ---
# old
# for i,size in enumerate(df.nx.unique()):
#     l = df.loc[df.nx == size, ['nx','ny','nz']].iloc[0]
#     label = l.to_csv(header=False, index=False).replace('\n','x')[:-1]
#     cray_label = "Cray " + label
#     cheyenne_label = "Cheyenne " + label
#     plt.plot(df[df.nx == size].np, df[df.nx == size].time, marker = '.',
#              label=cray_label)
#     plt.plot(df_c[df_c.nx == size].np, df_c[df_c.nx == size].time,
#              marker = 's', label=cheyenne_label)

graph_size=500
graph_size=2000


def plot_data(data_in, name):
    for i,size in enumerate(df.nx.unique()):
        if (size != graph_size):
            continue
        # label = (size**2 * 30) / 1000
        # label = name+' '+str(size)+'x'+str(size)+'x30'

        label = name

        data = data_in[(data_in.nx == size) & (data_in.n_particles == 0)]
        data_p = data_in[(data_in.nx == size) & (data_in.n_particles != 0)]

        if (not data.empty):
            plt.plot(data.np,   data.time, marker = '.', label=label)
        if (not data_p.empty):
            plt.plot(data_p.np, data_p.time, marker = 'x',
                     label=label+' w/ particles')

plot_data(df, 'Cray')
plot_data(df_c, 'SGI')
# plot_data(df_c, 'Cheyenne')




# plt.legend(title="Machine")
plt.legend()
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


filename="strong_scaling_"+str(graph_size)+"_loglog.png"
fig = plt.gcf()
fig.set_size_inches((4,3))
plt.tight_layout()
plt.savefig(filename, dpi=300)
plt.show()
print("Fin!")
