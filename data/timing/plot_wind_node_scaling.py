from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import matplotlib.scale as scale
import matplotlib.lines as lines
import numpy as np
import pandas as pd
import sys

plt.rc('font', family='serif')
plt.rc('text',usetex=True)
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2

f_cray = open('cray_wind_results.txt')

plot_title = ""
# plot_title = "250k particles"
plot_title = "Varying Wind Speed and Nodes, Constant Particle Count"


header = ['n_nodes', 'nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time', 'dry', 'wind_speed', 'num_communicated']

df = pd.read_csv(f_cray, sep='\s+',header=None, comment='#')
df.columns = header

# --- plot data ---
for i,nodes in enumerate(df.n_nodes.unique()):
    size = 113636
    d = df.loc[(df.n_particles == size) & (df.dry == 'T')
    & (df.n_nodes == nodes) & (df.wind_speed <= 8.0)]

    # prepare label
    t = r"$\sim$"
    num = int(round(size*44, -3)/1000)
    total_p = round(size * 44,-3)
    num = str(int(round(total_p,-3)/1000)) + "k"

    label = str(nodes) + " node"
    label += 's x ' if (nodes > 1) else ' x '
    label += str(int(44 / nodes)) + ' processes'

    marker = '.' if (nodes > 1) else 's'
    plt.plot(d.wind_speed, d.time, marker = marker, label=label)


plt.legend()
plt.xlabel("wind speed (meters per second)")
plt.ylabel("time (seconds)")
plt.title(plot_title)


plt.tight_layout()
plt.show()
print("Fin!")
