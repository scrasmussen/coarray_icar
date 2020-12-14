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
# --- be sure to add ---
# AND fig = plt.figure(figsize=(4,3))
# AND plt.tight_layout() before plt.show()

f_cray = open('cray_wind_results.txt')
# f_cheyenne = open('cheyenne_strong_scaling.txt')

# ---- read input data ----
plot_title = ""
if (len(sys.argv) > 1):
    if (sys.argv[1] == '-t' ):
        plot_title=sys.argv[2]
    else:
        plot_title="Cray Wind Scaling: 500x500x30"

plot_title="Cray: Particles Communicated"
plot_title=""


header = ['n_nodes','nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time', 'dry', 'wind_speed', 'num_communicated']

# df = pd.read_csv(open('cheyenne_results.txt'), sep='\s+',header=None)
# df.columns = header4
# sys.exit()
# df = pd.read_csv(f, sep='\s+',header=None, names=header)
df = pd.read_csv(f_cray, sep='\s+',header=None, comment='#')
df.columns = header
# df_c = pd.read_csv(f_cheyenne, sep='\s+',header=None, comment='#')

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

title_proxy = mpatches.Rectangle((0,0),0,0,color='w')
sat   = lines.Line2D([],[],marker='s',color='black')
unsat = lines.Line2D([],[],marker='.',color='black')
handles=[title_proxy]
labels=[r"\textbf{Number of Parcels}"]

df.num_communicated /= 1000000
df = df[df.n_nodes == 1]
# --- plot data ---
for i,size in enumerate(df.n_particles.unique()):

    d = df.loc[(df.n_particles == size) & (df.dry == 'T') ]
    p = plt.plot(d.wind_speed, d.num_communicated, marker = '.')

    d = df[(df.n_particles == size) & (df.dry == 'F') ]
    plt.plot(d.wind_speed, d.num_communicated, marker = 's', color=p[0].get_color())

    # prepare label
    t = r"$\sim$"
    num = int(round(size*44, -3)/1000)
    total_p = round(size * 44,-3)
    if (total_p / 1000000.0 < 1.0):
        num = str(int(round(total_p,-3)/1000)) + "k"
    else:
        num = str(int(round(total_p,-3)/1000000)) + " mil"
    label = t + num + " particles"
    patch = mpatches.Patch(color=p[0].get_color())
    labels.append(label)
    handles.append(patch)



plt.legend(handles,labels) #,title="Dimension and Optimization")
plt.xlabel("wind speed (meters per second)")
plt.ylabel("particles communicated (million)")
plt.title(plot_title)

# plt.yscale('log', base=2)
# plt.xscale('log', base=2)

# ax.set_xticklabels([])
# ax.set_yticklabels([])


# print(len(ax.get_yticklabels()))
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# sys.exit()
filename="particles_communicate_delta_wind.png"
fig = plt.gcf()
fig.set_size_inches((4,3))
plt.tight_layout()
plt.savefig(filename, dpi=300)
plt.show()
print("Fin!")
