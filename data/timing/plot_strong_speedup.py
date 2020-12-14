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


# ---- read input data ----
header = ['n_nodes', 'nx','nz','ny','np','x_images','y_images','n_particles',
          'timesteps', 'time',  'is_dry', 'wind_speed', 'scaling_run']

f_cray     = open('cray_strong_scaling.txt')
# f_cray     = open('cray_weak_scaling.txt') does this even make sense?
df   = pd.read_csv(f_cray, sep='\s+',header=None, comment='#')
df.columns = header

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

# --- plot data ---
for i,nx in enumerate(df.nx.unique()):
    print(nx)
    if (nx == 500):
        label = '500x500'
        color = 'r'
    else:
        label = '2000x2000'
        color = 'b'

    data = df[(df.nx == nx) & (df.n_particles == 0)]
    data_p = df[(df.nx == nx) & (df.n_particles != 0)]
    t1 = (data.loc[data.np == 1].time.values[0])
    t1_p = (data_p.loc[data_p.np == 1].time.values[0])

    # speedup = t_1 / t_n
    speedup = t1 / data.time
    speedup_p = t1_p / data_p.time

    data = df[(df.nx == nx) &
              (df.n_particles == 0)]
    data_p = df[(df.nx == nx) &
                (df.n_particles != 0)]

    plt.plot(data.np.values, speedup, marker = '.',
             label=label, color=color)
    plt.plot(data_p.np, speedup_p, marker = 'x',
             label=str(label)+' w/ particles', color=color)
             # color=discrete_cmap(i*4))

# --- plot ideal speedup ---
max_np = df.np.max()
plt.plot([1,max_np],[1,max_np],label="ideal speedup", linestyle='--', color='g')
# --- plot strong scaling speedup ---
amdahl_x=[]
amdahl=[]
if (f_cray.name == 'timing_results.txt'):
    s = 0.1
elif (f_cray.name == 'cray_strong_scaling.txt'):
    s = 0.03
else:
    s = 0.1
p = 1 - s
for n in range(1,max_np+1):
    amdahl_x.append(n)
    amdahl.append(1 / (s + p / n))
# plt.plot(amdahl_x,amdahl,label="Amdahl's Law s="+str(s), linestyle=':')


plt.legend(title="Dimensions")
plt.xlabel("number of images")
plt.ylabel("speedup")
# plt.title("Speedup")

# plt.xscale('log', base=2)
# plt.yscale('log', base=2)
# ax.set_xticklabels([])
# ax.set_yticklabels([])


# print(len(ax.get_yticklabels()))
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# sys.exit()

filename="strong_scaling_speedup.png"
fig = plt.gcf()
fig.set_size_inches((4,3))
plt.tight_layout()
plt.savefig(filename, dpi=300)
plt.show()
print("Fin!")
