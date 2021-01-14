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

f_cray = open('thread_scaling_loops.txt')
# f_cheyenne = open('cheyenne_strong_scaling.txt')

# ---- read input data ----
plot_title = ""
# plot_title="Strong Scaling"

header = ['n_nodes', 'nx','nz','ny','np','x_images','y_images','n_particles',
          'timesteps', 'time',  'is_dry', 'wind_speed', 'scaling_run',
          'loop_type', 'n_threads']

df = pd.read_csv(f_cray, sep='\s+',header=None)
# df_c = pd.read_csv(f_cheyenne, sep='\s+',header=None, comment='#')

df.columns = header
# df_c.columns= header

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

def get_label(loop):
    label=''
    if (loop == 'do'):
        label = 'Do'
    elif (loop == 'dc'):
        label = 'Do Concurrent, with Threads'
    elif (loop == 'dc_no_threads'):
        label = 'Do Concurrent, without Threads'
    elif (loop == 'omp_simd'):
        label = 'OpenMP SIMD'
    elif (loop == 'omp_parallel_simd'):
        label = 'OpenMP Parallel SIMD'
    # elif (loop == ''):
    #     label = ''

    return label

def plot_data(data_in, name):
    for i,loop in enumerate(df.loop_type.unique()):
        label = get_label(loop)

        data = data_in[(data_in.loop_type == loop)]


        if (not data.empty):
            plt.plot(data.n_threads,   data.time, marker = '.', label=label)


plot_data(df, 'Cray')
# plot_data(df_c, 'SGI')
# plot_data(df_c, 'Cheyenne')




# plt.legend(title="Machine")
plt.legend()
plt.xlabel("OMP_NUM_THREADS value")
plt.ylabel("time (seconds)")
plt.title(plot_title)

# plt.yscale('log', base=2)
# plt.xscale('log', base=2)

# ax.set_xticklabels([])
# ax.set_yticklabels([])


# print(len(ax.get_yticklabels()))
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# sys.exit()


# filename="strong_scaling_"+str(graph_size)+"_loglog.png"
fig = plt.gcf()
fig.set_size_inches((4,3))
plt.tight_layout()
# plt.savefig(filename, dpi=300)
plt.show()
print("Fin!")
