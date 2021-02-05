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

# f_cray = open('cray_particle_scaling.txt')
f_system76 = open('system76_particle_scaling.txt')
f_system76 = open('cray_particle_scaling.txt')

# ---- read input data ----
plot_title = ""
# plot_title="Strong Scaling"

header = ['n_nodes', 'nx','nz','ny','np','x_images','y_images','n_particles',
          'timesteps', 'time',  'is_dry', 'wind_speed', 'scaling_run',
          'loop_type', 'n_threads']

# df = pd.read_csv(f_cray, sep='\s+',header=None)
df_s = pd.read_csv(f_system76, sep='\s+',header=None, comment='#')

# df.columns = header
df_s.columns= header

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

def get_label(loop):
    # print(":"+loop+":")
    # sys.exit()
    print(loop)
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
    # NVidia
    elif (loop == 'omp_parallel_do_simd'):
        label = 'OpenMP Parallel SIMD???'
    elif (loop == 'omp_do_simd'):
        label = 'OpenMP SIMD???'
    elif (loop == 'dc_gpu'):
        label = 'Do Concurrent, GPU'
    elif (loop == 'dc_multicore'):
        label = 'Do Concurrent, Multicore'
    elif (loop == 'omp_do_simd_gpu'):
        label = 'OpenMP SIMD, GPU w/ Procedures'
    elif (loop == 'omp_do_simd_expand_gpu'):
        label = 'OpenMP SIMD, GPU'
    elif (loop == 'omp_do_simd_multicore'):
        label = 'OpenMP SIMD, Multicore w/ Procedures'
    elif (loop == 'omp_do_simd_expand_multicore'):
        label = 'OpenMP SIMD, Multicore'
    elif (loop == 'omp_parallel_do_simd_gpu'):
        label = 'OpenMP Parallel SIMD, GPU w/ Procedures'
    elif (loop == 'omp_parallel_do_simd_expand_gpu'):
        label = 'OpenMP Parallel SIMD, GPU'
    elif (loop == 'omp_parallel_do_simd_multicore'):
        label = 'OpenMP parallel_SIMD, Multicore w/ Procedures'
    elif (loop == 'omp_parallel_do_simd_expand_multicore'):
        label = 'OpenMP parallel_SIMD, Multicore '
    elif (loop == 'do_gpu_expanded'):
        label = 'Do GPU'
    elif (loop == 'do_multicore_expanded'):
        label = 'Do Multicore'
    elif (loop == 'do_multicore'):
        label = 'Do Multicore w/ Procedures'

    else:
        print("WARNING, LABEL NOT FOUND FOR "+loop)
        sys.exit()

    return label

def plot_data(data_in, name):
    unique_loop_list = \
        df_s[(df_s.n_particles==1000000) & (df_s.time < 1000) ].loop_type.unique()
    for i,loop_type in enumerate(unique_loop_list):
        label = get_label(loop_type)

        # data = data_in[(data_in.loop_type == loop_type) & (data_in[].time<500)]
        data = data_in[(data_in.loop_type == loop_type)]


        if (not data.empty):
            plt.plot(data.n_particles, data.time, marker = '.', label=label)


# plot_data(df, 'Cray')
# plot_data(df_g, 'SGI')
plot_data(df_s, 'System76')




# plt.legend(title="Machine")
plt.legend()
plt.xlabel("Number of Particles")
plt.ylabel("time (seconds)")
plt.title(plot_title)

# plt.yscale('log', base=2)
# plt.xscale('log', base=2)

# ax.set_xticklabels([])
# ax.set_yticklabels([])

# plt.ylim(0,600)
# print(len(ax.get_yticklabels()))
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# sys.exit()


# filename="strong_scaling_"+str(graph_size)+"_loglog.png"
fig = plt.gcf()
# fig.set_size_inches((4,3))
# plt.tight_layout()
# plt.savefig(filename, dpi=300)
plt.show()
print("Fin!")
