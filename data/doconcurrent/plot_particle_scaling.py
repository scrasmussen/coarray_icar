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
f_system76_g = open('system76_particle_scaling_gfortran.txt')
f_cray = open('cray_particle_scaling2.txt')

# ---- read input data ----
plot_title = ""
# plot_title="Strong Scaling"

header = ['n_nodes', 'nx','nz','ny','np','x_images','y_images','n_particles',
          'timesteps', 'time',  'is_dry', 'wind_speed', 'scaling_run',
          'loop_type', 'n_threads']

df_c = pd.read_csv(f_cray, sep='\s+',header=None, comment='#')
df_c.columns = header

df_s = pd.read_csv(f_system76, sep='\s+',header=None, comment='#')
df_s.columns = header

df_s_g = pd.read_csv(f_system76_g, sep='\s+',header=None, comment='#')
df_s_g.columns = header




# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

def get_label(loop):
    # print(loop)
    label=''
    # Cray loops, cleaner
    if (loop == 'do-func'):
        label = 'Do w/ functions'
    elif (loop == 'do-nfunc'):
        label = 'Do no functions'
    elif (loop == 'do-local-funcs'):
        label = 'Do w/ local functions'
    elif (loop == 'do-all-minus-dry-interop'):
        label = 'Do, local functions, one non-module one'
    elif (loop == 'do-all-minus-exner-satmr'):
        label = 'Do, no functions but two non-module ones'

    elif (loop == 'dc-func'):
        label = 'Do Concurrent w/ functions'
    elif (loop == 'dc-nfunc'):
        label = 'Do Concurrent no functions'
    elif (loop == 'dc-threads-func'):
        label = 'Do Concurrent, threads & functions'
    elif (loop == 'dc-threads-nfunc'):
        label = 'Do Concurrent, threads & no functions'

    elif (loop == 'omp-simd-func'):
        label = 'OpenMP SIMD w/ functions'
    elif (loop == 'omp-simd-nfunc'):
        label = 'OpenMP SIMD no functions'
    elif (loop == 'omp-parallel-simd-func'):
        label = 'OpenMP Parallel SIMD w/ functions'
    elif (loop == 'omp-parallel-simd-nfunc'):
        label = 'OpenMP Parallel SIMD no functions'



    # NVidia
    # elif (loop == 'omp_parallel_do_simd'):
    #     label = 'OpenMP Parallel SIMD???'
    # elif (loop == 'omp_do_simd'):
    #     label = 'OpenMP SIMD???'
    elif (loop == 'do'): # ??
        label = 'Do w/ functions???'
    elif (loop == 'omp_do_simd'): # ??
        label = 'OpenMP Do SIMD w/ functions???'
    elif (loop == 'omp_parallel_do_simd'): # ??
        label = 'OpenMP Parallel Do SIMD w/ functions???'

    elif (loop == 'dc_gpu'):
        label = 'Do Concurrent, GPU'
    elif (loop == 'dc_multicore'):
        label = 'Do Concurrent, Multicore'
    elif (loop == 'omp_do_simd_gpu'):
        label = 'OpenMP SIMD, GPU w/ functions'
    elif (loop == 'omp_do_simd_nfunc_gpu'):
        label = 'OpenMP SIMD, GPU no functions'
    elif (loop == 'omp_do_simd_multicore'):
        label = 'OpenMP SIMD, Multicore w/ functions'
    elif (loop == 'omp_do_simd_nfunc_multicore'):
        label = 'OpenMP SIMD, Multicore no functions'
    elif (loop == 'omp_parallel_do_simd_gpu'):
        label = 'OpenMP Parallel SIMD, GPU w/ functions'
    elif (loop == 'omp_parallel_do_simd_nfunc_gpu'):
        label = 'OpenMP Parallel SIMD, GPU no functions'
    elif (loop == 'omp_parallel_do_simd_multicore'):
        label = 'OpenMP Parallel SIMD, Multicore w/ functions'
    elif (loop == 'omp_parallel_do_simd_nfunc_multicore'):
        label = 'OpenMP Parallel SIMD, Multicore no functions'
    elif (loop == 'do_gpu_nfunc'):
        label = 'Do GPU no functions'
    elif (loop == 'do_multicore_nfunc'):
        label = 'Do Multicore no functions'
    elif (loop == 'do_multicore'):
        label = 'Do Multicore w/ Functions'

    else:
        print("WARNING, LABEL NOT FOUND FOR |"+loop+"|")
        sys.exit()

    return label



def plot_data(data_in, name):
    unique_loop_list = \
        data_in[(data_in.n_particles==1000000) &
                # (data_in.time < 100000) ].loop_type.unique()
                (data_in.time < 100) ].loop_type.unique()


                # (data_in.time < 100) & # System76
                # (data_in.time < 500) ].loop_type.unique() # System76
                # (data_in.time < 100) ].loop_type.unique() # Cray
    for i,loop_type in enumerate(unique_loop_list):
        label = get_label(loop_type)
        if '???' in label:
            continue

        # if 'GPU' not in label:
        #     continue


        data = data_in[(data_in.loop_type == loop_type)]



        # data = data_in[(data_in.loop_type == loop_type) &
        #                (data_in.time<600) &
        #                # (data_in.time<100) &
        #                (data_in.n_particles == 1000000)]

        # data = data_in[(data_in.loop_type == loop_type) &
        # (data_in[data_in['n_particles'] >= 1000000].time<500)]
        marker = '.'
        if 'GPU' in label:
            marker = 's'
        elif 'OpenMP' in label:
            marker = 'x'
        elif 'Do Concurrent' in label:
            marker = 'o'


        if (not data.empty):
            data = data_in[(data_in.loop_type == loop_type)]
            plt.plot(data.n_particles, data.time, marker = marker, label=label)


# plot_data(df_c, 'Cray')
# plot_data(df_g, 'SGI')
# plot_data(df_s, 'System76')
plot_data(df_s_g, 'System76 Gfortran')
# plot_data(df_c, 'Cray')




# plt.legend(title="Machine")
plt.legend()
plt.xlabel("Number of Particles")
plt.ylabel("time (seconds)")
plt.title(plot_title)

# plt.yscale('log', base=2)
plt.xscale('log', base=10)

# ax.set_xticklabels([])
# ax.set_yticklabels([])

bottom, upper = plt.ylim()
# plt.ylim(0,600)
# plt.ylim(30,upper)
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
