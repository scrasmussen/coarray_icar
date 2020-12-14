import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas as pd
import sys

plt.rc('font', family='serif')
plt.rc('text',usetex=True)
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
# --- be sure to add ---
# AND plt.tight_layout() before plt.show()


data_path=sys.argv[1]

df = pd.DataFrame(columns=['id','num_comm','num_p'])
f_size = ['250k','1mil','5mil']
h_max = 0
for size in f_size:
    f_name = data_path + '/particles-num-comm-'+size+'-8ws.txt'
    f = open(f_name)
    df_new = pd.read_csv(f, sep='\s+',header=None, names=['id','num_comm'])
    df_new.insert(2,'num_p',size)
    df = pd.concat([df,df_new])

    l_max = df_new.num_comm.max()
    if (l_max > h_max):
        h_max = l_max
    del(df_new)



# f_cheyenne = open('cheyenne_strong_scaling.txt')
# for size in df.num_p.unique():
#     plt.hist(df[df.num_p == size].num_comm)
hist_data = [df[df.num_p == '250k'].num_comm,
             df[df.num_p == '1mil'].num_comm,
             df[df.num_p == '5mil'].num_comm]
hist_labels = ['250k','1 mil.','5 mil.']
h_bins = list(range(0,h_max))
plt.hist(hist_data, density=True,label=hist_labels, align='left', bins=h_bins)
# plt.hist(hist_data, label=hist_labels)




# ---- read input data ----

plot_title="Cray: Number of Times Particles Communicated"
plot_title = ""
plt.legend(title="Number of particles")
plt.xlabel("times communicated")
plt.ylabel("probability distribution")
plt.title(plot_title)



filename="particles_communicated_histogram_delta_particle_count.png"
fig = plt.gcf()
fig.set_size_inches((4,3))
plt.tight_layout()
plt.savefig(filename, dpi=300)
plt.show()
print("Fin!")
