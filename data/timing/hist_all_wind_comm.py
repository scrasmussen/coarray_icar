import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
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

df = pd.DataFrame(columns=['num_comm'])
f_wind = [1,2,4,8,16,32,64]
max=0
# f_wind = [1]
for wind in f_wind:
    f_name = data_path + '/particles-num-comm-250k-' + str(wind) + 'ws.txt'
    f = open(f_name)
    df_new = pd.read_csv(f, sep='\s+',header=None, names=['num_comm'])
    norm = df_new.value_counts(normalize=True, sort=False)
    plt.plot([x[0] for x in norm.index], norm.values*100, '-', label=str(wind))

    df_max = df_new.num_comm.max()
    if (df_max > max):
        max = df_max
    del(df_new)


ax = plt.gca()

ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(1))

ax.grid(b=True, which='major', linestyle='-')
# ax.grid(b=True, which='minor', linestyle='--')

plt.xlim(left=0)
plt.ylim(bottom=0)

# ax.minorticks_on()
# plt.minorticks_on()
# data=(ticks,ticks))
# markevery=(4))



# f_cheyenne = open('cheyenne_strong_scaling.txt')
# for size in df.num_p.unique():
#     plt.hist(df[df.num_p == size].num_comm)
# hist_data = [df[df.wind == 1].num_comm,
#              df[df.wind == 2].num_comm,
#              df[df.wind == 4].num_comm,
#              df[df.wind == 8].num_comm,
#              df[df.wind == 16].num_comm,
#              df[df.wind == 32].num_comm,
#              df[df.wind == 64].num_comm]
# hist_labels = f_wind
# hist_bins = list(range(0,19))
# plt.hist(hist_data, density=True,label=hist_labels, bins=hist_bins)
# plt.hist(hist_data, label=hist_labels)




# ---- read input data ----

plot_title = ""
plot_title="Cray: Percentage of Times Particles are Communicated"
plt.legend(title="Wind Speed")
plt.xlabel("number of times parcel communicated")
plt.ylabel("percentage of parcels")
plt.title(plot_title)

plt.tight_layout()
plt.show()
print("Fin!")
