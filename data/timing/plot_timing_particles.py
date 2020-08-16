from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys


if (len(sys.argv) < 2):
    sys.exit("Error: too few arguments for `plot.py [graph_data.txt]`")

# ---- read input data ----
f = open(sys.argv[1])
header = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time']
df = pd.read_csv(f, sep='\s+',header=None, names=header)


# --- order data ---
df.sort_values('n_particles', inplace=True)

# --- animation ----
gif    = True
gif    = False
repeat = False if gif else True

# ax.scatter(particles.timestep, df.time,
#            cmap=discrete_cmap, c=particles.identifier)
#            marker='.', edgecolor='black')
#            color=cmap_c[0])

# --- setup colormap ---
discrete_cmap = plt.get_cmap('tab20b')

# --- plot data ---
for i,size in enumerate(df.nx.unique()):
    l = df.loc[df.nx == size, ['nx','ny','nz']].iloc[0]
    label = l.to_csv(header=False, index=False).replace('\n','x')[:-1]
    label = "With particles"
    plt.plot(df.n_particles, df.time, marker = '.',
             label=label)
             # color=discrete_cmap(i*4))

plt.plot([0,df.n_particles.max()], [943.596,943.596], marker = '.',
             label="Baseline: particles turned off", color='red')


# plt.legend(title="Dimensions")
plt.legend()
plt.xlabel("number of particles (million)")
plt.ylabel("time (seconds)")

def format_xaxis(value, tick_num):
    return int(value / 1000000.0)

ax = plt.gca()
ax.ticklabel_format(axis="x",style='plain')
ax.xaxis.set_major_formatter(plt.FuncFormatter(format_xaxis))
# plt.title("Particle Scaling, 500 timesteps, 16 images")

model = LinearRegression()
model.fit(df.n_particles.values.reshape(-1,1), df.time.values.reshape(-1,1))
print("Linear Regression: y = " +str(model.coef_[0][0])+"x + " +
      str(model.intercept_[0]))
r_squared = model.score(df.n_particles.values.reshape(-1,1),
                      df.time.values.reshape(-1,1))
print("Coefficient of determination R^2 = ", str(r_squared))
# plt.plot(df.n_particles, model.predict(df.n_particles.values.reshape(-1,1)), color='red')


a = model.coef_[0][0]
b = model.intercept_[0]
y = a * 40000000 + b
time = 12 * 60 * 60
print("y = ", str(y), "data = ", time, "error =",(time-y)/y)

# plt.xscale('log', basex=2)
# ax.set_xticklabels([])
# ax.set_yticklabels([])


# print(len(ax.get_yticklabels()))
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# sys.exit()

if (gif):
    print("Gif!")
    # ani.save('test.gif', writer=animation.PillowWriter, fps=None,dpi=20) # fps was 5
    # ani.save('test.gif', writer=animation.ImageMagickWriter, fps=None) # fps was 5
else:
    plt.show()
print("Fin!")
